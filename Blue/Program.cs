using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using WorkingDogsCore;

namespace Blue
{
    // Blue corrects a set of reads using a k-mer consensus table derived from a set of reads (a .cbt file created by Tessel).
    // 
    // usage: Blue <options> <k-mer file> <list of file names or patterns for reads files to be corrected>
    //                  [-r run] : tag inserted into corrected reads files, and name used for stats file. Default is "corrected"
    //                   -m minReps : min k-mer depth for OK when scanning reads (used only when dynamic calculation fails)
    //                  [-f file format] (fasta or fastq) : Blue will try to work the format out from the file so this parameter is optional
    //                  [-hp] : do additional checks for possible errors at the end of every homopolymer string - intended for 454 and IonTorrent
    //                  [-t #threads] : no. of threads to use for the correction. Default is 1 thread.
    //                  [-l length] : trim reads to this length before correction
    //                  [-fixed] : fixed length reads, corrected reads are extended or trimmed to be the same length as raw reads (default)
    //                  [-variable] : variable length reads, corrected reads may be trimmed or extended 
    //                  [-extend maxbases] : extend trimmed reads by up to maxbases bases.
    //                  [-good nn%] : discard reads that are not at least this long after correction/trimming. Default is 70%.
    //                  [-paired] : handle pairs of read files as paired (both need to be good/corrected)
    //                  [-unpaired] : don't treat pairs of files as pair sets
    //                  [-o outputDir] : directory for corrected reads etc  
    //
    // The first filename-like parameter (not preceded by a '-') is the k-mer consensus table (.cbt file).
    // This is followed by a list of file names or patterns for the sequence data files to be corrected. 
    // Blue accepts either fastq or fasta input files and writes corrected reads in the same format. 
    // Fastq files include qual scores. 
    // Fasta files may have their quals in a separate .qual file (and these are expected to be lists of comma-separated integers)
    //
    // -r hg -m 50 -good 80 -t 8 Cspor_25.cbt  s_1_?_sequence.fastq
    // 
    // V1.2 introduces adapter trimming, read trimming/extension and discarding erroneous reads (previous default was to keep such reads unchanged)
    //

    // types of corrections made to a k-mer
    public enum FixTypes
    {
        fixNone = 0,                        // no change (default)
        fixSub = 1,                         // substitution 
        fixDel = 2,                         // fix deletion by inserting
        fixIns = 3,                         // fix insertion by deleting
        fixN = 4,                           // fix by replacing N
        fixAbandon = 5,                     // abandoned attempt at fixing mer
        maxFixTypes = 6
    }

    // results of checking whether read should be healed (and what happened when we tried healing it)
    public enum ReadState
    {
        readUnknown,                        // initial state
        readOK,                             // read seems to contain no errors
        readNeedsChecking,                  // read seems to be OK but still needs checking (depth drop etc)
        readIsBroken,                       // read (still) contains bad kMers/pairs
        readCorrected,                      // read was broken but is now OK
        readTooDeep,                        // read was deeper than requested max depth (parameter)
        readNotLongEnough                   // read didn't pass the minimum good kMers test (trimmmed too much//too many corrections)
    }

    // reasons for calling TryHealingMer
    enum kMerState
    {
        merUnknown,                         // not yet determined
        merOK,                              // not checking this mer at all
        merCheck,                           // kMer looks to be ok but needs checking
        merUnsure,                          // below OK or asymmetric counts
        merBad                              // below minDepth
    }

    enum PairDirection
    {
        backwards,                          // m points to the start of the k-mer containing the second fragment
        forwards                            // m points to the start of the k-mer containing the first fragment
    }

    enum PairWanted
    {
        backOnly = 0,                       // only try for 'back' pairs
        either = 1                          // try for backwards first, then forwards
    }

    enum SaveOption
    {
        saveGoodOnly = 1,                   // only good reads are saved
        saveProblemsAsWell = 2              // save problem reads as well - to separate 'problems' file
    }

    enum ErrorModel
    {
        mostlySubs = 1,                     // largely sub errors and the occasional rare indel
        indelsCommon = 2                    // indels common after homopolymer runs so do extra checks  
    }

    enum Tracing
    {
        traceOff = 0,                       // no tracing
        traceChanges = 1,                   // tracing healing at a high level (default trace)
        traceChoices = 2,                   // trace choices made (and variants to choose from)
        traceRead = 3,                      // detailed tracing for a specified read (set by check for a read)
        traceFollowers = 4                  // trace out the recursive followers code (set manually if needed)
    }

    enum VariantWanted
    {
        varyLast,                           // vary last base only
        varyAnyOne,                         // allow any single base to vary
        varyAnyTwo                          // allow any two bases to vary
    }

    enum ProgramPhase                       // states for monitor thread reporting
    {
        loading = 0,                        // loading kMers and pairs
        healing = 1                         // correcting reads
    }

    public enum AbandonReason
    {
        notAbandoned,
        tooManyNs,
        rewriting,
        treeSize,
        noNextMer
    }

    class Program
    {
        const string version = "2.2.0";             // version - update on every significant change

        // internal performance analysis options
        const bool perfTrace = false;               // save performance stats for each read

        // monitoring/reporting progress
        const int monitorInterval = 100;            // monitor thread progress at this interval
        const int reportInterval = 60000;           // report back at this interval (ms)
        static bool stopReporter = false;           // tell monitor thread to finish
        static ManualResetEvent signalReporter = new ManualResetEvent(false);
        // rough (non-mutexed) stats for progress reporting
        static long progressReads = 0;              // how many reads have been processed so far
        static long progressHealedReads = 0;        // how many reads were healed
        static long progressOKReads = 0;            // how many reads looked OK
        static ProgramPhase programPhase = ProgramPhase.loading;

#if STATS
        static long smdCalls = 0;
        static long spdCalls = 0;
        static long gpdCalls = 0;
        static long gpdFast = 0;
        static long gpdRead = 0;
        static long gpdNoRoom = 0;
        //static int[] hdubLocs = new int[200];
        static long depthBadButRedeemed = 0;
        static Dictionary<string, int> MIHCCauses = new Dictionary<string, int>(10);
        static long gdcCalls = 0;
        static long gds1Calls = 0;
        static long gds4Calls_THR = 0;
        static long gds4Calls_Peek = 0;
        static long gds4Calls_StartIsBad = 0;
        static long gds4Calls_Alts = 0;
        static long gds4Calls_Next = 0;
        static long gds4Calls_TrimAbs = 0;
        static long gds4Calls_FPV = 0;
#endif

        //static Stopwatch smdTimer = new Stopwatch();
        //static Stopwatch spdTimer = new Stopwatch();
        //static Stopwatch gmdTimer = new Stopwatch();
        //static Stopwatch gpdTimer = new Stopwatch();

        const int maxMerSize = 32;                  // max bases that can fit in a single ulong
        const int batchSize = 1000;                 // how many reads are fetched by each healing thread for a batch
        const int defaultReadLength = 200;          // default char[] length (resized as needed)
        const int defaultHeaderLength = 50;         // default header length

        const int maxNs = 3;                        // max no. of N's in a single read (avoid combinatoric explosion)
        const int maxFollowerRepairs = 5;           // max no. of repairs that can be made while trying to find followers
        const int maxFollowerRepairsTail = 2;       // max no. of repairs when counting followers in error-rich tail
        const int maxTHMAllowed = 1000;             // doesn't need to be too big, just decreases performance and doesn't add to correction rate
        const int almostEnd = 2;                    // how many bases from the end is 'close'
        const int goodMerRunLength = 4;             // how many consecutive OK kMers are considered a 'good run'
        static int[] gapCosts = new int[] { 2, 1, 1, 0, 0 }; // costs related to gaps between fixes - close fixes are penalised
        const int maxConsecutiveNFixesAllowed = 10; // max no.of consecutive Ns that can be fixed
        const int maxGap = 5;                       // no. of consecutive Ins (deletions) allowed in GenerateInsVariants
        const int highDepthFactor = 10;             // high depth reads are this many times the average
        const int veryHighDepthFactor = 100;        // and very deep reads are... 
        static int rewriteRegion = 20;              // length of region being checked for rewriting

        static ErrorModel errorModel = ErrorModel.mostlySubs;    // what type of data are we dealing with
        static bool onlySubsAllowed = false;        // never try for indels (should be slightly less accurate but faster)

        // (default) reads are always trimmed of bad bases after correction and may be extended (while the next base is unambiguous) to either the starting length or longer
        static bool readsFixedLength = false;       // true -> reads are fixed length and should be the same length when written out. false -> reads lengths can change  
        static bool readsFixedLengthPadded = false; // not only fixed length but padded with Ns if necessary  
        static int extendBases = 0;                 // extend reads by up to this many additional bases (after trimming)
                                                    // only sufficiently 'good' reads are written to the the corrected reads file (not too many fixes, not too short after trimming)
        static SaveOption saveReads = SaveOption.saveGoodOnly; // default is to save only good enough reads
        static int saveGoodMinLengthPct = 70;       // % of initial read that must have been 'good' after correction
        static int wantedReadLength = -1;           // trim reads to this length before correction
        static int minOKQual = 30;                  // poor/OK qual boundary used when determining whether we're now in the error-prone tail
        static int balancedFactor = 10;             // how much a read pair need to be out of balance before being considered unbalanced

        static int readsFormat = SeqFiles.formatNone; // what format are the reads files (see MerString.formatXXX enumeration)
        // this will be set automatically by looking at the file names/content if not explicitly set by a parameter
        static bool pairingSpecified = false;       // pairing parameter set explicitly - so don't set by default
        static bool pairedReads = false;            // do we correct reads in sets or singly? (true by default if we have multiple input files)
        static bool fullQualHeaders = false;        // do the reads use full FASTQ qual headers or just "+"?
        static int qualBase = 0;                    // qual offset for fastq reads. Set to 33 or 64 once format is known. Set to 0 for Fasta files
        static char replacementQual = (char)35;     // qual value used for fixed/inserted bases. (0-40)

        const bool findBiggest = true;              // look for biggest in list (for FindBestValue)
        const bool findSmallest = false;            // look for smallest in list (for FindBestValue)

        static StreamWriter perf = null;
        //static StreamWriter trace = new StreamWriter("trace.txt");
        //static StreamWriter depthsTrace = new StreamWriter("readDepths.txt");

        static string myProcessNameAndArgs;

        // precise stats (accumulated locally & merged under lock protection)
        static HealingStats stats = new HealingStats();
        static Stopwatch timeSpentReading = new Stopwatch();
        static Stopwatch[] threadWaitingRead;
        static Stopwatch[] threadWaitingWrite;
        static Stopwatch[] threadCorrecting;
        static double totalReadingTime = 0.0;
        static double totalWritingTime = 0;
        static double totalCorrectingTime = 0;

#if STATS
        static double timeInChecking = 0.0;
        static double timeInTHRFirstPass = 0.0;
        static double timeInTHRTHPass = 0.0;
        static double timeInTHRReversePass = 0.0;
        [ThreadStatic]
        static Stopwatch threadInChecking = new Stopwatch();
        [ThreadStatic]
        static Stopwatch threadInTHRFirstPass = new Stopwatch();
        [ThreadStatic]
        static Stopwatch threadInTHRTHPass = new Stopwatch();
        [ThreadStatic]
        static Stopwatch threadInTHRReversePass = new Stopwatch();
#endif

        static MerComparerClass merVariantComparer;

        const int traceReadLimit = 10000;                       // stop trace after this many reads have been healed
        static StreamWriter traces;                             // trace file
        static string tracesFN;                                 // and its file name
        static Dictionary<string, string> refSequences;         // possible set of known-to-be-correct sequences correspondong to those being corrected
        static Tracing tracing = Tracing.traceOff;              // level of tracing wanted
        static List<string> traceClump = new List<string>(100); // traces are written to this list before being written to the trace file.
        static bool tracingThisRead = false;                    // flag to say we're in the trapped read

        static Object traceLock = new Object();                 // lock object used to single-thread writing out traces
        static Object statsLock = new Object();                 // lock object used when merging local stats with global ones

        // trace string consts - always keep in sync with corresponding enums
        static char[] fixTypes = new char[] { '-', 'S', 'D', 'I', 'N', 'A' };
        static string[] fixNames = new string[] { "None", "Sub", "Del", "Ins", "N", "Abn" };

        // shared application-wide data. 
        // ----------------------------
        // These variables are read-only after they've been initialised
        // They are declared static to avoid passing them as (unchanging) parameters to a number of highly-called methods.
        static int noHealingThreads = 2;            // how many concurrent healing threads are wanted
        static kMerTable kMersTable = null;         // the collection of tiled reads from the .cbt file
        static int merSize = 0;                     // how big the mers are in the tiled files
        static int averageLoadedMerDepth = 0;       // and the average depth of coverage for all loaded mers
        static int requestedMinDepth = 0;           // (param) default reps count between bad & poor
        static int requestedMaxDepth = int.MaxValue;// (param) ignore reads with an average depth greater than this
        static int minLoadDepth = 0;                // min rep count needed before mer will be loaded into uniqueMers table (derived from minReps)
        static PairTable pairsTable = null;         // the set of kMer pairs read from the .prs file
        const int merPairFragmentLength = 16;       // pairs of 16-mers
        static int merPairGap = 0;                  // gap between the two 16-mers in a pair
        static int merPairLength = 0;               // length to go backwards to the start of a read
        static int averageLoadedPairDepth = 0;      // average depth of loaded pairs
        static bool amplicons = false;              // data is (PacBio) amplicons. No HDUB, no balanced checks

        const int lowComplexityBases = 6;           // change this to change the number of HP bases considered low complexity
        static HashSet<ulong> lowComplexityTrap;    // small hash set of low-complexity fragments. N bases only. Used to detect failure in LC regions
        static ulong lowComplexityMask;             // masks just the low N bases of the kMer


        [ThreadStatic]
        static Queue<kMerProperties> freeMerProperties;     // recyclable merProperties objects - one per thread to avoid locking overhead
        [ThreadStatic]
        static List<kMerProperties> allocatedMerProperties; // used to return merProperties objects after a THM exception
        [ThreadStatic]
        static List<Depth[]> freeDepths;                    // recyclable depth arrays objects - one per thread to avoid locking overhead
        [ThreadStatic]
        static List<Depth[]> allocatedDepths;               // return from here on exception
        [ThreadStatic]
        static Queue<List<kMerProperties>> freeVariantSets; // recyclable variantSets (THM)
        [ThreadStatic]
        static List<List<kMerProperties>> allocatedVariantSets;          // return from here on exception
        [ThreadStatic]
        static List<MerVariant[]> freeMerVariants;          // recyclable merVariant arrays (THM) - may contain arrays of differing lengths
        [ThreadStatic]
        static List<MerVariant[]> allocatedMerVariants;          // return from here on exception
        [ThreadStatic]
        static Queue<Sequence> freeSequences;               // recyclable sequences (reads and quals) for read contexts and RC healing

        //static int mpAllocated = 0;
        //static int seqAllocated = 0;
        //static int depthsAllocated = 0;
        //static int mvAllocated = 0;
        //static int vsAllocated = 0;

        public static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                WriteUsage();
                return;
            }

            string runName = null;               // tag used for stats/trace/corrected files created by this run
            string cbtFN = null;                 // name of .cbt file (and used for create .prs pairs file name)
            string statsFN = null;               // name of stats file (default is cbtFN + "stats.txt")
            string outputDir = null;             // output directory (default is current directory)
            List<string> FNParams = new List<string>();    // the .cbt name and the set of file names or patterns

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-h" || args[p] == "-help")
                    {
                        WriteUsageFull();
                        return;
                    }

                    if (args[p] == "-r" || args[p] == "-run")
                    {
                        if (!CheckForParamValue(p, args.Length, "run name string expected after -r|-run"))
                            return;
                        runName = args[p + 1];
                        p++;
                        //Console.WriteLine("-run " + runName);
                        continue;
                    }

                    if (args[p] == "-s" || args[p] == "-stats")
                    {
                        if (!CheckForParamValue(p, args.Length, "stats file name string expected after -s|-stats"))
                            return;
                        statsFN = args[p + 1];
                        p++;
                        //Console.WriteLine("-stats " + statsFN);
                        continue;
                    }

                    if (args[p] == "-m" || args[p] == "-min" || args[p] == "-minReps")
                    {
                        if (!CheckForParamValue(p, args.Length, "minReps number expected after -m|-min"))
                            return;
                        try
                        {
                            requestedMinDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -m|-min parameter: " + args[p + 1]);
                            return;
                        }
                        if (requestedMinDepth < 1)
                        {
                            Console.WriteLine("requested minimum depth must be >= 1");
                            return;
                        }
                        p++;
                        //Console.WriteLine("-min " + requestedMinDepth);
                        continue;
                    }

                    if (args[p] == "-max" || args[p] == "-maxDepth")
                    {
                        if (!CheckForParamValue(p, args.Length, "maxDepth number expected after -max|-maxDepth"))
                            return;
                        try
                        {
                            requestedMaxDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -max|-maxDepth parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        //Console.WriteLine("-max " + requestedMaxDepth);
                        continue;
                    }

                    if (args[p] == "-b" || args[p] == "-balanced")
                    {
                        if (!CheckForParamValue(p, args.Length, "balanced ratio expected after -b|-balanced"))
                            return;
                        try
                        {
                            balancedFactor = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -b|-balanced parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        //Console.WriteLine("-balanced " + balancedFactor);
                        continue;
                    }

                    if (args[p] == "-l" || args[p] == "-length")
                    {
                        if (!CheckForParamValue(p, args.Length, "trimmed length expected after -l|-length"))
                            return;
                        try
                        {
                            wantedReadLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -l|-length parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        //Console.WriteLine("-length " + wantedReadLength);
                        continue;
                    }

                    if (args[p] == "-f" || args[p] == "-format")
                    {
                        if (!CheckForParamValue(p, args.Length, "reads format expected after -f|-format"))
                            return;
                        string readsFormatParam = args[p + 1].ToLower();
                        if (readsFormatParam == "fna")
                            readsFormat = SeqFiles.formatFNA;
                        else if (readsFormatParam == "fasta")
                            readsFormat = SeqFiles.formatFNA;
                        else if (readsFormatParam == "fa")
                            readsFormat = SeqFiles.formatFNA;
                        else if (readsFormatParam == "fastq")
                            readsFormat = SeqFiles.formatFASTQ;
                        else if (readsFormatParam == "fq")
                            readsFormat = SeqFiles.formatFASTQ;
                        else
                        {
                            Console.WriteLine("reads format must be fasta or fastq: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-hp")
                    {
                        errorModel = ErrorModel.indelsCommon;
                        //Console.WriteLine("-hp");
                        continue;
                    }

                    if (args[p] == "-subsonly")
                    {
                        onlySubsAllowed = true;
                        //Console.WriteLine("-subs");
                        continue;
                    }

                    if (args[p] == "-paired")
                    {
                        pairingSpecified = true;
                        pairedReads = true;
                        //Console.WriteLine("-paired");
                        continue;
                    }

                    if (args[p] == "-unpaired")
                    {
                        pairingSpecified = true;
                        pairedReads = false;
                        //Console.WriteLine("-unpaired");
                        continue;
                    }

                    if (args[p] == "-fixed")
                    {
                        readsFixedLength = true;
                        //Console.WriteLine("-fixed");
                        continue;
                    }

                    if (args[p] == "-fixedPadded" || args[p] == "fp")
                    {
                        readsFixedLength = true;
                        readsFixedLengthPadded = true;
                        //Console.WriteLine("-fixedPadded");
                        continue;
                    }

                    if (args[p] == "-v" || args[p] == "-variable")
                    {
                        readsFixedLength = false;
                        //Console.WriteLine("-variable");
                        continue;
                    }

                    if (args[p] == "-amplicons")
                    {
                        amplicons = true;
                        //Console.WriteLine("-amplicons");
                        continue;
                    }

                    if (args[p] == "-g" || args[p] == "-good")
                    {
                        if (!CheckForParamValue(p, args.Length, "%good expected after -g|-good"))
                            return;
                        int saveGoodParam = 0;
                        try
                        {
                            saveGoodParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -g|-good parameter: " + args[p + 1]);
                            return;
                        }

                        if (saveGoodParam > 100)
                            saveGoodParam = 100;
                        if (saveGoodParam < 0)
                            saveGoodParam = 0;
                        saveGoodMinLengthPct = saveGoodParam;

                        p++;
                        //Console.WriteLine("-good " + saveGoodParam);
                        continue;
                    }

                    if (args[p] == "-mq" || args[p] == "-minqual")
                    {
                        if (!CheckForParamValue(p, args.Length, "min qual score expected after -mq|-minqual"))
                            return;
                        int saveMQParam = 0;
                        try
                        {
                            saveMQParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -mq|-minqual parameter: " + args[p + 1]);
                            return;
                        }

                        if (saveMQParam > 40)
                            saveMQParam = 40;
                        if (saveMQParam < 0)
                            saveMQParam = 0;
                        minOKQual = saveMQParam;

                        p++;
                        //Console.WriteLine("-mq " + saveMQParam);
                        continue;
                    }

                    if (args[p] == "-extend")
                    {
                        if (!CheckForParamValue(p, args.Length, "#bases expected after -extend"))
                            return;
                        int extendParam = 0;
                        try
                        {
                            extendParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number after the -extend parameter: " + args[p + 1]);
                            return;
                        }

                        if (extendParam > 100)
                            extendParam = 100;
                        if (extendParam < 0)
                            extendParam = 0;

                        extendBases = extendParam;
                        continue;
                    }

                    if (args[p] == "-trace")
                    {
                        tracing = Tracing.traceChanges;
                        continue;
                    }
                    if (args[p] == "-tracechanges")
                    {
                        tracing = Tracing.traceChanges;
                        continue;
                    }
                    if (args[p] == "-tracechoices")
                    {
                        tracing = Tracing.traceChoices;
                        continue;
                    }

                    if (args[p] == "-problems")
                    {
                        saveReads = SaveOption.saveProblemsAsWell;
                        continue;
                    }

                    if (args[p] == "-t" || args[p] == "-threads")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -t|-threads"))
                            return;
                        try
                        {
                            noHealingThreads = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -t|-threads parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        //Console.WriteLine("-threads " + noHealingThreads);
                        continue;
                    }

                    if (args[p] == "-o" || args[p] == "-output")
                    {
                        if (!CheckForParamValue(p, args.Length, "directory name expected after -o|-output"))
                            return;
                        outputDir = args[p + 1];
                        p++;
                        //Console.WriteLine("-output " + outputDir);
                        continue;
                    }

                    Console.WriteLine("unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                FNParams.Add(args[p]);
            }

            if (runName == null)
                runName = "corrected_" + requestedMinDepth;

            if (requestedMinDepth == 0)
            {
                Console.WriteLine("no minimum k-mer depth specified (-m)");
                return;
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("expected a cbt file name and at least one reads file name or pattern");
                return;
            }

            if (noHealingThreads < 1)
            {
                Console.WriteLine("noThreads must be at least 1 - " + noHealingThreads);
                return;
            }

            // validate the output directory & set the output prefix string
            string fnSeparator = Path.DirectorySeparatorChar.ToString();  // \ for Windows; / for Unix/Linux
            if (outputDir != null)
            {
                try
                {
                    // add a trailing \ if the output directory name doesn't already have one
                    if (!Directory.Exists(outputDir))
                        Directory.CreateDirectory(outputDir);
                    if (!outputDir.EndsWith(fnSeparator))
                        outputDir += fnSeparator;
                    string testOutputFN = outputDir + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                    StreamWriter testTemp = new StreamWriter(testOutputFN);
                    testTemp.Close();
                    File.Delete(testOutputFN);
                }
                catch
                {
                    Console.WriteLine("Output directory: " + outputDir + " was invalid");
                    return;
                }
            }

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names or patterns");
                return;
            }

            // take the cbt file name from the start of the non-option list
            cbtFN = FNParams[0];
            FNParams.RemoveAt(0);

            if (FNParams.Count == 0)
            {
                Console.WriteLine("did not find any reads file names");
                return;
            }

            // calculate the min load depth from the min reps depth - don't need to load all of the singletons and other errors into memory
            minLoadDepth = requestedMinDepth / 5;
            if (minLoadDepth <= 1)
                minLoadDepth = 2;

            // find out who we are so we can track what program & args produced the result files
            Process myProcess = Process.GetCurrentProcess();
            myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;
            myProcessNameAndArgs += " (" + version + ") " + DateTime.Now.ToString();

            Console.WriteLine(myProcessNameAndArgs);

            if (!File.Exists(cbtFN))
            {
                Console.WriteLine("k-mer consensus (.cbt) file not found: " + cbtFN);
                return;
            }

            List<string> readsFileNames = new List<string>(FNParams.Count);
            List<string> readsFilePaths = new List<string>(FNParams.Count);
            foreach (string readsFNP in FNParams)
            {
                string readsFileName;
                string readsFilePath;
                GetPathFN(readsFNP, out readsFilePath, out readsFileName);
                readsFilePaths.Add(readsFilePath);
                readsFileNames.Add(readsFileName);
            }

            List<string> expandedReadsFNs = new List<string>();
            for (int f = 0; f < FNParams.Count; f++)
            {
                string[] matchedReadsFNs = Directory.GetFiles(readsFilePaths[f], readsFileNames[f], SearchOption.TopDirectoryOnly);
                foreach (string matchedReadsFN in matchedReadsFNs)
                    expandedReadsFNs.Add(matchedReadsFN);
            }

            // make sure there aren't any duplicates in the file list (seems to be a bug on the Cherax SGI HPC system and it returns each file name twice)
            List<string> distinctReadsFNs = new List<string>();
            foreach (string fn in expandedReadsFNs)
                if (!distinctReadsFNs.Contains(fn))
                    distinctReadsFNs.Add(fn);

            // finally... the set of fully qualified, distinct reads files
            string[] readsFNs;
            readsFNs = distinctReadsFNs.ToArray();

            int noOfReadsFiles = distinctReadsFNs.Count;
            if (noOfReadsFiles == 0)
            {
                Console.WriteLine("No matching reads files found");
                return;
            }

            if (statsFN == null)
            {
                // construct a stats file name if one wasn't given
                int dotIdx = readsFileNames[0].LastIndexOf('.');
                if (dotIdx >= 0)
                {
                    statsFN = readsFileNames[0].Substring(0, readsFileNames[0].LastIndexOf('.'));
                    statsFN = statsFN.Replace('?', '_');
                    statsFN = statsFN.Replace('*', '_');
                    statsFN = statsFN.Replace('/', '_');
                    statsFN = statsFN.Replace('\\', '_');
                    statsFN = statsFN.Replace("__", "_");
                    statsFN = statsFN.Replace("__", "_");
                }
                else
                    statsFN = readsFileNames[0];
                statsFN = statsFN + "_" + runName + "_stats.txt";
                statsFN = statsFN.Replace("__", "_");
            }

            DateTime startHealing = DateTime.Now;

            tracesFN = runName + "_trace_" + merSize + "_" + requestedMinDepth + ".txt";
            if (tracing > Tracing.traceOff)
            {
                traces = new StreamWriter(outputDir + tracesFN);
                traces.WriteLine(myProcessNameAndArgs);
            }

            // load the .cbt file into a kMerTable
            //Console.WriteLine(GC.GetTotalMemory(false) + " memory before loadCBT");
            kMersTable = new kMerTable(cbtFN, minLoadDepth, 1000);
            merSize = kMersTable.merSize;
            averageLoadedMerDepth = kMersTable.averageDepthLoaded;
            //Console.WriteLine(GC.GetTotalMemory(false) + " memory after loadCBT");
            if (amplicons)
            {
                kMersTable.hdubFilter.Clear();
                rewriteRegion = 30;
            }

            if (kMersTable.kMersTable == null)
            {
                Console.WriteLine("kMer table load failed");
                return;
            }

            if (merSize == 0)
            {
                Console.WriteLine("bad k-mer size found at start of .cbt file");
                return;
            }

            if (tracing > Tracing.traceOff)
            {
                traces.WriteLine("loaded " + kMersTable.distinctMersLoaded + "/" + kMersTable.totalMersLoaded + " " + merSize + "-mers" + " avg depth = " + kMersTable.averageDepthLoaded);
                traces.WriteLine();
            }

            // and load the pairs file if it exists
            string pairsFN = cbtFN.Replace(".cbt", ".prs");

            if (File.Exists(pairsFN))
            {
                //pairsTable = new PairTable(pairsFN, minLoadDepth / 2);
                pairsTable = new PairTable(pairsFN, minLoadDepth);
                merPairGap = pairsTable.pairGap;
                merPairLength = pairsTable.pairFullLength;
                averageLoadedPairDepth = pairsTable.averageDepthLoaded;

                if (pairsTable.pairsTable == null)
                {
                    Console.WriteLine("pairs table load failed");
                    return;
                }

                if (tracing > Tracing.traceOff)
                {
                    traces.WriteLine("loaded " + pairsTable.totalPairsLoaded + " pairs");
                    traces.WriteLine();
                }
            }

            // load the low-complexity trap
            lowComplexityTrap = new HashSet<ulong>();
            lowComplexityMask = 0xffffffffffffffff >> (merSize - lowComplexityBases) * 2; // will be and'ed with kmer with zero trailing bits in ulong
            AddToLowComplexityTrap("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".Substring(0, merSize), lowComplexityTrap, lowComplexityMask);
            AddToLowComplexityTrap("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".Substring(0, merSize), lowComplexityTrap, lowComplexityMask);
            AddToLowComplexityTrap("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".Substring(0, merSize), lowComplexityTrap, lowComplexityMask);
            AddToLowComplexityTrap("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".Substring(0, merSize), lowComplexityTrap, lowComplexityMask);

            merVariantComparer = new MerComparerClass();

            // start the monitor/synchronising thread
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();
            programPhase = ProgramPhase.healing;

            if (perfTrace)
            {
                perf = new StreamWriter(outputDir + runName + "_perf.txt");
                perf.WriteLine("read#\tfaHeader\ttime\tcacheSize\treverse\tresult\tunhealedRead\tthmCalls\tpoor");
            }

#if STATS
            MIHCCauses.Add("bad", 0);
            MIHCCauses.Add("bad ubhp alt", 0);
            MIHCCauses.Add("unbalanced alt", 0);
            MIHCCauses.Add("check middling alt", 0);
            MIHCCauses.Add("OK middling no alt", 0);
            MIHCCauses.Add("unsure middling", 0);
            MIHCCauses.Add("check depthDrop alt", 0);
            MIHCCauses.Add("OK depthDrop no alt", 0);
            MIHCCauses.Add("check m==0 alt", 0);
            MIHCCauses.Add("check tryHarder alt", 0);
            MIHCCauses.Add("check peek alt", 0);
#endif

            // if we weren't told what format the reads are in, use the format from the first file 
            if (readsFormat == SeqFiles.formatNone)
                readsFormat = SeqFiles.DetermineFileFormat(readsFNs[0]);

            // resolve the FASTQ qual ambiguity by reading through quals until one is encountered that can only come from either of the alternative sets
            if (readsFormat == SeqFiles.formatFASTQ)
                qualBase = SeqFiles.ResolveFastqQualAmbiguity(readsFNs[0], out fullQualHeaders);

            // and check whether we've got Unix data so we can write out the corrected files in the same format
            //string lfConvention = SeqFiles.LFConvention(readsFNs[0]);

            // if we have an even number of files, assume they're paired by default
            if (noOfReadsFiles % 2 == 0 && !pairingSpecified)
                pairedReads = true;
            if (noOfReadsFiles > 1 && noOfReadsFiles % 2 == 1 && pairingSpecified && pairedReads)
            {
                Console.WriteLine("found " + noOfReadsFiles + " reads files but asked to process them in pairs");
                return;
            }

            int filesInSet = (pairedReads && noOfReadsFiles > 1) ? 2 : 1;

            for (int fp = 0; fp < noOfReadsFiles; fp += filesInSet)
            {
                StreamReader[] readsFiles = new StreamReader[filesInSet];
                StreamReader[] qualFiles = new StreamReader[filesInSet];            // will stay null for all but 454 reads+quals in FASTA format
                BufferedReader[] bufferedReadsFiles = new BufferedReader[filesInSet];

                StreamWriter[] healedReads = new StreamWriter[filesInSet];
                StreamWriter[] healedQuals = new StreamWriter[filesInSet];          // will stay null for all but 454+quals reads in FASTA format
                BufferedWriter[] bufferedHealedReads = new BufferedWriter[filesInSet];
                BufferedWriter[] bufferedHealedQuals = new BufferedWriter[filesInSet];

                StreamWriter[] singleReads = new StreamWriter[filesInSet];
                StreamWriter[] singleQuals = new StreamWriter[filesInSet];
                BufferedWriter[] bufferedSingleReads = new BufferedWriter[filesInSet];
                BufferedWriter[] bufferedSingleQuals = new BufferedWriter[filesInSet];

                StreamWriter[] problemReads = null;
                StreamWriter[] problemQuals = null;
                BufferedWriter[] bufferedProblemReads = new BufferedWriter[filesInSet];
                BufferedWriter[] bufferedProblemQuals = new BufferedWriter[filesInSet];
                if (saveReads == SaveOption.saveProblemsAsWell)
                    problemReads = new StreamWriter[filesInSet];

                // open original and healed streams for each read file (and qual if necessary)
                for (int f = 0; f < filesInSet; f++)
                {
                    string fullReadsFN = readsFNs[fp + f];
                    string readsPath;
                    string readsFN;
                    GetPathFN(fullReadsFN, out readsPath, out readsFN);
                    string fileSuffix = readsFN.Substring(readsFN.LastIndexOf('.'));
                    string fileWithoutSuffix = readsFN.Substring(0, readsFN.LastIndexOf('.'));

                    readsFiles[f] = SeqFiles.OpenSeqStream(fullReadsFN);
                    bufferedReadsFiles[f] = new BufferedReader(readsFormat, readsFiles[f], qualFiles[f]);

                    Console.WriteLine("correcting " + fullReadsFN);

                    if (readsFormat == SeqFiles.formatFNA)
                    {
                        string qualFN = fileWithoutSuffix + ".qual";
                        if (File.Exists(qualFN))
                            qualFiles[f] = SeqFiles.OpenSeqStream(qualFN);
                    }

                    string outputPath = outputDir == null ? readsPath + fnSeparator : outputDir;

                    int initialBufferSize = batchSize * (defaultHeaderLength + 1 + defaultReadLength * 2 + 4 * 2); // typical fastq

                    healedReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + fileSuffix, false, readsFiles[f].CurrentEncoding, 10000);
                    bufferedHealedReads[f] = new BufferedWriter(healedReads[f], initialBufferSize, noHealingThreads, "healed");

                    if (pairedReads)
                    {
                        singleReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_singles" + fileSuffix);
                        bufferedSingleReads[f] = new BufferedWriter(singleReads[f], initialBufferSize/5, noHealingThreads, "singles");
                    }

                    if (saveReads == SaveOption.saveProblemsAsWell)
                    {
                        problemReads[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_problems" + fileSuffix);
                        bufferedProblemReads[f] = new BufferedWriter(problemReads[f], initialBufferSize/5, noHealingThreads, "problems");
                    }

                    if (qualFiles[f] != null)
                    {
                        healedQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + ".qual");
                        bufferedHealedQuals[f] = new BufferedWriter(healedQuals[f], initialBufferSize/5, noHealingThreads, "healedQuals");

                        if (pairedReads)
                        {
                            singleQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_singles" + ".qual");
                            bufferedSingleQuals[f] = new BufferedWriter(singleQuals[f], initialBufferSize/5, noHealingThreads, "singleQuals");
                        }

                        if (saveReads == SaveOption.saveProblemsAsWell)
                        {
                            problemQuals[f] = new StreamWriter(outputPath + fileWithoutSuffix + "_" + runName + "_problems" + ".qual");
                            bufferedProblemQuals[f] = new BufferedWriter(problemQuals[f], initialBufferSize/5, noHealingThreads, "problemQuals");
                        }
                    }
                }

                // load reference sequences corresponding to reads if these exist
                refSequences = new Dictionary<string, string>();
                foreach (string readFN in readsFNs)
                {
                    string suffix = readFN.Substring(readFN.LastIndexOf('.'));
                    string refFN = readFN.Replace(suffix, ".ref");
                    if (File.Exists(refFN))
                        LoadRefSequences(refFN, refSequences);
                }

                HealingThreadParams[] healingParams = new HealingThreadParams[noHealingThreads];
                Thread[] healingThreads = new Thread[noHealingThreads];

                threadWaitingRead = new Stopwatch[noHealingThreads];
                threadWaitingWrite = new Stopwatch[noHealingThreads];
                threadCorrecting = new Stopwatch[noHealingThreads];

                // ready a new thread for each parallel healer
                for (int b = 0; b < noHealingThreads; b++)
                {
                    healingParams[b] = new HealingThreadParams();
                    healingParams[b].threadNumber = b + 1;
                    healingParams[b].bufferedReadsFiles = bufferedReadsFiles;
                    healingParams[b].bufferedHealedReads = bufferedHealedReads;
                    healingParams[b].bufferedProblemReads = bufferedProblemReads;
                    healingParams[b].bufferedSingleReads = bufferedSingleReads;
                    healingParams[b].bufferedHealedQuals = bufferedHealedQuals;
                    healingParams[b].bufferedProblemQuals = bufferedProblemQuals;
                    healingParams[b].bufferedSingleQuals = bufferedSingleQuals;
                    healingThreads[b] = new Thread(new ParameterizedThreadStart(Program.HealingThread));
                    healingThreads[b].Priority = ThreadPriority.BelowNormal;
                    healingThreads[b].Name = b.ToString();
                    healingThreads[b].Start(healingParams[b]);
                    Thread.Sleep(100);
                    if (tracing > Tracing.traceOff)
                    {
                        // single-thread when tracing to give reproducible ordering in trace files
                        healingThreads[b].Join();
                        //Console.WriteLine("finished healing thread " + b);
                    }
                }

                // and wait for all threads to finish (no need to wait if tracing as this is single-threaded)
                if (tracing == Tracing.traceOff)
                    for (int b = 0; b < noHealingThreads; b++)
                    {
                        healingThreads[b].Join();
                        healingThreads[b] = null;
                        //Console.WriteLine("finished healing thread " + b);
                    }

                foreach (BufferedWriter bw in bufferedHealedReads)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedHealedQuals)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedHealedQuals)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedProblemReads)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedProblemQuals)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedSingleReads)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (BufferedWriter bw in bufferedSingleQuals)
                    if (bw != null)
                        bw.CloseBufferedWriter();

                foreach (StreamWriter r in healedReads)
                    if (r != null)
                        r.Close();
                if (healedQuals != null)
                    foreach (StreamWriter q in healedQuals)
                        if (q != null)
                            q.Close();
                if (problemReads != null)
                    foreach (StreamWriter r in problemReads)
                        if (r != null)
                            r.Close();
                if (problemQuals != null)
                    foreach (StreamWriter q in problemQuals)
                        if (q != null)
                            q.Close();
                if (singleReads != null)
                    foreach (StreamWriter s in singleReads)
                        if (s != null)
                            s.Close();
                if (singleQuals != null)
                    foreach (StreamWriter q in singleQuals)
                        if (q != null)
                            q.Close();
                if (perf != null)
                    perf.Close();

                for (int b = 0; b < noHealingThreads; b++)
                {
                    totalReadingTime += threadWaitingRead[b].Elapsed.TotalSeconds;
                    totalWritingTime += threadWaitingWrite[b].Elapsed.TotalSeconds;
                    totalCorrectingTime += threadCorrecting[b].Elapsed.TotalSeconds;
                }

            } // end of file (pair) processing loop

            BufferedWriter.FinishBufferWriter();

            DateTime endHealing = DateTime.Now;
            StreamWriter statsFile = new StreamWriter(outputDir + statsFN);

            statsFile.WriteLine(myProcessNameAndArgs);
            statsFile.WriteLine();
            StatsWrite(statsFile, stats.readsRead + "\treads in initial files");
            StatsWrite(statsFile, stats.OKReadsWritten + "\treads OK");
            StatsWrite(statsFile, stats.correctedReadsWritten + "\treads corrected");
            StatsWrite(statsFile, stats.discardedBroken + "\treads dropped (uncorrectable, too deep or too short after trimming)");
            StatsWrite(statsFile, stats.shortReadsFound + "\treads dropped (too short to correct)");
            StatsWrite(statsFile, stats.discardedOK + "\treads OK but dropped (pair not good enough)");
            statsFile.WriteLine();

            StatsWrite(statsFile, stats.OKReads + "\treads OK on first inspection");
            StatsWrite(statsFile, stats.OKReadsChecked + "\treads needed further checking but were OK");
            StatsWrite(statsFile, stats.healedFirstPass + "\treads healed in first pass");
            StatsWrite(statsFile, stats.healedAltPass + "\treads healed in second pass");
            StatsWrite(statsFile, stats.healedRCPass + "\treads healed in reverse pass");
            StatsWrite(statsFile, stats.readsNotHealedAtAll + "\treads broken but not healed");
            StatsWrite(statsFile, stats.tooDeepReads + "\treads skipped (too deep)");
            StatsWrite(statsFile, stats.tooShortReads + "\treads dropped (too short/poor after trimming");
            StatsWrite(statsFile, stats.trimmed + "\treads trimmed");
            StatsWrite(statsFile, stats.extended + "\treads extended");

            StatsWrite(statsFile, stats.abandonedNs + "\thealing passes abandoned (too many Ns)");
            StatsWrite(statsFile, stats.abandonedRewriting + "\thealing passes abandoned (rewriting)");
            StatsWrite(statsFile, stats.abandonedTree + "\thealing passes abandoned (tree size limit)");
            StatsWrite(statsFile, stats.abandonedFutile + "\thealing passes abandoned (uncorrectable kMer)");
            StatsWrite(statsFile, stats.hdubReads + "\tadapter-like reads found (very deep & unbalanced)");
            statsFile.WriteLine();

            StatsWrite(statsFile, stats.mers + "\tkMers in initial files");
            StatsWrite(statsFile, stats.replacedMers + "\tkMers fixed");
            StatsWrite(statsFile, stats.fixesByType[(int)FixTypes.fixSub] + "\tsubs");
            StatsWrite(statsFile, stats.fixesByType[(int)FixTypes.fixDel] + "\tdels");
            StatsWrite(statsFile, stats.fixesByType[(int)FixTypes.fixIns] + "\tins");
            statsFile.WriteLine();

            StatsWrite(statsFile, (endHealing - startHealing).TotalSeconds.ToString("#.0") + "\tsecs elapsed");
            StatsWrite(statsFile, totalCorrectingTime.ToString("#.0") + "\tsecs correcting reads");
            StatsWrite(statsFile, totalReadingTime.ToString("#.0") + "\tsecs reading reads");
            StatsWrite(statsFile, totalWritingTime.ToString("#.0") + "\tsecs writing reads");

#if STATS
            StatsWrite(statsFile, smdCalls + "\tSetMerDepths calls");
            StatsWrite(statsFile, spdCalls + "\tSetPairDepths calls");
            StatsWrite(statsFile, gpdCalls + "\tGetPairDepth calls");
            StatsWrite(statsFile, gpdNoRoom + "\tGetPairDepth(no room) calls");
            StatsWrite(statsFile, gpdRead + "\tGetPairDepth(read) calls");
            StatsWrite(statsFile, gpdFast + "\tGetPairDepth(fast) calls");
            StatsWrite(statsFile, gdcCalls + "\tkmt:GetDepthCounts calls");
            StatsWrite(statsFile, gds1Calls + "\tkmt:GetDepthSum(1) calls");
            StatsWrite(statsFile, gds4Calls_THR + "\tTHR:GetDepthSum(4)  calls");
            StatsWrite(statsFile, gds4Calls_Peek + "\tPeek:GetDepthSum(4) calls");
            StatsWrite(statsFile, gds4Calls_StartIsBad + "\tStartIsBad:GetDepthSum(4) calls");
            StatsWrite(statsFile, gds4Calls_Alts + "\tAlts:GetDepthSum(4) calls");
            StatsWrite(statsFile, gds4Calls_Next + "\tNext:GetDepthSum(4) calls");
            StatsWrite(statsFile, gds4Calls_TrimAbs + "\tTrimAbs:GetDepthSum(4) calls");
            StatsWrite(statsFile, gds4Calls_FPV + "\tFPV:GetDepthSum(4) calls");
            StatsWrite(statsFile, gpdCalls + "\tkpt:GetPairDepth calls");
            StatsWrite(statsFile, depthBadButRedeemed + "\tdepthBadButRedeemed");

            foreach (KeyValuePair<string, int> kvp in MIHCCauses)
                StatsWrite(statsFile, kvp.Value + "\t" + kvp.Key);

            //StatsWrite(statsFile, gpdTimer.Elapsed.TotalSeconds.ToString("F4") + "\tgdp in SetMerDepths");
            //StatsWrite(statsFile, gpdTimer.Elapsed.TotalSeconds.ToString("F4") + "\tgpd in SetPairDepths");
            //StatsWrite(statsFile, smdTimer.Elapsed.TotalSeconds.ToString("F4") + "\tSetMerDepths");
            //StatsWrite(statsFile, spdTimer.Elapsed.TotalSeconds.ToString("F4") + "\tSetPairDepths");

            StatsWrite(statsFile, timeInChecking.ToString("F4") + "\tInitialScan");
            StatsWrite(statsFile, timeInTHRFirstPass.ToString("F4") + "\tFirst THR");
            StatsWrite(statsFile, timeInTHRTHPass.ToString("F4") + "\tTH THR");
            StatsWrite(statsFile, timeInTHRReversePass.ToString("F4") + "\tReverse THR");

            //for (int i = 0; i < 200; i++)
            //    if (hdubLocs[i] != 0)
            //        Console.Write(hdubLocs[i] + "@" + i + ", ");
            //Console.WriteLine();
#endif

            if (tracing > Tracing.traceOff)
                traces.Close();
            statsFile.Close();

            //depthsTrace.Close();

            stopReporter = true;
            signalReporter.Set();
            monitorProgress.Join();

            //Console.WriteLine("seqs " + seqAllocated + ", " + "mp " + mpAllocated + ", depths " + depthsAllocated + ", vs " + vsAllocated + ", mv " + mvAllocated);
            //trace.Close();
        }

        private static void StatsWrite(StreamWriter statsFile, string statsLine)
        {
            statsFile.WriteLine(statsLine);
            statsLine = statsLine.Replace('\t', ' ');
            Console.WriteLine(statsLine);
        }

        private static void AddToLowComplexityTrap(string lowComplexitySeq, HashSet<ulong> lowComplexityTrap, ulong lowComplexityMask)
        {
            string lowComplexityPrefix = lowComplexitySeq.Substring(0, merSize - 1);
            foreach (char b in new char[] { 'A', 'C', 'G', 'T' })
            {
                ulong lowComplexityMer = kMers.CondenseMer(lowComplexityPrefix+b, merSize);
                lowComplexityMer = lowComplexityMer & lowComplexityMask;
                if (!lowComplexityTrap.Contains(lowComplexityMer))
                    lowComplexityTrap.Add(lowComplexityMer);
            }
        }

        private static void GetPathFN(string readsFN, out string readsPath, out string readsFNP)
        {
            char FSC = Path.DirectorySeparatorChar;
            string FSS = new string(FSC, 1);
            readsPath = null;
            if (readsFN.Contains(FSS))
            {
                readsPath = readsFN.Substring(0, readsFN.LastIndexOf(FSC));
                readsFNP = readsFN.Substring(readsFN.LastIndexOf(FSC) + 1);
            }
            else
            {
                readsPath = Directory.GetCurrentDirectory();
                readsFNP = readsFN;
            }
        }

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Blue [-help] [-r run] -m minReps [-f fasta|fastq] [-hp] [-t threads] [-l length] [-fixed] [-variable] [-good nn%] [-problems] [-extend nn]" +
                                    " [-output dir] [-paired] [-unpaired] k-merFN readsFNs or patterns (" + version + ")");
        }

        private static void WriteUsageFull()
        {
            Console.WriteLine("usage: Blue <options> <k-mer file> <list of file names or patterns for reads files to be corrected>\n" +
                                "\t[-h|help] : display this usage help and exit\n" +
                                "\t[-r|run name] : tag inserted into corrected reads files, and name used for stats file. Default is corrected_<depth>\n" +
                                "\t -m|min minReps : min k-mer depth for OK when scanning reads (used only when dynamic calculation fails)\n" +
                                "\t[-f|format fasta or fastq] : Blue will try to work out the format from the file so this parameter is optional\n" +
                                "\t[-hp] : do additional checks for possible errors at the end of every homopolymer string - intended for 454 and IonTorrent\n" +
                                "\t[-t|threads #threads] : no. of threads to use for the correction. Default is 1 thread.\n" +
                                "\t[-l length] : trim reads to this length before correction.\n" +
                                "\t[-fixed] : fixed length reads, corrected reads are extended/trimmed to the same length as original reads (default).\n" +
                                "\t[-v|variable] : variable length reads, may be trimmed or extended during healing (default if -hp set).\n" +
                                "\t[-g|good nn%] : only save reads that are at least this long after correction and trimming. Default is 80%.\n" +
                                "\t[-problems] : save reads that couldn't be corrected to a 'problems' file.\n" +
                                "\t[-extend nn] : extend reads by up to nn bases where possible.\n" +
                                "\t[-paired]: reads from both files in a pair must be good/corrected (default with pairs of files)\n" +
                                "\t[-unpaired] : don't treat pairs of files as pair sets. Each file will be corrected separately.\n" +
                                "\t[-s|stats name] : name of stats file (constructed from first seq file name if not provided).\n" +
                                "\t[-o|output outputDir] : directory for corrected reads etc" +
                                "\t k-mer file : .cbt file created by Tessel. The same name is used to find the pairs file (.prs) if it exists\n");
        }

        private static bool CheckForParamValue(int p, int argsLength, string msg)
        {
            if (p == argsLength - 1)
            {
                Console.WriteLine(msg);
                return false;
            }
            return true;
        }

        private static void LoadRefSequences(string refFN, Dictionary<string, string> refSequences)
        {
            //>FUOCO2J01D9UB4 rank=0000217 x=1635.0 y=834.0 length=120
            //=TGCCGCACCTGCTGCGACAGGGTCGGCTGCGAAACGTGCAGCGCCTCGGCGGCGCGGGTGAAGTTGCGATGCTCGGCGACAGCGAGCAGGTAGCGCAGGTGGCGAAGCAGCATGGAAACATCC 

            //>FUOCO2J01EDV9U rank=0000223 x=1681.0 y=1552.0 length=105
            //=CCGCCAGCAACTTGGCTTCCAGCACTTCGTTGCTGTCGTAGACGTCGTAGACGACCTTGATCCCGGTTTCCTTGGTGAACTTCTCCAGGGTGTCCGGCGCGATGTAG 

            StreamReader refs = new StreamReader(refFN);
            bool EOF = false;
            string line;
            string[] parsedLine;
            char[] blankSeparator = new char[] { ' ' };

            while (!EOF)
            {
                line = refs.ReadLine();
                if (line == null)
                    break;

                if (line.Length > 1 && line[0] == '>')
                {
                    // found next header line
                    parsedLine = line.Split(blankSeparator);
                    string readID = parsedLine[0].Substring(1);     // strip off leading > or @
                    string refSequence = "";
                    bool within = true;

                    while (within)
                    {
                        // skip any lines other than the reference (=), until we hit it or are about to hit the next group
                        if (refs.Peek() == '>')
                            break;
                        refSequence = refs.ReadLine();
                        if (refSequence.Length > 1 && refSequence[0] == '=')
                        {
                            // found the first reference sequence for this read
                            refSequences.Add(readID, refSequence);
                            break;
                        }
                    }
                }
            }
        }

        private static void HealingThread(object threadParams)
        {
            //Console.WriteLine(GC.GetTotalMemory(false) + " memory on entry to HealingThread");
            HealingThreadParams theseParams = (HealingThreadParams)threadParams;
            int threadNumber = theseParams.threadNumber;                                // which healing thread is this one?
            BufferedReader[] bufferedReadsFiles = theseParams.bufferedReadsFiles;       // the buffered read layer over the (shared) reads files (and qual files if they exist)
            BufferedWriter[] bufferedHealedReads = theseParams.bufferedHealedReads;     // and the buffered healed reads (and quals)
            BufferedWriter[] bufferedSingleReads = theseParams.bufferedSingleReads;     // good but unpaired reads
            BufferedWriter[] bufferedProblemReads = theseParams.bufferedProblemReads;   // (if saving) problem reads
            BufferedWriter[] bufferedHealedQuals = theseParams.bufferedHealedQuals;     // these arrays are always allocated but will  
            BufferedWriter[] bufferedSingleQuals = theseParams.bufferedSingleQuals;     // contain null entries if not in-use
            BufferedWriter[] bufferedProblemQuals = theseParams.bufferedProblemQuals;

            // per-thread stats counters - added to global counters under lock control at end
            HealingStats threadStats = new HealingStats();
            // and the per-thread set of available (recycled) merProperties and depths objects
            freeMerProperties = new Queue<kMerProperties>(100);
            allocatedMerProperties = new List<kMerProperties>(100);
            freeDepths = new List<Depth[]>(100);
            allocatedDepths = new List<Depth[]>(100);
            freeMerVariants = new List<MerVariant[]>(100);
            allocatedMerVariants = new List<MerVariant[]>(100);
            freeVariantSets = new Queue<List<kMerProperties>>(100);
            allocatedVariantSets = new List<List<kMerProperties>>(100);
            freeSequences = new Queue<Sequence>(100);
            // per-thread readProperties structure - same one used for all reads
            ReadProperties readProps = new ReadProperties(500);

            threadWaitingRead[threadNumber - 1] = new Stopwatch();
            threadWaitingWrite[threadNumber - 1] = new Stopwatch();
            threadCorrecting[threadNumber - 1] = new Stopwatch();
#if STATS
            threadInChecking = new Stopwatch();
            threadInTHRFirstPass = new Stopwatch();
            threadInTHRTHPass = new Stopwatch();
            threadInTHRReversePass = new Stopwatch();
#endif
            int noReadFNs = bufferedReadsFiles.Length;
            bool[] fileActive = new bool[noReadFNs];                        // have not yet reached EOF on this reads file
            bool[][] readValid = new bool[noReadFNs][];                     // is this read valid or were we past EOF?
            int filesStillActive = noReadFNs;                               // how many active reads files are still around

            // allocate buffer arrays, but pointers to buffers may be null if corresponding files are not in use
            BufferedWriter.WriteBuffer[] healedReadsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];
            BufferedWriter.WriteBuffer[] healedQualsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];
            BufferedWriter.WriteBuffer[] singleReadsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];
            BufferedWriter.WriteBuffer[] singleQualsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];
            BufferedWriter.WriteBuffer[] problemReadsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];
            BufferedWriter.WriteBuffer[] problemQualsWriteBuffers = new BufferedWriter.WriteBuffer[noReadFNs];

            Sequence[][] headerSet = new Sequence[noReadFNs][];             // a set of read headers
            Sequence[][] readSet = new Sequence[noReadFNs][];               // a set of reads, possibly one from each file (length+char[])
            Sequence[][] qualHeaderSet = new Sequence[noReadFNs][];         // set of headers for the quals
            Sequence[][] qualsSet = new Sequence[noReadFNs][];              // set of quals (in canonical form)
            Sequence[][] healedReadSet = new Sequence[noReadFNs][];         // and the returned corrected reads
            Sequence[][] healedQualSet = new Sequence[noReadFNs][];         // the corresponding set of quals in canonical form (returned from healing)
            ReadState[][] readStatus = new ReadState[noReadFNs][];          // state of read after HealARead
            bool[][] readPossiblyChanged = new bool[noReadFNs][];           // read was possibly altered - use list form of quals just in case
            bool[][] readIsGood = new bool[noReadFNs][];                    // read found to be good (possibly after correction)  - and restored to its initial state
            bool[] allReadsGood = new bool[batchSize];                      // all reads in the 'pair' are good

            for (int f = 0; f < noReadFNs; f++)
            {
                fileActive[f] = true;                                       // stays true until EOF
                headerSet[f] = new Sequence[batchSize];
                readSet[f] = new Sequence[batchSize];
                qualHeaderSet[f] = new Sequence[batchSize];
                qualsSet[f] = new Sequence[batchSize];
                healedReadSet[f] = new Sequence[batchSize];
                healedQualSet[f] = new Sequence[batchSize];
                readStatus[f] = new ReadState[batchSize];
                readIsGood[f] = new bool[batchSize];
                readValid[f] = new bool[batchSize];

                for (int b = 0; b < batchSize; b++)
                {
                    headerSet[f][b] = new Sequence(defaultHeaderLength);
                    readSet[f][b] = new Sequence(defaultReadLength);
                    qualHeaderSet[f][b] = new Sequence(defaultHeaderLength);
                    qualsSet[f][b] = new Sequence(defaultReadLength);
                    healedReadSet[f][b] = new Sequence(defaultReadLength);
                    healedQualSet[f][b] = new Sequence(defaultReadLength);
                    readValid[f][b] = true;
                }
            }

            //Console.WriteLine(GC.GetTotalMemory(false) + " memory after HealingThread initialisation " + healerNumber);

            // get the next set of reads and correct them 
            int readsRead = 0;
            while (filesStillActive > 0)
            {
                threadWaitingRead[threadNumber - 1].Start();

                // read the next read from each of the files together - giving a pair of reads if we have paired reads - and do this for a batch to reduce locking overhead
                lock (bufferedReadsFiles)
                {
                    // fetch the next batch of reads 
                    timeSpentReading.Start();

                    for (int f = 0; f < noReadFNs; f++)
                    {
                        if (!fileActive[f])
                        {
                            for (int b = 0; b < batchSize; b++)
                                readValid[f][b] = false;
                            continue;
                        }

                        readsRead = bufferedReadsFiles[f].ReadReads(batchSize, headerSet[f], readSet[f], qualHeaderSet[f], qualsSet[f]);
                        if (readsRead != batchSize)
                        {
                            fileActive[f] = false;
                            for (int b = readsRead; b < batchSize; b++)
                                readValid[f][b] = false;
                            filesStillActive--;
                        }
                        threadStats.readsRead += readsRead;
                        progressReads += readsRead;
                    }

                    timeSpentReading.Stop();

                } // lock to ensure synchronised reading from paired reads files

                threadWaitingRead[threadNumber - 1].Stop();
                threadCorrecting[threadNumber - 1].Start();

                // process the just-acquired batch of reads
                for (int b = 0; b < readsRead; b++)
                {
                    // now have a set of reads (the n'th read from each file, if they exist. So try healing each one in turn. 
                    allReadsGood[b] = true;
                    for (int f = 0; f < noReadFNs; f++)
                        if (readValid[f][b])
                        {
                            if (wantedReadLength > 0)
                            {
                                int readLength = readSet[f][b].Length;
                                if (readLength > wantedReadLength)
                                {
                                    readSet[f][b].Length = wantedReadLength;
                                }
                                int qualsLength = qualsSet[f][b].Length;
                                if (qualsLength > wantedReadLength)
                                {
                                    qualsSet[f][b].Length = wantedReadLength;
                                }
                            }
                            if (readsFormat == SeqFiles.formatFASTQ)
                                SeqFiles.ConvertQualsToCanonicalForm(qualsSet[f][b], readsFormat, qualBase);      // FASTA quals were converted to numbers when they were read in

                            ReadState healingResult = HealARead(headerSet[f][b], readSet[f][b], qualsSet[f][b], readProps, healedReadSet[f][b], healedQualSet[f][b], threadStats);

                            readStatus[f][b] = healingResult;
                            readIsGood[f][b] = healingResult == ReadState.readOK || healingResult == ReadState.readCorrected;
                            allReadsGood[b] = allReadsGood[b] & readIsGood[f][b];
                        }
                }

                // have now healed this batch, so write out each set (unless one of them is 'bad' and we're only writing good read sets

                threadCorrecting[threadNumber - 1].Stop();
                threadWaitingWrite[threadNumber - 1].Start();

                // allocate write buffers 
                if (readsRead > 0)
                    for (int f = 0; f < noReadFNs; f++)
                    {
                        healedReadsWriteBuffers[f] = bufferedHealedReads[f].GetWriteBuffer();
                        if (bufferedHealedQuals[f] != null)
                            healedQualsWriteBuffers[f] = bufferedHealedQuals[f].GetWriteBuffer();
                        if (bufferedSingleReads[f] != null)
                            singleReadsWriteBuffers[f] = bufferedSingleReads[f].GetWriteBuffer();
                        if (bufferedSingleQuals[f] != null)
                            singleQualsWriteBuffers[f] = bufferedSingleQuals[f].GetWriteBuffer();
                        if (bufferedProblemReads[f] != null)
                            problemReadsWriteBuffers[f] = bufferedProblemReads[f].GetWriteBuffer();
                        if (bufferedProblemQuals[f] != null)
                            problemQualsWriteBuffers[f] = bufferedProblemQuals[f].GetWriteBuffer();
                    }

                // fill all the write buffers (and write any singles or problems)
                for (int b = 0; b < readsRead; b++)
                {
                    for (int f = 0; f < noReadFNs; f++)
                        if (readValid[f][b])
                        {
                            bool saveThisRead = pairedReads ? allReadsGood[b] : readIsGood[f][b];

                            // write this read to healed reads file if good enough
                            if (saveThisRead)
                            {
                                WriteReadAndQual(bufferedHealedReads[f], healedReadsWriteBuffers[f], bufferedHealedQuals[f], healedQualsWriteBuffers[f],
                                                 headerSet[f][b], healedReadSet[f][b], qualHeaderSet[f][b], healedQualSet[f][b]);

                                if (readStatus[f][b] == ReadState.readOK)
                                    threadStats.OKReadsWritten++;
                                if (readStatus[f][b] == ReadState.readCorrected)
                                    threadStats.correctedReadsWritten++;

                                continue;
                            }

                            // if we get to here, we may have a pair with one good read and one bad one 
                            if (readIsGood[f][b])
                            {
                                // this is an acceptable read so write it out to the appropriate singles file
                                WriteReadAndQual(bufferedSingleReads[f], singleReadsWriteBuffers[f], bufferedSingleQuals[f], singleQualsWriteBuffers[f],
                                                 headerSet[f][b], healedReadSet[f][b], qualHeaderSet[f][b], healedQualSet[f][b]);

                                threadStats.discardedOK++;
                            }
                            else
                            {
                                // this is a bad read
                                threadStats.discardedBroken++;

                                // and save it if we're saving problems
                                if (saveReads == SaveOption.saveProblemsAsWell)
                                {
                                    WriteReadAndQual(bufferedProblemReads[f], problemReadsWriteBuffers[f], bufferedProblemQuals[f], problemQualsWriteBuffers[f],
                                                     headerSet[f][b], readSet[f][b], qualHeaderSet[f][b], qualsSet[f][b]);
                                }
                            }
                        }
                }

                // release all the filled write buffers for writing - as paired sets as needed to keep read pairings
                if (readsRead > 0)
                {
                    BufferedWriter.WritePairedBufferSet(healedReadsWriteBuffers, healedQualsWriteBuffers);
                    BufferedWriter.WriteBufferSet(singleReadsWriteBuffers, singleQualsWriteBuffers);
                    BufferedWriter.WriteBufferSet(problemReadsWriteBuffers, problemQualsWriteBuffers);
                }

                threadWaitingWrite[threadNumber - 1].Stop();

                if (tracing > Tracing.traceOff && threadStats.readsRead > traceReadLimit)                  // only trace this much 
                    break;

            } // end of reading/healing loop

            //Console.WriteLine("read and healed " + threadStats.reads + " on thread " + healerNumber);
            //Console.WriteLine(freeReadContexts.Count + " pooled read contexts");

            lock (statsLock)
            {
                stats.readsRead += threadStats.readsRead;
                stats.OKReadsWritten += threadStats.OKReadsWritten;
                stats.correctedReadsWritten += threadStats.correctedReadsWritten;
                stats.brokenReadsFound += threadStats.brokenReadsFound;
                stats.shortReadsFound += threadStats.shortReadsFound;
                stats.discardedOK += threadStats.discardedOK;
                stats.OKReads += threadStats.OKReads;
                stats.OKReadsChecked += threadStats.OKReadsChecked;
                stats.healedReads += threadStats.healedReads;
                stats.readsNotHealedAtAll += threadStats.readsNotHealedAtAll;
                stats.healedFirstPass += threadStats.healedFirstPass;
                stats.healedAltPass += threadStats.healedAltPass;
                stats.healedRCPass += threadStats.healedRCPass;
                stats.hdubReads += threadStats.hdubReads;
                stats.tooDeepReads += threadStats.tooDeepReads;
                stats.tooShortReads += threadStats.tooShortReads;
                stats.abandonedReads += threadStats.abandonedReads;
                stats.abandonedRewriting += threadStats.abandonedRewriting;
                stats.abandonedTree += threadStats.abandonedTree;
                stats.abandonedNs += threadStats.abandonedNs;
                stats.abandonedFutile += threadStats.abandonedFutile;
                stats.mers += threadStats.mers;
                stats.replacedMers += threadStats.replacedMers;
                stats.discardedBroken += threadStats.discardedBroken;
                stats.trimmed += threadStats.trimmed;
                stats.extended += threadStats.extended;
                for (int i = 0; i < (int)FixTypes.maxFixTypes; i++)
                    stats.fixesByType[i] += threadStats.fixesByType[i];
#if STATS
                timeInChecking += threadInChecking.Elapsed.TotalSeconds;
                timeInTHRFirstPass += threadInTHRFirstPass.Elapsed.TotalSeconds;
                timeInTHRTHPass += threadInTHRTHPass.Elapsed.TotalSeconds;
                timeInTHRReversePass += threadInTHRReversePass.Elapsed.TotalSeconds;
#endif
            }
        }

        private static ReadState HealARead(Sequence readHeader, Sequence originalRead, Sequence originalQuals, ReadProperties readProps, Sequence read, Sequence quals, HealingStats threadStats)
        {
            bool readNeedsHealing = false;                      // does the read (now) need correction or can it pass through unscathed
            bool readWasChanged = false;                        // last attempt at healing changed something
            bool readIsFaulty = false;                          // read contains uncorrectable errors (or is too short after trimming)
            bool readWasCorrected = false;                      // read was corrected/trimmed at some time
            int minLengthForRead = 0;                           // reads that end up shorter than this are abandoned
            int maxChangedBases = 0;                            // as are reads that will have more than this many changes after correction/extension
            Tracing wasTracing = tracing;                       // remember previous tracing state so it can be reset (used when breaking on a read)

            bool initiallyOK = false;                           // |
            bool OKButChecked = false;                          // |
            bool healedFirstPass = false;                       // | used to defer stats updates until after length is determined to be OK
            bool healedSecondPass = false;                      // | 
            bool healedReversePass = false;                     // |

            Stopwatch healReadTimer = null;
            if (tracing > Tracing.traceOff)
            {
                healReadTimer = new Stopwatch();
                healReadTimer.Start();
            }

            const bool breakOnRead = false;            
            if (breakOnRead && (originalRead.Matches("TAAATCTAATTGGTATAATATAATATATCAAACTAATAATTTAATATTAATTATTAGTTCGATTGTAGTAGTAGGTTCTGTAACAATTTTAGGATATAAA") ||
                                readHeader.Matches("@SRR617721.7393796 SOLEXA4:20:D0V34ACXX:1:2114:5208:49282 length=101")))
            {
                // force tracing for this read if we're not already tracing (but changing trace state is not threadsafe...)
                wasTracing = tracing;
                if (tracing == Tracing.traceOff && noHealingThreads == 1)
                {
                    traces = new StreamWriter(tracesFN);
                    traces.WriteLine(myProcessNameAndArgs);
                    traces.WriteLine();
                    healReadTimer = new Stopwatch();
                    healReadTimer.Start();
                    tracing = Tracing.traceChanges;
                }
                // and change the tracing level as desired
                if (tracing != Tracing.traceOff)
                {
                    tracing = Tracing.traceChoices;
                    //tracing = Tracing.traceRead;
                    tracing = Tracing.traceFollowers;
                    tracingThisRead = true;
                }
#if (DEBUG)
                if (breakOnRead)
                    Debugger.Break();
#endif
            }

            // make a copy of the read (this is the one that will be modified). The original read is kept in case we need to revert back to it.
            originalRead.CopyTo(read);
            originalQuals.CopyTo(quals);

            int startOfPoorQualTail = FindPoorQualTail(quals, minOKQual);
            minLengthForRead = read.Length * saveGoodMinLengthPct / 100;
            maxChangedBases = originalRead.Length - minLengthForRead;
            threadStats.mers += (read.Length - merSize + 1);

            if (minLengthForRead < merSize)
                minLengthForRead = merSize;
            if (read.Length < merSize)                          // ignore reads shorter than a single k-mer
            {
                threadStats.shortReadsFound++;
                return ReadState.readNotLongEnough;             // and exit early to avoid problems coming from having no kMers at all
            }

#if STATS
            threadInChecking.Start();
#endif
            // initialise the status for this new read 
            readProps.Initialise();
            // get the depths of the reads we're about to heal. The depths & kMers will be updated by TryHealingRead, but the pairs will not.
            SetMerDepths(read, readProps);
            SetPairDepths(read, readProps);
            // and calculate the target depths for this read 
            SetDepthThreshholds(readProps, readProps.merCount);
            SetPairDepthThreshholds(readProps, readProps.merCount);
            // and see if needs healing
            readProps.initialReadState = CheckReadForProblems(read, readProps);
            // and set the marker for the noisy tail (where we don't try so hard when counting followers)
            readProps.startOfNoisyTail = startOfPoorQualTail >= merSize ? startOfPoorQualTail - merSize + 1 : 0;

#if STATS
            threadInChecking.Stop();
#endif

            //lock (depthsTrace)
            //{
            //    depthsTrace.WriteLine(readHeader.ToString() + "\t" + readProps.OKDepth + "\t" + readProps.minDepth + "\t" + readProps.OKPairDepth + "\t" +
            //                          readProps.minPairDepth + "\t" + readProps.hmZeroPresent + "\t" + readProps.unbalancedRead); 
            //}

            // should we try to correct this read?
            readNeedsHealing = readProps.initialReadState == ReadState.readIsBroken || readProps.initialReadState == ReadState.readNeedsChecking;

            //readNeedsHealing = false;
            //readProps.initialReadState = ReadState.readOK;

            if (tracing > Tracing.traceOff && readProps.initialReadState == ReadState.readOK)
            {
                TraceHeader(readHeader, readProps);
                TraceRead("*", read);
            }

            if (tracing > Tracing.traceOff && readProps.initialReadState == ReadState.readTooDeep)
            {
                TraceHeader(readHeader, readProps);
                TraceRead("!", read);
            }

            // does the read need to have adapter sequences trimmed?
            if (readNeedsHealing && readProps.deepUnbalancedPresent)
            {
                int newStartingBase;
                int newLength;
                int basesTrimmed;

                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("Trimming HDUB");

                basesTrimmed = TrimExtend.TrimUnbalancedHighDepth(kMersTable, read, quals, readProps.minDepth, readProps.mers, out newStartingBase, out newLength);

                if (basesTrimmed > 0)
                {
                    // trim the read/qual if we found any adapter kMers
                    read.Trim(newStartingBase, newLength);
                    if (quals.Length > 0)
                        quals.Trim(newStartingBase, newLength);

                    if (read.Length < merSize)                          // ignore reads now shorter than a single k-mer
                    {
                        threadStats.shortReadsFound++;
                        return ReadState.readNotLongEnough;
                    }

                    if (tracing > Tracing.traceOff)
                    {
                        TraceClumpAdd("Trimmed " + basesTrimmed + " HDUB bases");
                        if (readHeader != null && readHeader.Length != 0)
                            TraceClumpAdd(readHeader.ToString());
                        TraceRead("-", read);
                    }

                    // recalculate depths etc for the trimmed read
                    readProps.Initialise();
                    SetMerDepths(read, readProps);
                    SetPairDepths(read, readProps);
                    SetDepthThreshholds(readProps, readProps.merCount);
                    SetPairDepthThreshholds(readProps, readProps.merCount);
                    readProps.initialReadState = CheckReadForProblems(read, readProps);

                    // and remember for the stats
                    threadStats.hdubReads++;
                }
            }

            // reads that did not need checking/healing at all (perhaps state changed after trimming)
            if (readProps.initialReadState == ReadState.readOK)
            {
                readNeedsHealing = false;
                initiallyOK = true;
                readProps.finalReadState = ReadState.readOK;
                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("Read OK");
            }

            // read look like it could need correction, so start a full left-to-right healing pass
            if (readNeedsHealing)
            {
#if STATS
                threadInTHRFirstPass.Start();
#endif
                // try healing the read in a forward direction. The readProperties are set before the HDUB test and updated by TryHealingRead
                readWasChanged = TryHealingRead(readHeader, read, quals, readProps, 0, threadStats);
#if STATS
                threadInTHRFirstPass.Stop();
#endif
                readProps.finalReadState = (readProps.remainingBadMers > 0 || readProps.healingAbandoned) ? ReadState.readIsBroken : ReadState.readOK;

                // if we had to recalculate the depth thresholds, it's possible that the read before the recalc needs further fixes
                if (readProps.depthsRecalculated && readProps.finalReadState == ReadState.readOK)
                {
                    // check read again with the updated thresholds
                    SetPairDepths(read, readProps);     // pairs were not maintained during correction
                    SetPairDepthThreshholds(readProps, readProps.merCount);
                    readProps.finalReadState = CheckReadForProblems(read, readProps);
                    // but only concerned with errors... not putative correction opportunities
                    if (readProps.finalReadState == ReadState.readNeedsChecking)
                        readProps.finalReadState = ReadState.readOK;
                    if (tracing > Tracing.traceOff)
                        TraceClumpAdd("Read OK after correction, but checking as depths were adjusted");
                }

                if (readProps.finalReadState == ReadState.readOK)
                {
                    readNeedsHealing = false;
                    if (readWasChanged)
                    {
                        // broken reads that were completely healed in the forward pass
                        healedFirstPass = true;
                        readWasCorrected = true;

                        if (tracing > Tracing.traceOff)
                            TraceClumpAdd("Read healed (forward)");
                    }
                    else
                    {
                        // reads that we tried to heal but turned out OK after all
                        OKButChecked = true;

                        if (tracing > Tracing.traceOff)
                            TraceClumpAdd("Read OK (but checked over)");
                    }
                }
            }

            // Illumina sometimes produces unbalanced but good reads, often with low coverage regions
            // These reads will be unbalanced, and no better alternatives will have been found for any of the kMers.
            // Rescue any such reads that have consistent kMer or pair coverage.
            if (!readWasChanged && readProps.finalReadState == ReadState.readIsBroken && readProps.unbalancedRead)
            {
                bool incompleteCoverage = false;

                // if we have pairs, check that the covered-by-pairs region of the read is covered by either good-enough kMers or pairs
                if (readProps.pairsCount > 0)
                {
                    // pair depths aren't maintained so recreate them if the read has changed
                    if (readWasChanged)
                    {
                        SetPairDepths(read, readProps);
                        SetPairDepthThreshholds(readProps, readProps.merCount);
                    }

                    for (int i = 0; i < readProps.pairsCount; i++)
                        if (readProps.depths[i] == 0 || readProps.pairDepths[i] == 0)
                        {
                            incompleteCoverage = true;
                            break;
                        }
                }

                // and check that the end of the read (not covered in the above pairs+kMer depths check) also has good enough depth
                if (!incompleteCoverage)
                {
                    for (int i = readProps.pairsCount; i < readProps.merCount; i++)
                        if (readProps.depths[i] == 0)
                        {
                            incompleteCoverage = true;
                            break;
                        }
                }

                if (!incompleteCoverage)
                {
                    readProps.finalReadState = ReadState.readOK;
                    OKButChecked = true;
                    readNeedsHealing = false;
                    if (tracing > Tracing.traceOff)
                        TraceClumpAdd("Pass-through of unchanged, fully-covered unbalanced read");
                }
            }

            // This attempt at healing wasn't (completely) successful. If the start of the read is bad, but it then comes good,
            // we'll reverse the corrected read and try correcting the errant region from the good section. This will allow us to correct more bad
            // bases at the start of the read than could be done in the normal left-to-right scan. Don't do this if the previous correction attempt
            // was abandoned as this indicates there may an uncorrected failure somewhere in the middle of the read which will have resulted in broken pairs.
            if (readNeedsHealing && readProps.firstGoodMer > 0 && StartOfReadIsBad(readProps) &&
                     !(readProps.healingAbandoned && (readProps.abandonReason == AbandonReason.tooManyNs || readProps.abandonReason == AbandonReason.treeSize)))
            {
#if STATS
                threadInTHRReversePass.Start();
#endif
                if (readProps.healingAbandoned)
                {
                    originalRead.CopyTo(read);
                    originalQuals.CopyTo(quals);
                }

                Sequence reversedRead = AllocateSequence();
                read.CopyTo(reversedRead);
                Sequence.ReverseComplement(reversedRead);
                Sequence reversedQuals = AllocateSequence();
                quals.CopyTo(reversedQuals);
                if (reversedQuals.Length > 0)
                    reversedQuals.Reverse();

                int changedMersForward = readProps.changedMers;
                readProps.Reset();
                // get the depths of the reversed read
                SetMerDepths(reversedRead, readProps);
                // thresholds will be unchanged

                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("Trying reverse direction");

                //int startingPointForScan = read.Length - (readProps.firstGoodMer - 1) - merSize;

                readWasChanged = TryHealingRead(readHeader, reversedRead, reversedQuals, readProps, 0, threadStats);

                if (readWasChanged)
                {
                    // restore read direction
                    Sequence.ReverseComplement(reversedRead);
                    if (reversedQuals.Length > 0)
                        reversedQuals.Reverse();

                    // copy reversed-reversed read back to healed read and restore depths etc
                    reversedRead.CopyTo(read);
                    if (reversedQuals.Length > 0)
                        reversedQuals.CopyTo(quals);
                    SetMerDepths(read, readProps);
                }
                else
                {
                    originalRead.CopyTo(read);
                    originalQuals.CopyTo(quals);
                    SetMerDepths(read, readProps);
                }

                // can't do abandon trimmimg after a reverse pass - abandonAtM will be from the reversed read, and the trimming will be on for forward read
                readProps.abandonedAtM = -1;

                ReturnSequence(reversedRead);
                ReturnSequence(reversedQuals);

                readProps.changedMers += changedMersForward;
                readProps.finalReadState = (readProps.remainingBadMers > 0 || readProps.healingAbandoned) ? ReadState.readIsBroken : ReadState.readOK;

                if (readWasChanged && readProps.finalReadState == ReadState.readOK)
                {
                    // broken reads that were completely healed after trying in the reverse direction
                    readNeedsHealing = false;
                    healedReversePass = true;
                    readWasCorrected = true;

                    if (tracing > Tracing.traceOff)
                    {
                        TraceRead("<", read);
                        TraceClumpAdd("Read healed (reverse)");
                    }
                }

#if STATS
                threadInTHRReversePass.Stop();
#endif
            } // third (reverse) pass

            // faulty kMers may still be present in the read so trim it back to good kMers from both ends
            if (readNeedsHealing)
            {
                // when we get here, the kMers and depths will be set correctly but pair depths may be out-of-date
                if (readWasChanged)
                {
                    SetPairDepths(read, readProps);
                    SetPairDepthThreshholds(readProps, readProps.merCount);
                }

                if (readProps.healingAbandoned && readProps.abandonedAtM > 0)
                {
                    int lengthBeforeTrimming = read.Length;
                    TrimAbandonedRead(read, quals, readProps);
                    if (tracing > Tracing.traceOff)
                    {
                        TraceClumpAdd("trimmed end of abandoned read by " + (lengthBeforeTrimming - read.Length) + " bases");
                        TraceRead("-", read);
                    }
                }

                int basesTrimmedFromStart = 0;
                if (read.Length > merSize)
                {
                    if (readProps.depths[0] < readProps.minDepth || (readProps.pairsCount > 0 && readProps.pairDepths[0] < readProps.minPairDepth))
                        basesTrimmedFromStart = TrimExtend.TrimReadStart(read, quals, readProps.merCount, readProps.pairsCount, false, readProps.depths, readProps.minDepth, readProps.pairDepths, readProps.minPairDepth, 0, read.Length);
                }
                int lengthAfterTrimming = read.Length - basesTrimmedFromStart;
                if (tracing > Tracing.traceOff && basesTrimmedFromStart > 0)
                {
                    TraceClumpAdd("trimmed start by " + basesTrimmedFromStart + " bases");
                    if (lengthAfterTrimming > 0)
                        TraceRead("-", read.SubSeq(basesTrimmedFromStart, lengthAfterTrimming));
                }

                int basesTrimmedFromEnd = 0;
                if (lengthAfterTrimming > merSize)
                {
                    basesTrimmedFromEnd = TrimExtend.TrimReadEnd(read, quals, readProps.merCount, readProps.pairsCount, false,
                                                                 readProps.mers, readProps.depths, readProps.minDepth, readProps.pairDepths, readProps.minPairDepth,
                                                                 basesTrimmedFromStart, lengthAfterTrimming, merSize, kMersTable);
                    lengthAfterTrimming = lengthAfterTrimming - basesTrimmedFromEnd;
                    if (tracing > Tracing.traceOff && basesTrimmedFromEnd > 0)
                    {
                        TraceClumpAdd("trimmed end by " + basesTrimmedFromEnd + " bases");
                        if (lengthAfterTrimming > 0)
                            TraceRead("-", read.SubSeq(basesTrimmedFromStart, lengthAfterTrimming));
                    }
                }

                // finished trimming so (in-place) Trim the read to remove any unwanted ends
                int totalBasesTrimmed = basesTrimmedFromStart + basesTrimmedFromEnd;
                if (totalBasesTrimmed > 0)
                {
                    read.Trim(basesTrimmedFromStart, lengthAfterTrimming);
                    if (quals.Length > 0)
                        quals.Trim(basesTrimmedFromStart, lengthAfterTrimming);

                    // trim kMer & pair depths as well
                    readProps.merCount -= totalBasesTrimmed;
                    for (int i = 0; i < readProps.merCount; i++)
                    {
                        readProps.mers[i] = readProps.mers[basesTrimmedFromStart + i];
                        readProps.depths[i] = readProps.depths[basesTrimmedFromStart + i];
                    }
                    if (readProps.pairsCount > 0)
                    {
                        readProps.pairsCount -= totalBasesTrimmed;
                        for (int i = 0; i < readProps.pairsCount; i++)
                            readProps.pairDepths[i] = readProps.pairDepths[basesTrimmedFromStart + i];
                    }
                }

                // read has been trimmed either here or when it was abandoned during healing
                if (originalRead.Length - read.Length > 0)
                {
                    threadStats.trimmed++;
                    readProps.finalReadState = CheckReadForProblems(read, readProps);
                    if (readProps.finalReadState == ReadState.readNeedsChecking)
                        readProps.finalReadState = ReadState.readOK;
                    if (readProps.finalReadState == ReadState.readOK)
                        readWasCorrected = true;
                }

            } // end of trimming of potentially broken reads

            // extend reads if requested (and we're not skipping the read because it's too deep to be of interest, and we don't think it's broken)
            if (readProps.initialReadState != ReadState.readTooDeep)
            {
                int basesToAdd = 0;
                // if we're trying to restore the original read length, see how many bases we'd have to add to get back to there
                if (readsFixedLength && read.Length < originalRead.Length)
                    basesToAdd = originalRead.Length - read.Length;
                // asked to extend read by N bases, so go for this unless we're already trying to extend for more than this
                if (extendBases > basesToAdd)
                    basesToAdd = extendBases;

                // only extend reads that meet the goodness criteria (changed bases + added bases < allowed altered/added bases)
                if (basesToAdd + readProps.changedMers > maxChangedBases)
                {
                    basesToAdd = 0;
                    if (tracing > Tracing.traceOff)
                        TraceClumpAdd("read too short/altered to extend: " + read.Length + "/" + originalRead.Length);
                }

                if (basesToAdd > 0)
                {
                    if (readProps.finalReadState != ReadState.readIsBroken)
                    {
                        int basesAdded = TrimExtend.ExtendRead(read, quals, readProps.merCount, readProps.pairsCount, readProps.mers,
                                                                readProps.depths, readProps.minDepth, readProps.pairDepths, readProps.minPairDepth,
                                                                kMersTable, pairsTable, basesToAdd);
                        if (basesAdded > 0)
                        {
                            threadStats.extended++;
                            if (tracing > Tracing.traceOff)
                            {
                                TraceClumpAdd("extended by " + basesAdded + "/" + basesToAdd + " bases");
                                TraceRead("+", read);
                            }
                        }
                    }
                    else
                    {
                        if (tracing > Tracing.traceOff)
                            TraceClumpAdd("read broken - not extended");
                    }
                }
            }

            // read too short (possibly after trimming)
            if (read.Length < minLengthForRead && readProps.finalReadState != ReadState.readIsBroken)
            {
                threadStats.tooShortReads++;
                readProps.finalReadState = ReadState.readNotLongEnough;
            }

            // too many changed kMers?
            if (readProps.changedMers > maxChangedBases)
            {
                threadStats.tooShortReads++;
                readProps.finalReadState = ReadState.readNotLongEnough;
            }

            // if the corrected/trimmed read is shorter than allowed, or still broken, restore read to its starting state and declare it abandoned
            if (readProps.finalReadState == ReadState.readIsBroken || readProps.finalReadState == ReadState.readNotLongEnough)
            {
                readIsFaulty = true;
                originalRead.CopyTo(read);
                originalQuals.CopyTo(quals);

                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("uncorrected read: " + TraceBrokenReason(readProps));
            }

            if (readProps.finalReadState == ReadState.readOK)
            {
                if (initiallyOK)
                {
                    threadStats.OKReads++;
                    progressOKReads++;
                }
                if (OKButChecked)
                {
                    threadStats.OKReadsChecked++;
                    progressOKReads++;
                }
                if (readWasCorrected)
                {
                    threadStats.healedReads++;
                    progressHealedReads++;
                    if (healedFirstPass)
                        threadStats.healedFirstPass++;
                    if (healedSecondPass)
                        threadStats.healedAltPass++;
                    if (healedReversePass)
                        threadStats.healedRCPass++;
                }
            }

            // add trailing Ns to fixed length reads if requested (unless the read is abandoned)
            if (readsFixedLengthPadded && read.Length < originalRead.Length && read.Length >= merSize && !readIsFaulty)
            {
                int NsToAdd = originalRead.Length - read.Length;
                for (int i = 0; i < NsToAdd; i++)
                {
                    read.Append('N');
                    if (quals.Length > 0)
                        quals.Append((char)1);
                }
            }

            if (readProps.finalReadState == ReadState.readIsBroken && !readWasChanged)
            {
                threadStats.readsNotHealedAtAll++;
                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("Broken but no corrections made");
            }

            if (readProps.initialReadState == ReadState.readTooDeep)
            {
                readIsFaulty = true;
                threadStats.tooDeepReads++;
                if (tracing > Tracing.traceOff)
                    TraceClumpAdd("Read too deep to correct");
            }

            if (tracing > Tracing.traceOff)
            {
                healReadTimer.Stop();
                float timeTakenToHeal = (float)(healReadTimer.ElapsedMilliseconds) / 1000;
                bool slowHealingRead = timeTakenToHeal > 1.0;
                TraceClumpAdd("healing took " + timeTakenToHeal.ToString("F3") + (slowHealingRead ? " (slow)" : ""));
                WriteTraceClump(traceClump);
                traces.WriteLine();
            }

            if (tracing > Tracing.traceOff)                       // reset one-off trace-in-depth toggle  
            {
                tracing = wasTracing;
                tracingThisRead = false;
                if (wasTracing == Tracing.traceOff && traces != null)
                    traces.Close();
            }

            //if (read.Length != quals.Length)
            //    Debugger.Break();

            if (readProps.finalReadState == ReadState.readOK && readWasCorrected)
                readProps.finalReadState = ReadState.readCorrected;
            return readProps.finalReadState;

        }  // for a single read 

        private static void WriteReadAndQual(BufferedWriter bufferedReadsWriter, BufferedWriter.WriteBuffer readsBuffer, BufferedWriter bufferedQualsWriter, BufferedWriter.WriteBuffer qualsBuffer, Sequence readHeader, Sequence read, Sequence qualHeader, Sequence quals)
        {
            WriteRead(bufferedReadsWriter, readsBuffer, readHeader, read);
            if (readsFormat == SeqFiles.formatFASTQ)
                WriteQual(bufferedReadsWriter, readsBuffer, qualHeader, quals);
            if (readsFormat == SeqFiles.formatFNA && quals.Length > 0)
                WriteQual(bufferedQualsWriter, qualsBuffer, qualHeader, quals);
        }

        private static void WriteRead(BufferedWriter bufferedWriter, BufferedWriter.WriteBuffer buffer, Sequence readHeader, Sequence readSeq)
        {
            // @1:1:0:686#0/1
            // NTGGAGAATTCTGGATCCTCGGACTAAAACAATAGCAGTTGATTCGCTCACAGTTCTGGAGGCTAGAGGTATGAAA
            // +1:1:0:686#0/1
            // @hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\hhhhhhhhhhhhhhhhVg[hhhU^hhWhfgVhc^Dhh`_V
            //
            // >FUOCO2J02IPX2K rank=0000260 x=3458.0 y=3626.0 length=70
            // GAGCAGCTCCTCAAGCAACTGCAACTGGATTGAGCAAGGAATGTTCGCAGCTACCCGACT
            // GACCCGTCTT

            // write out the header if it exists
            if (readHeader != null)
                bufferedWriter.WriteLine(buffer, readHeader);

            if (readsFormat == SeqFiles.formatFNA)
            {
                // if we concatenated shorter source lines (e.g. from 454 reads) when the data was read,
                // we'll now split it again as it's being written
                int m = 0;
                int hrLen = readSeq.Length;
                while (m < hrLen)
                {
                    int wLen = Math.Min(60, hrLen - m);
                    bufferedWriter.WriteLine(buffer, readSeq, m, wLen);
                    m += 60;
                }
            }
            else
            {
                bufferedWriter.WriteLine(buffer, readSeq);
            }
        }

        private static void WriteQual(BufferedWriter bufferedWriter, BufferedWriter.WriteBuffer buffer, Sequence qualHeader, Sequence healedQuals)
        {
            // The quals will come in text form (healedQualsString) if they are unchanged.
            // Altered quals will always come in List form. 

            // no current qual, but there were some at some stage... happens with regression test 454 data but should otherwise never happen
            if (healedQuals.Length == 0)
                return;

            if (readsFormat == SeqFiles.formatFNA)
            {
                if (qualHeader.Length != 0)
                    bufferedWriter.WriteLine(buffer, qualHeader);

                // easier to just convert to list and then write out 60 quals/line than to parse the concatenated qual string 
                // if there were known to be no changes to the text-form quals
                List<int> healedQualsList = new List<int>(500);
                for (int i = 0; i < healedQuals.Length; i++)
                    healedQualsList.Add((int)healedQuals.Bases[i]);

                // if we concatenated shorter qual lines (e.g. from 454 reads) when they were being read,
                // we'll now reformat and split them again as they're being written
                int m = 0;
                int qualCount = healedQualsList.Count;
                while (m < qualCount)
                {
                    int wLen = Math.Min(60, qualCount - m);
                    for (int i = m; i < m + wLen; i++)
                    {
                        int nextQual = healedQualsList[i];
                        bufferedWriter.Write(buffer, nextQual.ToString() + ' ');
                    }
                    bufferedWriter.WriteLine(buffer);
                    m += 60;
                }
            }

            if (readsFormat == SeqFiles.formatFASTQ)
            {
                bufferedWriter.WriteLine(buffer, qualHeader);

                // convert the canonical quals back to either Sanger or Solexa format
                for (int i = 0; i < healedQuals.Length; i++)
                    healedQuals.Bases[i] += (char)qualBase;

                bufferedWriter.WriteLine(buffer, healedQuals);
            }
        }

        private static int FindPoorQualTail(Sequence quals, int minOKQual)
        {
            if (quals.Length == 0)
                return -1;

            const int windowLength = 10;
            const int passesNeededInWindow = windowLength * 3 / 4;
            int basesTrimmed = 0;
            char[] qualCharArray = quals.Bases;
            int goodBaseIdx = 0;
            //int[] qualInts = new int[qualCharArray.Length];
            //for (int i = 0; i < qualCharArray.Length; i++)
            //    qualInts[i] = qualCharArray[i];

            for (int i = quals.Length - 1; i >= 0; i--)
                if (qualCharArray[i] >= minOKQual)
                {
                    goodBaseIdx = i;
                    break;
                }

            basesTrimmed = quals.Length - goodBaseIdx;

            for (int i = goodBaseIdx; i > windowLength; i--)
            {
                int passedInWindow = 0;
                for (int w = 0; w < windowLength; w++)
                {
                    int qualInt = (int)qualCharArray[i - w];
                    if (qualInt >= minOKQual)
                        passedInWindow++;
                }
                if (passedInWindow < passesNeededInWindow)
                    basesTrimmed++;
                else
                    break;
            }

            return quals.Length - basesTrimmed;
        }

        private static void WriteTraceClump(List<string> lines)
        {
            if (lines.Count == 0)
                return;

            lock (traceLock)
            {
                foreach (string l in lines)
                    traces.WriteLine(l);
                traces.Flush();
            }
            lines.Clear();
        }

        private static void TraceClumpAdd(string s)
        {
            lock (traceLock)
            {
                traceClump.Add(s);
            }
        }

        private static void TraceHeader(Sequence readHeader, ReadProperties readProps)
        {
            if (readHeader != null && readHeader.Length != 0)
                TraceClumpAdd(readHeader.ToString() + " OK=" + readProps.OKDepth + " min=" + readProps.minDepth + " pOK=" + readProps.OKPairDepth +
                                " pmd=" + readProps.minPairDepth + " ubr=" + readProps.unbalancedRead + " hmz=" + readProps.hmZeroPresent + " snt=" + readProps.startOfNoisyTail);
        }

        private static void TraceRead(string tag, Sequence read)
        {
            string traceLine = "";
            int merCount = 0;
            int[] merPlusDepths;
            int[] merRCDepths;

            traceLine += tag + new string(read.Bases, 0, read.Length) + " (" + read.Length + ")";
            TraceClumpAdd(traceLine);
            if (read.Length > 0)
            {
                traceLine = "";
                GetMerDepthsForTrace(read, out merPlusDepths, out merRCDepths);
                merCount = merPlusDepths.Length;

                for (int m = 0; m < merCount; m++)
                    traceLine += (merPlusDepths[m] + "/" + merRCDepths[m] + "\t");
                TraceClumpAdd(traceLine);
            }
        }

        private static string Indent(int n)
        {
            return new string(' ', n * 5);
        }

        private static void GetMerDepthsForTrace(Sequence read, out int[] merPlusDepths, out int[] merRCDepths)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            merPlusDepths = new int[mersInRead];
            merRCDepths = new int[mersInRead];

            // read too short to tile for mers
            if (readLength < merSize)
                return;

            bool merIsValid = false;
            ulong lastMer = 0;

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                    merIsValid = Sequence.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = Sequence.CondenseMer(read, i, merSize, out lastMer);
                if (merIsValid)
                {
                    int plusCount = 0;
                    int rcCount = 0;
                    if (kMersTable.GetDepthCounts(lastMer, out plusCount, out rcCount))
                    {
                        merPlusDepths[i] = plusCount;
                        merRCDepths[i] = rcCount;
                    }
                }
                // depth arrays are zero by default so don't bother setting zero counts in all other cases
            }
        }

        // Makes a single pass through a read, trying to heal it whenever it finds that the mer depth
        // has fallen below the 'minReps' target or the +/RC values are unbalanced or whenever else it seems sensible.
        //
        // Returns true if the read was corrected at all - and updates various bits of read status in readProps
        //
        // If the healing is being done in reverse, the string has already been RCed and only the start (now end) of the read is corrected
        // 
        private static bool TryHealingRead(Sequence readHeader, Sequence read, Sequence quals, ReadProperties readProps, int startingPointForScan, HealingStats threadStats)
        {
            Stopwatch healReadTimer = null;
            if (perfTrace)
            {
                healReadTimer = new Stopwatch();
                healReadTimer.Start();
            }

            ulong kMer;                                             // the current kMer being examined/corrected
            string referenceRead = null;                            // reference read (for tracing)
            int minDepth = readProps.minDepth;                      // depths lower than this are generally regarded as 'bad'
            int OKDepth = readProps.OKDepth;                        // and lower than this are 'poor'
            int minPairDepth = readProps.minPairDepth;              // and pairs (if present) must be at least this deep
            int OKPairDepth = readProps.OKPairDepth;                // and good ones are this deep
            kMerState checkReason = kMerState.merOK;                // why are we looking for a replacement for the current kMer?
            ReadState readAfterHealing = ReadState.readUnknown;
            bool readWasChanged = false;                            // assume that the read is not altered (until we alter it)
            int firstGoodMer = -1;                                  // index of first 'good' kMer
            int changedMers = 0;                                    // how many kMers were changed/fixed

            int initialMerCount = readProps.merCount;               // initial read length in kMers (used to stop reads growing)
            int initialReadLength = read.Length;                    // and starting length in bases (used to trim uncorrected ends)
            int tailOfRead = Math.Max((initialMerCount - 10), (initialMerCount * 90) / 100);   // tail of read is last 10% or 10
            bool subFixesOnly = onlySubsAllowed;                    // set whenever it seems sensible to restrict fix types for performance reasons (e.g. high depth reads)

            int currentMerCount = initialMerCount;                  // how many mers are in the current being-healed read
            int previousMerDepth = -1;                              // depth of the previous mer examined/fixed in the scanning loop
            int previousPairDepth = -1;                             // depth of the previous pair
            int consecutiveNs = 0;                                  // count of consecutive N fixes - stop rewriting long N regions
            int nextMerToCheck = 0;                                 // next kMer (idx) to check

            firstGoodMer = -1;                                      // no good mers to start
            int thmCalls = 0;                                       // start counting again for each THR call

            // trace read before we start healing it
            if (tracing > Tracing.traceOff)
            {
                if (readHeader != null && readHeader.Length != 0)
                {
                    string faHeaderString = readHeader.ToString();
                    TraceHeader(readHeader, readProps);
                    int readIDEnd = faHeaderString.IndexOf(' ');
                    if (readIDEnd < 0)
                        readIDEnd = faHeaderString.Length;
                    string readID = faHeaderString.Substring(1, readIDEnd - 1);
                    if (refSequences.ContainsKey(readID))
                    {
                        referenceRead = refSequences[readID];
                        TraceClumpAdd(referenceRead);
                    }
                }
                TraceRead(">", read);
            }

            readProps.remainingBadMers = 0;                         // start again for each pass over the read

            // something may be broken in this read so move along it progressively, trying to fix each mer and 
            // rebuilding the read string as we make fixes. 

            for (int m = startingPointForScan; m < currentMerCount; m++)
            {
                int mersToFirstChoice = 1;                              // how many (checked) bases can be skipped after a check/fix
                FixTypes fixType = FixTypes.fixNone;                    // assume no change is made to this mer
                char followingBase = ' ';                               // base following this mer

                thmCalls = 0;                                           // restart the thmCalls counter afresh

                if (m < currentMerCount - 1)
                    followingBase = read.Bases[m + merSize];
                else
                    followingBase = '\0';

                if (!Sequence.CondenseMer(read, m, merSize, out kMer))             // next mer must have contained Ns if this failed
                {
                    string merWithNs = null;                                    // saved mer including Ns (for tracing only)

                    if (FindBestReplacementForNs(read, m, OKDepth, out kMer))   // so try to find the best possible replacement
                    {
                        if (tracing > Tracing.traceOff)
                            merWithNs = read.ToString(m, merSize);
                        string replacedNsMer = kMers.ExpandMer(kMer, merSize);
                        for (int r = 0; r < merSize; r++)                       // and replace the N-containing mer
                            read.Bases[m + r] = replacedNsMer[r];
                        readWasChanged = true;
                        readProps.merChangeCost[m] = 1;
                        fixType = FixTypes.fixN;
                        consecutiveNs++;
                        // are we just rewriting too much?
                        if (consecutiveNs > maxConsecutiveNFixesAllowed && m < tailOfRead)
                        {
                            threadStats.abandonedNs++;
                            if (tracing > Tracing.traceOff)
                                TraceClumpAdd("read abandoned - rewriting Ns");

                            readProps.healingAbandoned = true;
                            readProps.abandonReason = AbandonReason.tooManyNs;
                            readProps.abandonedAtM = m > 0 ? m - 1 : 0;
                            break;                                      // and give up on it
                        }

                        if (tracing > Tracing.traceOff)
                        {
                            string traceLine = "";
                            traceLine += fixTypes[(int)fixType];
                            for (int p = 0; p < m; p++)
                                traceLine += " ";
                            traceLine += kMers.ExpandMer(kMer, merSize) + " for " + merWithNs + " @" + m;
                            TraceClumpAdd(traceLine);
                        }
                    }
                    else
                    {
                        // can't find a suitable replacement mer, so move on - and it could be fixed in a later reverse pass
                        // (could be a combination of Ns and errors)
                        readProps.remainingBadMers++;
                        continue;
                    }
                }
                else
                    consecutiveNs = 0;

                if (tracing >= Tracing.traceRead)
                {
                    bool breakOnMer = false;
                    string traceMer = kMers.ExpandMer(kMer, merSize);
                    if (breakOnMer && traceMer == "TCCGGCCAAGGGCCAGGCCAAGGCA")
                    {
                        TraceClumpAdd("Found target mer " + traceMer + " @" + m);
                        Debugger.Break();
                    }
                }

                // look up next mer and retrieve counts
                bool unbalancedMer;
                bool tiltedMer;
                int merDepth = kMersTable.GetDepthSum(kMer, minDepth, out unbalancedMer, out tiltedMer);
                readProps.mers[m] = kMer;
                readProps.depths[m] = merDepth;
                int pairDepth = GetPairDepth(null, kMer, readProps.mers, readProps.depths, m, minDepth);
#if STATS
                Interlocked.Increment(ref gds4Calls_THR);
                Interlocked.Increment(ref gpdCalls);
#endif
                ulong startingMer = kMer;
                checkReason = kMerState.merOK;

                // see if this mer is a candidate for healing -- either broken or worth checking out
                if (m == nextMerToCheck && MerIsHealingCandidate(read, readProps, m, kMer, unbalancedMer, tiltedMer, followingBase, merDepth, previousMerDepth, pairDepth, previousPairDepth, true, out checkReason))
                {
                    kMerProperties replacementMer = null;                   // status etc from TryHealingMer call
                    bool merWasChanged = false;                             // did TryHealingMer make a change?
                    int costOfChange = 0;                                   // cost of the change (if the kMer was changed)
                    int followerRepairs = 0;                                // starting value for the number of follower repairs allowed (decremented on recursive calls)

                    // allow repairs to go to the end of the read if we're in the tail
                    if (m > tailOfRead)
                    {
                        if (errorModel == ErrorModel.mostlySubs)
                            subFixesOnly = true;
                        followerRepairs = currentMerCount - m;
                    }

                    if (followerRepairs == 0)
                    {
                        // still in a good part of the read so allow full recursion depth
                        if (m < readProps.startOfNoisyTail)
                            followerRepairs = maxFollowerRepairs;
                        else
                        // reduce allowed recursive depth if we're in the noisy tail (and limit fixes to subs when appropriate)
                        {
                            followerRepairs = maxFollowerRepairsTail;
                            if (errorModel == ErrorModel.mostlySubs)
                                subFixesOnly = true;
                        }
                    }

                    bool thmOK = TryHealingMer(checkReason, subFixesOnly, false, kMer, merDepth, pairDepth, previousMerDepth, previousPairDepth, unbalancedMer,
                                                   read, readProps, m, currentMerCount, followerRepairs, FixTypes.fixNone, 0, out replacementMer, true, 0, ref thmCalls);
                    
                    if (!thmOK)
                    {
                        // must be tree explosion error as that's the only way THM can fail
                        if (tracing > Tracing.traceOff)
                            TraceClumpAdd("tree explosion - read abandoned @" + m);
                        threadStats.abandonedTree++;

                        readProps.healingAbandoned = true;
                        readProps.abandonReason = AbandonReason.treeSize;
                        readProps.abandonedAtM = m > 0 ? m - 1 : 0;

                        List<kMerProperties> mpToReturn = new List<kMerProperties>(allocatedMerProperties);
                        foreach (kMerProperties mp in mpToReturn)
                            ReturnMerProperty(mp);
                        allocatedMerProperties.Clear();
                        List<Depth[]> depthsToReturn = new List<Depth[]>(allocatedDepths);
                        foreach (Depth[] d in depthsToReturn)
                            ReturnDepths(d);
                        allocatedDepths.Clear();
                        List<MerVariant[]> mvToReturn = new List<MerVariant[]>(allocatedMerVariants);
                        foreach (MerVariant[] mv in mvToReturn)
                            ReturnMerVariants(mv);
                        allocatedMerVariants.Clear();
                        List<List<kMerProperties>> vsToReturn = new List<List<kMerProperties>>(allocatedVariantSets);
                        foreach (List<kMerProperties> vs in vsToReturn)
                            ReturnVariantSet(vs);
                        allocatedVariantSets.Clear();
                        break;
                    }

                    if (replacementMer != null)
                    {
                        merWasChanged = replacementMer.fixType != FixTypes.fixNone;
                        mersToFirstChoice = replacementMer.mersToFirstChoice;
                    }

                    if (merWasChanged)
                    {
                        // ensure that this correction wasn't one too many for this read
                        int adjustedCost = 1;
                        int maxCost = /*m >= tailOfRead ? rewriteRegion : */rewriteRegion / 2;
                        int startingMerDepth = merDepth;
                        fixType = replacementMer.fixType;

                        if (fixType == FixTypes.fixIns)
                            adjustedCost = -replacementMer.lengthDelta;
                        if (errorModel == ErrorModel.mostlySubs && (fixType == FixTypes.fixDel || fixType == FixTypes.fixIns))
                            adjustedCost++;
                        readProps.merChangeCost[m] = adjustedCost;

                        costOfChange = CheckForRewriting(m, readProps.merChangeCost, m >= tailOfRead);

                        if (costOfChange > maxCost && replacementMer.mersToNextFix <= goodMerRunLength)
                        {
                            threadStats.abandonedRewriting++;
                            if (tracing > Tracing.traceOff)
                                TraceClumpAdd("read abandoned - rewriting: " + "cost=" + costOfChange + " @" + m);

                            readProps.healingAbandoned = true;
                            readProps.abandonReason = AbandonReason.rewriting;
                            readProps.abandonedAtM = RetreatToGoodMers(m, readProps.merChangeCost);

                            ReturnMerProperty(replacementMer);
                            break;                                                  // and give up on it
                        }
                        
                        // haven't failed the 'are we rewriting' test so go ahead and accept the change
                        currentMerCount = UpdateRead(fixType, read, m, replacementMer.variant, replacementMer.lengthDelta, currentMerCount, initialMerCount);
                        if (quals.Length > 0)
                            UpdateQuals(fixType, quals, m, replacementMer.lengthDelta);
                        merDepth = replacementMer.depth;
                        kMer = replacementMer.variant;
                        pairDepth = replacementMer.pairDepth;
                        readProps.zeroStrand[m] = false;
                        UpdateMersDepths(fixType, read, m, readProps, kMer, merDepth, !replacementMer.unbalanced);
                        if (fixType != FixTypes.fixIns)
                            changedMers++;
                        else
                            changedMers += -(replacementMer.lengthDelta);
                        // if we used more repairs than we were allowed, must have increased them for a second try at getting a good look downstream
                        if (replacementMer.fixes > followerRepairs)
                            followerRepairs = replacementMer.fixes;

                        // revise mindepths if correction has resulted in a deeper kMer
                        if (merDepth > OKDepth * 2 && ((startingMerDepth < minDepth) || (unbalancedMer && MerIsHomopolymer(startingMer))))
                        {
                            //int saveMerCount = readProps.merCount;      // in case Dels have resulted in extra bases
                            //int savePairCount = readProps.pairsCount;   // that should not normally be used
                            readProps.minDepth = 0;
                            //SetPairDepths(read, readProps);
                            //readProps.merCount = saveMerCount;
                            //readProps.pairsCount = savePairCount;
                            SetDepthThreshholds(readProps, m+1);
                            SetPairDepthThreshholds(readProps, readProps.merCount);
                            OKDepth = readProps.OKDepth;
                            minDepth = readProps.minDepth;
                            minPairDepth = readProps.minPairDepth;
                            OKPairDepth = readProps.OKPairDepth;
                            if (OKDepth > readProps.initialOKDepth*2)
                                readProps.depthsRecalculated = true;
                            if (tracing > Tracing.traceChoices)
                                TraceClumpAdd("recalculating depths and thresholds @" + m + " d=" + minDepth + "/" + OKDepth +  
                                                 " pd=" + minPairDepth + "/" + OKPairDepth + " ubr=" + readProps.unbalancedRead);
                        }

                        threadStats.replacedMers++;
                        threadStats.fixesByType[(int)fixType]++;
                    } // kMer was changed
                    
                    if (!merWasChanged && checkReason == kMerState.merBad)
                    {
                        // and check that we haven't hit a 'bad' kMer that we can't correct
                        readProps.remainingBadMers++;

                        // give up if we have an uncorrected bad kMer (only if we're not at the start of the read - leave for reverse pass)
                        if (checkReason == kMerState.merBad && firstGoodMer != -1)
                        {
                            threadStats.abandonedFutile++;
                            readProps.healingAbandoned = true;
                            readProps.abandonReason = AbandonReason.noNextMer;
                            readProps.abandonedAtM = RetreatToGoodMers(m, readProps.merChangeCost);

                            if (tracing > Tracing.traceOff)
                                TraceClumpAdd("uncorrectable errors - read abandoned @ " + readProps.abandonedAtM + "/" + m + " " + kMers.ExpandMer(kMer, merSize));

                            break;
                        }
                    }

                    if (tracing > Tracing.traceOff && merWasChanged)
                    {
                        string traceLine = "";
                        traceLine += fixTypes[(int)fixType];
                        for (int p = 0; p < m; p++)
                            traceLine += " ";
                        string originalMer;
                        string healedMer;
                        originalMer = kMers.ExpandMer(startingMer, merSize);
                        healedMer = kMers.ExpandMer(replacementMer.variant, merSize);
                        traceLine += healedMer + " (fo=" + replacementMer.goodFollowers + "/" + replacementMer.allFollowers + "/" +
                                        replacementMer.maxFollowers + "/" + replacementMer.merCount +
                                        ", fx=" + replacementMer.fixes + ", s=" + replacementMer.sum + ", d=" + replacementMer.lengthDelta + ") for " +
                                        originalMer + " @" + m + " (" + checkReason + ") " +
                                        " thm=" + thmCalls + " cc=" + costOfChange + " ";
                        int startScan = m - 20;
                        if (startScan < 0)
                            startScan = 0;
                        for (int i = startScan; i <= m; i++)
                            traceLine += readProps.merChangeCost[i];
                        TraceClumpAdd(traceLine);
                        if (replacementMer != null && replacementMer.breadcrumbs != null)
                            TraceClumpAdd(replacementMer.breadcrumbs);
                    }

                    readWasChanged = readWasChanged | merWasChanged;
                    if (replacementMer != null)
                        ReturnMerProperty(replacementMer);

                    //if (tracing > Tracing.traceOff)
                    //    TraceClumpAdd("     fmp=" + freeMerProperties.Count + " " + "fs=" + freeSequences.Count +
                    //    " mpAlloc=" + mpAllocated + " mpRet=" + mpReturned + " seqAlloc=" + seqAllocated + " seqRet=" + seqReturned);

                } // mer was healing candidate

                // mer was either OK or we tried to fix it (and either succeeded or failed)
                previousMerDepth = merDepth;
                previousPairDepth = pairDepth;
                if (firstGoodMer < 0 && merDepth >= OKDepth && (pairDepth == -1 || pairDepth >= OKPairDepth))  
                    firstGoodMer = m;

                // should we improve the quals score of a passable but unchanged kMer?
                if (quals.Length > 0 && (fixType == FixTypes.fixN || (fixType == FixTypes.fixNone && merDepth >= OKDepth && quals.Bases[m + merSize - 1] < replacementQual)))
                {
                    UpdateQuals(FixTypes.fixNone, quals, m, 0);
                }

                // move 'm' on to avoid re-scanning up to the next known choice point
                if (mersToFirstChoice > 0)
                    nextMerToCheck = m + mersToFirstChoice;
                else
                    nextMerToCheck = m + 1;

            } // for each mer in the read

            //int unreturned = 0;
            //foreach (kMerProperties mp in allocatedMerProperties)
            //    unreturned++;
            //if (unreturned > 0)
            //    Debugger.Break();
            allocatedMerProperties.Clear();

            // trim any unchecked/uncorrected bases from the end of the read (bases could have been added)
            int correctReadLength = currentMerCount + merSize - 1;
            if (read.Length > correctReadLength)
            {
                read.Length = correctReadLength;
                if (quals.Length > 0)
                    quals.Length = correctReadLength;
            }

            if (tracing > Tracing.traceOff)
            {
                TraceRead("<", read);
                if (referenceRead != null)
                    TraceClumpAdd(referenceRead);
                //string traceLine = "";
                //for (int i = 0; i < merCount; i++)
                //    traceLine += (healedTypes[i] + "\t");
                //TraceClumpAdd(traceLine);
            }

            //if (tracing > Tracing.traceOff)
            //{
            //    TraceClumpAdd("fmp=" + freeMerProperties.Count + " " + "fs=" + freeSequences.Count +
            //                    " mpAlloc=" + mpAllocated + " mpRet=" + mpReturned + " seqAlloc=" + seqAllocated + " seqRet=" + seqReturned);
            //}
            readProps.firstGoodMer = firstGoodMer;
            readProps.merCount = currentMerCount;
            readProps.changedMers = changedMers;        // won't be quite right if read abandoned but better to remember some changes were made than to say that there were none 

            if (perfTrace)
            {
                if (healReadTimer.ElapsedMilliseconds >= 50)
                {
                    string healingResult = readAfterHealing == ReadState.readIsBroken ? "still broken" : "OK";
                    TracePerf(read.ToString(), readHeader.ToString(), healingResult, healReadTimer, thmCalls);
                }
            }

            return readWasChanged;
        }

        private static Depth[] AllocateDepths(int noMerVariants)
        {
            // freeMerProperties is thread-local so no locking is needed
            Depth[] depths = null;

            // try to find a suitable length array
            int di = 0;
            for (int i = 0; i < freeDepths.Count; i++)
                if (freeDepths[i].Length >= noMerVariants)
                {
                    di = i;
                    break;
                }

            if (freeDepths.Count == 0)
            {
                depths = new Depth[noMerVariants];
            }
            else
            {
                depths = freeDepths[di];
                freeDepths.RemoveAt(di);
            }

            if (depths.Length < noMerVariants)
                Array.Resize<Depth>(ref depths, noMerVariants);

            allocatedDepths.Add(depths);
            //Interlocked.Increment(ref depthsAllocated);
            return depths;
        }

        private static void ReturnDepths(Depth[] depths)
        {
            freeDepths.Add(depths);
            allocatedDepths.Remove(depths);
            //Interlocked.Decrement(ref depthsAllocated);
        }

        private static MerVariant[] AllocateMerVariants(int noMerVariants)
        {
            // freeMerVariants is thread-local so no locking is needed
            MerVariant[] merVariants = null;

            // try to find a suitable length array
            int mvi = 0;
            for (int i = 0; i < freeMerVariants.Count; i ++)
                if (freeMerVariants[i].Length >= noMerVariants)
                {
                    mvi = i;
                    break;
                }         

            if (freeMerVariants.Count == 0)
            {
                merVariants = new MerVariant[noMerVariants];
            }
            else
            {
                merVariants = freeMerVariants[mvi];
                freeMerVariants.RemoveAt(mvi);
            }

            if (merVariants.Length < noMerVariants)
                Array.Resize<MerVariant>(ref merVariants, noMerVariants);

            allocatedMerVariants.Add(merVariants);
            //Interlocked.Increment(ref mvAllocated);
            return merVariants;
        }

        private static void ReturnMerVariants(MerVariant[] merVariants)
        {
            freeMerVariants.Add(merVariants);
            allocatedMerVariants.Remove(merVariants);
            //Interlocked.Decrement(ref mvAllocated);
        }

        private static List<kMerProperties> AllocateVariantSet()
        {
            // freeVariantSets is thread-local so no locking is needed
            List<kMerProperties> variantSet = null;

            if (freeVariantSets.Count == 0)
            {
                variantSet = new List<kMerProperties>(10);
            }
            else
            {
                variantSet = freeVariantSets.Dequeue();
            }

            allocatedVariantSets.Add(variantSet);
            //Interlocked.Increment(ref vsAllocated);
            return variantSet;
        }

        private static void ReturnVariantSet(List<kMerProperties> variantSet)
        {
            variantSet.Clear();
            freeVariantSets.Enqueue(variantSet);
            allocatedVariantSets.Remove(variantSet);
            //Interlocked.Decrement(ref vsAllocated);
        }

        private static Sequence AllocateSequence()
        {
            // freeReadContexts is thread-local so no locking is needed
            // return after thm exception is handled through kMerProperties
            Sequence thisSequence = null;

            if (freeSequences.Count == 0)
            {
                thisSequence = new Sequence(defaultReadLength);
            }
            else
            {
                thisSequence = freeSequences.Dequeue();
            }

            //Interlocked.Increment(ref seqAllocated);
            return thisSequence;
        }

        private static void ReturnSequence(Sequence thisSequence)
        {
            //if (thisSequence == null)
            //    Debugger.Break();
            freeSequences.Enqueue(thisSequence);
            //Interlocked.Decrement(ref seqAllocated);
        }

        private static kMerProperties AllocateMerProperty()
        {
            // freeMerProperties is thread-local so no locking is needed
            kMerProperties merProperty = null;

            if (freeMerProperties.Count == 0)
            {
                merProperty = new kMerProperties();
            }
            else
            {
                merProperty = freeMerProperties.Dequeue();
            }
            allocatedMerProperties.Add(merProperty);
            //Interlocked.Increment(ref mpAllocated);
            return merProperty;
        }

        private static void ReturnMerProperty(kMerProperties merProperty)
        {
            if (merProperty.fixType != FixTypes.fixNone)
            {
                ReturnSequence(merProperty.fixContext);
                merProperty.fixContext = null;
            }

            // don't need to remove from allocated list - saves on linear search in main path
            allocatedMerProperties.Remove(merProperty);
            merProperty.validVariant = true;
            merProperty.markedVariant = false;
            merProperty.fixType = FixTypes.fixNone;
            merProperty.nextFixType = FixTypes.fixNone;
            merProperty.mersToNextFix = 0;
            merProperty.mersToFirstChoice = 0;
            merProperty.perfectFix = false;
            merProperty.variant = 0;
            merProperty.depth = 0;
            merProperty.unbalanced = false;
            merProperty.pairDepth = -1;
            merProperty.sum = 0;
            merProperty.goodFollowers = 0;
            merProperty.allFollowers = 0;
            merProperty.maxFollowers = 0;
            merProperty.merCount = 0;
            merProperty.fixes = 0;
            merProperty.breadcrumbs = null;

            freeMerProperties.Enqueue(merProperty);
            //Interlocked.Decrement(ref mpAllocated);
        }

        private static void TracePerf(string unhealedRead, string faHeader, string result, Stopwatch healReadTimer, int thmCalls)
        {
            healReadTimer.Stop();
            lock (perf)
            {
                perf.WriteLine(progressReads + "\t" + faHeader + "\t" + healReadTimer.ElapsedMilliseconds + "\t" +
                               result + "\t" + unhealedRead + "\t" + thmCalls);
            }
        }

        // Takes a mer string containing one or more Ns and returns the packed form of the
        // variant that has the highest score. Returns false if no such variant could be found.
        private static bool FindBestReplacementForNs(Sequence read, int m, int OKDepth, out ulong bestPackedMer)
        {
            int bestDepth = -1;                     // best depth of all the alternatives

            // don't try replacing the Ns if there are too many of them in the mer - gets very combinatoric
            int nCount = 0;
            for (int i = 0; i < merSize; i++)
                if (read.Bases[m + i] == 'N')
                    nCount++;
            if (nCount > maxNs)
            {
                bestPackedMer = 0;
                bestDepth = 0;
                return false;
            }

            string mer = read.ToString(m, merSize);
            List<string> merVariants = new List<string>();
            bool stillMoreNs = true;
            merVariants.Add(mer);
            Sequence nextMerB = new Sequence(merSize);

            // generate mers with all replacements for the Ns
            while (stillMoreNs)
            {
                string nextMer = merVariants[0];
                nextMerB.CopyFrom(nextMer);
                int nextN = nextMer.IndexOf('N');
                if (nextN < 0)
                    break;
                foreach (char b in "ACGT")
                {
                    nextMerB.Bases[nextN] = b;
                    merVariants.Add(nextMerB.ToString());
                }
                merVariants.RemoveAt(0);

            }

            nextMerB = null;
            bestDepth = -1;
            bestPackedMer = 0;
            int depth;
            ulong variantMer;

            foreach (string mv in merVariants)
            {
                kMers.CondenseMer(mv, merSize, out variantMer);

                depth = kMersTable.GetDepthSum(variantMer);
#if STATS
                Interlocked.Increment(ref gds1Calls);
#endif

                if (depth > bestDepth)
                {
                    bestDepth = depth;
                    bestPackedMer = variantMer;
                }
            }

            return bestDepth > 0;
        }

        // Generate a set of mers from a read and calculate their read depths. 
        private static void SetMerDepths(Sequence read, ReadProperties readProps)
        {
#if STATS
            Interlocked.Increment(ref smdCalls);
#endif
            //smdTimer.Start();

            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            int pairsInRead = 0;
            if (pairsTable != null)
                pairsInRead = read.Length - pairsTable.pairFullLength + 1;
            if (pairsInRead < 0)
                pairsInRead = 0;
            readProps.pairsCount = pairsInRead;

            int balancedCount = 0;

            int firstHDUBIdx = -1;              // where did the first hdub kMer appear
            int hdubCount = 0;                  // how many did we see

            if (mersInRead > readProps.mers.Length)
            {
                Array.Resize<ulong>(ref readProps.mers, mersInRead + 100);
                Array.Resize<int>(ref readProps.depths, mersInRead + 100);
                Array.Resize<int>(ref readProps.pairDepths, mersInRead + 100);      // need extra space for debug use
                Array.Resize<int>(ref readProps.merChangeCost, mersInRead + 100);
                Array.Resize<bool>(ref readProps.zeroStrand, mersInRead + 100);
                Array.Resize<bool>(ref readProps.balancedDepths, mersInRead + 100);
            }

            readProps.merCount = mersInRead;
            readProps.deepUnbalancedPresent = false;

            bool merIsValid = false;
            ulong lastMer = 0;

            // read too short to tile for mers
            if (readLength < merSize)
            {
                readProps.merCount = 0;
                readProps.pairsCount = 0;
                return;
            }

            for (int i = 0; i < mersInRead; i++)
            {
                int plusCount = 0;
                int rcCount = 0;

                readProps.zeroStrand[i] = false;
                readProps.balancedDepths[i] = false;

                if (merIsValid)
                    merIsValid = Sequence.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = Sequence.CondenseMer(read, i, merSize, out lastMer);

                if (merIsValid)
                {
                    //gmdTimer.Start();
                    kMersTable.GetDepthCounts(lastMer, out plusCount, out rcCount);
                    //gmdTimer.Stop();
#if STATS
                    Interlocked.Increment(ref gdcCalls);
#endif
                    int summedDepth = plusCount + rcCount;
                    readProps.depths[i] = summedDepth;
                    readProps.mers[i] = lastMer;

                    if (amplicons)
                    {
                        int maxCount = Math.Max(plusCount, rcCount);
                        plusCount = maxCount;
                        rcCount = maxCount;
                    }

                    if (plusCount <= 1 || rcCount <= 1)
                    {
                        readProps.zeroStrand[i] = true;
                        if (!readProps.hmZeroPresent)
                            readProps.hmZeroPresent = MerIsHomopolymer(lastMer);
                    }

                    if (summedDepth >= requestedMinDepth && (plusCount < rcCount * balancedFactor && rcCount < plusCount * balancedFactor))
                    {
                        readProps.balancedDepths[i] = true;
                        balancedCount++;
                    }

                    if (kMersTable.hdubFilter.Contains(lastMer))
                    {
                        //if (firstHDUBIdx < 0)
                        //    firstHDUBIdx = i;
                        hdubCount++;
                        //lock (hdubLocs)
                        //{
                        //    hdubLocs[i]++;
                        //}
                    }
                }
                else
                {
                    readProps.depths[i] = 0;
                    readProps.mers[i] = 0;
                }
            }

            int hdubThreshhold = Math.Min(mersInRead - firstHDUBIdx, 5);
            readProps.deepUnbalancedPresent = hdubCount >= hdubThreshhold;
            readProps.balancedMerCount = balancedCount;

            //smdTimer.Stop();
        }

        private static void ResetSetMerDepths(Sequence read, int endingM, ReadProperties readProps)
        {
            bool merIsValid = false;
            ulong lastMer = 0;
            int[] depths = readProps.depths;
            ulong[] mers = readProps.mers;
            bool[] balanced = readProps.balancedDepths;
            int balancedCount = 0;

            for (int i = 0; i < endingM; i++)
            {
                int plusCount = 0;
                int rcCount = 0;

                if (merIsValid)
                    merIsValid = Sequence.CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                else
                    merIsValid = Sequence.CondenseMer(read, i, merSize, out lastMer);

                if (merIsValid)
                {
                    kMersTable.GetDepthCounts(lastMer, out plusCount, out rcCount);
#if STATS
                    Interlocked.Increment(ref gdcCalls);
#endif
                    int summedDepth = plusCount + rcCount;
                    //if (tracingThisRead)
                    //    TraceClumpAdd("resetting " + kMers.ExpandMer(mers[i], merSize) + " with " + kMers.ExpandMer(lastMer, merSize) + " @" + i);

                    depths[i] = summedDepth;
                    mers[i] = lastMer;
                    if (summedDepth >= requestedMinDepth && (plusCount < rcCount * 4 && rcCount < plusCount * 4))
                    {
                        balanced[i] = true;
                        balancedCount++;
                    }
                }
                else
                {
                    depths[i] = 0;
                    mers[i] = 0;
                }
            }

            readProps.balancedMerCount = balancedCount;

        }

        private static void SetPairDepths(Sequence read, ReadProperties readProps)
        {
            int pairsInRead = readProps.pairsCount;

            if (pairsInRead == 0)
                return;

#if STATS
            Interlocked.Increment(ref spdCalls);
#endif
            //spdTimer.Start();

            ulong pair = 0;
            bool pairValid = false;
            int[] pairDepths = readProps.pairDepths;

            for (int i = 0; i < pairsInRead; i++)
            {
                if (pairValid)
                    pairValid = kMerPairs.ConstructPairIncremental(read, i, merPairGap, pair, out pair);
                else
                    pairValid = kMerPairs.ConstructPair(read, i, merPairGap, out pair);

                if (pairValid)
                {
                    //gpdTimer.Start();
                    pairDepths[i] = pairsTable.GetPairDepth(pair);
                    //gpdTimer.Stop();
                }
                else
                    pairDepths[i] = 0;
            }

            //spdTimer.Stop();
        }

        //private static void CullPairDepths(int pairsInRead, int[] pairDepths, int[] depths, int minDepthWanted, int maxDepthWanted)
        //{
        //    // index if starting base of the kMer at the end of the pair
        //    int secondMerOffset = pairsTable.pairFullLength - merSize;

        //    for (int i = 0; i < pairsInRead; i++)
        //    {
        //        if (pairDepths[i] == 0)
        //            continue;

        //        int startingMerDepth = depths[i];
        //        int endingMerDepth = depths[i + secondMerOffset];

        //        // don't count pairs that are based on kMers whose depth is too high or too low
        //        if (startingMerDepth < minDepthWanted || endingMerDepth < minDepthWanted || startingMerDepth > maxDepthWanted || endingMerDepth > maxDepthWanted)
        //            pairDepths[i] = 0;
        //    }
        //}

        // Determines whether the read needs healing and sets the initial 'good' (OKDepth) and 'passable' (minDepth) depth levels.
        private static void SetDepthThreshholds(ReadProperties readProps, int lengthToScan)
        {
            int noOKMers = 0;                       // no. of mers that are better than OK depth
            int minDepth = 0;                       // minimum allowable depth (for 'any' follower
            int OKMerMean = 0;                      // average of counts for these merss)
            int OKDepth = 0;                        // depth for 'good' followers and detecting places to try error correction
            int lowestBalancedDepth = int.MaxValue; // lowest balanced kMer depth

            bool readUnbalanced = readProps.balancedMerCount < 1;
            int[] merDepths = readProps.depths;

            // calculate the harmonic mean of the depth of the passable mers (for depth capping)
            double invDepthSum = 0.0;
            int countedMers = 0;

            for (int m = 0; m < lengthToScan; m++)
            {
                int depth = merDepths[m];
                if (readUnbalanced)
                {
                    if (depth > requestedMinDepth)
                    {
                        invDepthSum += 1.0f / (double)depth;
                        countedMers++;
                    }
                }
                else
                {
                    if (depth > 0 && readProps.balancedDepths[m] && !readProps.zeroStrand[m])
                    {
                        invDepthSum += 1.0f / (double)depth;
                        countedMers++;

                        if (depth < lowestBalancedDepth)
                            lowestBalancedDepth = depth;
                    }
                }
            }

            int averageDepth = averageLoadedMerDepth; 
            if (countedMers > 0)
                averageDepth = (int)((double)countedMers / invDepthSum);

            // set noise floor at 5% of average (only important for unbalanced reads)
            int lowestWantedDepth = averageDepth / 6;
            // or the lowest balanced
            if (!readUnbalanced)
                lowestWantedDepth = lowestBalancedDepth / 6;

            // calculate the harmonic mean of the depths below the average to get a better estimation of the lower part of the depth range
            // this is intended to reduce the effect of repeat regions and droop. 
            invDepthSum = 0.0;
            for (int m = 0; m < lengthToScan; m++)
            {
                int depth = merDepths[m];
                // cap depths to trim peaks
                if (depth > averageDepth)
                    depth = averageDepth;

                if (depth >= lowestWantedDepth)
                {
                    if (readUnbalanced || readProps.balancedDepths[m])
                    {
                        noOKMers++;
                        invDepthSum += 1.0f / (double)depth;
                    }
                }
            }

            if (noOKMers > 0)
                OKMerMean = (int)((double)noOKMers / invDepthSum);
            else
                OKMerMean = averageDepth;
            if (OKMerMean == 0)
                OKMerMean = averageDepth;

            OKDepth = OKMerMean / 3;                                    // OK depth to the average of the lower depths
            minDepth = OKDepth / 2;                                     // and 'min' to 1/2 of this OK level *** check ***
            if (minDepth > lowestBalancedDepth)
                minDepth = lowestBalancedDepth - 1;                     // never set acceptable/bad boundary higher than the lowest balanced kMer

            readProps.minDepth = minDepth;
            readProps.OKDepth = OKDepth;
            if (readProps.initialOKDepth < 0)
                readProps.initialOKDepth = OKDepth;

            readProps.unbalancedRead = readUnbalanced;
        }

        private static void SetPairDepthThreshholds(ReadProperties readProps, int lengthToScan)
        {
            int OKPairMean = 0;                     // and the same for the pairs
            int noOKPairs = 0;                      // no. of pairs that are better than OK depth
            int minPairDepth = 0;                   // minimum depth allowed for a pair
            int OKPairDepth = 0;                    // 'good' pair depth

            int pairsInRead = readProps.pairsCount;
            int[] merDepths = readProps.depths;
            int[] pairDepths = readProps.pairDepths;

            if (pairsInRead > 0)
            {
                int lowestWantedDepth = readProps.minDepth / 2;
                int maxWantedDepth = readProps.OKDepth * 100;
                //// remove any erroneous/dubious pair depths now that minDepth has been calculated
                //CullPairDepths(pairsInRead, pairDepths, merDepths, lowestWantedDepth, maxWantedDepth);

                // calculate the mean for the passable pairs
                long sumPairDepths = 0;
                long nonZeroPairCount = 0;
                for (int m = 0; m < pairsInRead; m++)
                {
                    int pairDepth = pairDepths[m];
                    if (pairDepth > 0)
                    {
                        sumPairDepths += pairDepth;
                        nonZeroPairCount++;
                    }
                }
                int averagePairDepth = averageLoadedPairDepth; 
                if (nonZeroPairCount > 0)
                    averagePairDepth = (int)(sumPairDepths / nonZeroPairCount);
                averagePairDepth = Math.Max(averagePairDepth, readProps.OKPairDepth * 2);
                int lowestWantedPairDepth = averagePairDepth * 5 / 100;
                if (lowestWantedPairDepth < lowestWantedDepth / 2)
                    lowestWantedPairDepth = lowestWantedDepth / 2;
                if (lowestWantedPairDepth == 0)
                    lowestWantedPairDepth = 1;

                // calculate the harmonic mean of the depths below the average to get a better estimation of the lower part of the depth range
                // this is intended to reduce the effect of repeat regions and droop
                double invSumPairs = 0.0f;
                for (int m = 0; m < pairsInRead; m++)
                {
                    int pairDepth = pairDepths[m];
                    if (pairDepth > averagePairDepth)
                        pairDepth = averagePairDepth;
                    if (pairDepth >= lowestWantedPairDepth)
                    {
                        noOKPairs++;
                        invSumPairs += 1.0f / (double)pairDepth;
                    }
                }

                if (noOKPairs > 0)
                    OKPairMean = (int)((double)noOKPairs / invSumPairs);            // harmonic mean of all 'non-zero' pairs
                else
                    OKPairMean = readProps.OKDepth;                                 // no pairs deep enough so go with the kMer depths numbers
            }

            OKPairDepth = OKPairMean / 2;                                           
            OKPairDepth = Math.Max(OKPairDepth, readProps.OKPairDepth);             // could be updating pairs after change
            minPairDepth = OKPairDepth / 2;
            if (minPairDepth > readProps.minDepth)
                minPairDepth = readProps.minDepth;

            if (pairsInRead > 0 && minPairDepth < 1)
                minPairDepth = minLoadDepth;
            minPairDepth = Math.Max(minPairDepth, readProps.minPairDepth);

            //OKPairDepth = OKDepth * 2 / 3;
            //minPairDepth = OKPairDepth / 2;
            //if (minPairDepth == 0)
            //    minPairDepth = 1;

            readProps.minPairDepth = minPairDepth;
            readProps.OKPairDepth = OKPairDepth;
        }

        private static ReadState CheckReadForProblems(Sequence read, ReadProperties readProps)
        {
            // now that the min and OK depths for this read have been established, scan through the depths once again looking for problems
            int minDepth = readProps.minDepth;
            int minPairDepth = readProps.minPairDepth;
            int mersInRead = readProps.merCount;
            int pairsInRead = readProps.pairsCount;

            if (readProps.OKDepth > requestedMaxDepth)
                return(ReadState.readTooDeep);

            int zeroCount = 0;                      // how many below-min depths
            int zeroPairCount = 0;                  // how many below-min-depth pairs
            int zeroStrandCount = 0;                // how many 0 depths on a strand
            int deepestDepth = 0;                   // deepest depth encountered

            ReadState readNeedsHealing = ReadState.readOK;

            int prevMerDepth = 0;
            for (int m = 0; m < mersInRead; m++)
            {
                int thisDepth = readProps.depths[m];
                if (thisDepth < minDepth)                               // trigger healing if drop below minimum for read
                    zeroCount++;
                else
                {
                    if (errorModel == ErrorModel.indelsCommon && thisDepth < prevMerDepth / 2 && MerIsHomopolymerEnd(readProps.mers[m]))  // or we've hit a sudden drop after an HP run
                        readNeedsHealing = ReadState.readNeedsChecking;
                }
                if (readProps.zeroStrand[m])
                    zeroStrandCount++;
                prevMerDepth = thisDepth;
                if (thisDepth > deepestDepth)
                    deepestDepth = thisDepth;
            }

            for (int p = 0; p < pairsInRead; p++)
            {
                if (readProps.pairDepths[p] < minPairDepth)
                    zeroPairCount++;
            }

            readProps.remainingBadMers = zeroCount;

            if (zeroStrandCount > 0)
            {
                //if (readUnbalanced)
                readNeedsHealing = ReadState.readNeedsChecking;
                //else
                //    readNeedsHealing = ReadState.readIsBroken;
            }
            if (zeroCount > 0 | (pairsInRead > 0 && zeroPairCount > 0))
                readNeedsHealing = ReadState.readIsBroken;
            if (zeroCount == 0 && (pairsInRead > 0 && zeroPairCount > 0))
                readNeedsHealing = ReadState.readNeedsChecking;

            return readNeedsHealing;
        }

        private static string TraceBrokenReason(ReadProperties readProps)
        {
            string reason = "";

            if (readProps.finalReadState == ReadState.readNotLongEnough)
                reason = "too short/too many corrections.";
            if (readProps.finalReadState == ReadState.readIsBroken)
                reason += " still broken.";
            if (readProps.abandonReason == AbandonReason.noNextMer)
                reason += " no next kMer.";
            if (readProps.abandonReason == AbandonReason.rewriting)
                reason += " rewriting.";
            if (readProps.abandonReason == AbandonReason.tooManyNs)
                reason += " too many Ns.";
            if (readProps.abandonReason == AbandonReason.treeSize)
                reason += " tree size.";

            return reason;
        }

        private static bool MerIsHealingCandidate(Sequence read, ReadProperties readProps, int m, ulong packedMer, bool unbalancedMer, bool tiltedMer, char followingBase, int depth, int previousDepth, int pairDepth, int previousPairDepth, bool initialiseDepths, out kMerState checkReason)
        {
            checkReason = kMerState.merUnknown;
            bool unbalancedRead = readProps.unbalancedRead;
            // for calls from CountFollowers, we just want to keep the depth thresholds we have already
            int minDepth = readProps.minDepth;
            int OKDepth = readProps.OKDepth;
            int minPairDepth = readProps.minPairDepth;
            int OKPairDepth = readProps.OKPairDepth;

            // adjust depths if it looks we're being affected by an HP artefact
            if (readProps.hmZeroPresent && tiltedMer)
            {
                minDepth = minDepth * 2 / 3;
                if (minDepth < 1)
                    minDepth = 1;
                minPairDepth = minPairDepth * 2 / 3;
                if (minPairDepth < 1)
                    minPairDepth = 1;
                OKDepth = OKDepth * 2 / 3;
                OKPairDepth = OKPairDepth * 2 / 3;
            }           

            // be careful with unbalanced kMers in a balanced read
            bool merUnbalanced = unbalancedMer && !unbalancedRead && depth < previousDepth * 3 / 4;
            // could we be at the end of an HP run (for indel-common data)
            //bool checkHP = errorModel == ErrorModel.indelsCommon && previousDepth >= OKDepth && depth < previousDepth / 2 && MerIsHomopolymerEnd(packedMer, followingBase);
            // a sudden drop is worth checking out as well (look for those places where we're coming out of a repeat and hit an error that takes us down the wrong path)
            bool checkDepthDrop = previousDepth >= OKDepth && depth < previousDepth * 2 / 3;
            // bad if both are below minimum depth for a read (with some caveats)
            bool merDepthBad = MerIsBad(depth, minDepth, OKDepth, pairDepth, minPairDepth, merUnbalanced, previousDepth, previousPairDepth, m); 
            // unsure if either are below the OK level for the read
            bool merDepthMiddling = depth < OKDepth || (pairDepth != -1 && pairDepth < OKPairDepth);    

            //// tried using quals as error clues but this halved performance for no real improvement
            //bool qualPoor = false;
            //if (healingQuals != null && healingQuals.Length > 0)
            //    qualPoor = healingQuals.Bases[m + merSize - 1] < lowestAcceptableQual;

            bool merHasAlternativesKnown = false;
            bool merHasAlternatives = false;

            int minDepthForAlternative = minDepth;
            // coming from a very deep kMer, so be a bit picky about the minimum depth required for a viable alternative. 
            // Should be much the same as the previous depth but allow for forks
            if (previousDepth > averageLoadedMerDepth * 10)
                minDepthForAlternative = previousDepth / 5;
            if (minDepthForAlternative > OKDepth)
                minDepthForAlternative = minDepth;

            // 'bad' kMers 
            if (merDepthBad)
            {
                checkReason = RatchetMerState(checkReason, kMerState.merBad);
#if STATS
                MIHCCauses["bad"]++;
#endif
            }

            if (merUnbalanced && checkReason == kMerState.merUnknown)
            {
                if (!merHasAlternativesKnown)
                {
                    merHasAlternatives = MerHasAlternatives(packedMer, minDepthForAlternative, unbalancedMer, VariantWanted.varyLast);
                    merHasAlternativesKnown = true;
                }

                // unexpected unbalanced homopolymer kMer - likely to be an Illumina poor-quality artefact
                if (merHasAlternatives && MerIsHomopolymer(packedMer))
                {
                    checkReason = RatchetMerState(checkReason, kMerState.merBad);
#if STATS
                    MIHCCauses["bad ubhp alt"]++;
#endif
                }

                if (checkReason == kMerState.merUnknown && merHasAlternatives)
                {
                    checkReason = RatchetMerState(checkReason, kMerState.merUnsure);
#if STATS
                    MIHCCauses["unbalanced alt"]++;
#endif
                }
            }

            // depth 'unsure' - either kMer or pair is below desired depth
            if (merDepthMiddling && checkReason == kMerState.merUnknown)
            {
                // be quite forgiving if either of the depths is OK (or the depth is much the same as previously) and there is no alternative
                if ((depth >= OKDepth || pairDepth >= OKPairDepth || depth > previousDepth * 3 / 4)) 
                {
                    if (!merHasAlternativesKnown)
                    {
                        merHasAlternatives = MerHasAlternatives(packedMer, minDepthForAlternative, unbalancedMer, m == 0 ? VariantWanted.varyAnyOne : VariantWanted.varyLast);
                        merHasAlternativesKnown = true;
                    }
                    if (merHasAlternatives)
                    {
                        checkReason = RatchetMerState(checkReason, kMerState.merCheck);
#if STATS
                        MIHCCauses["check middling alt"]++;
#endif
                    }
                    else
                    {
                        checkReason = RatchetMerState(checkReason, kMerState.merOK);
#if STATS
                        MIHCCauses["OK middling no alt"]++;
#endif
                    }
                }

                if (checkReason == kMerState.merUnknown)
                {
                    checkReason = kMerState.merUnsure;
#if STATS
                    MIHCCauses["unsure middling"]++;
#endif
                }
            }

            //if (checkHP)
            //    checkReason = RatchetMerState(checkReason, kMerState.merCheck);

            if (checkDepthDrop && checkReason == kMerState.merUnknown)
            {
                if (!merHasAlternativesKnown)
                {
                    merHasAlternatives = MerHasAlternatives(packedMer, minDepthForAlternative, unbalancedMer, VariantWanted.varyLast);
                    merHasAlternativesKnown = true;
                }
                if (merHasAlternatives)
                {
                    checkReason = RatchetMerState(checkReason, kMerState.merCheck);
#if STATS
                    MIHCCauses["check depthDrop alt"]++;
#endif
                }
                else
                {
                    checkReason = RatchetMerState(checkReason, kMerState.merOK);
#if STATS
                    MIHCCauses["OK depthDrop no alt"]++;
#endif
                }
            }

            //if (qualPoor)
            //    checkReason = ChangeMerState(checkReason, kMerState.merPoor);

            // look for alternatives for seemingly good initial kMers if things seem to go pear-shaped soon afterwards
            if (m == 0 && checkReason == kMerState.merUnknown)
            {
                bool peekUnbalanced;
                bool peekTilted;
                int peekDepth = kMersTable.GetDepthSum(PeekMer(packedMer, m, read, readProps.merCount + merSize - 1, merSize-1), minDepth, out peekUnbalanced, out peekTilted);

                if (peekDepth < minDepth)
                {
                    // are there viable alternatives to this seemingly good kMer?
                    if (!merHasAlternativesKnown)
                    {
                        merHasAlternatives = MerHasAlternatives(packedMer, minDepthForAlternative, unbalancedMer, VariantWanted.varyAnyOne);
                        merHasAlternativesKnown = true;
                    }
                    if (merHasAlternatives)
                    {
                        checkReason = RatchetMerState(checkReason, kMerState.merCheck);
#if STATS
                        MIHCCauses["check m==0 alt"]++;
#endif
                    }
                }
            }

            // and peek ahead to see if we're about to get into trouble (and we've just done this for m=0)
            if (checkReason == kMerState.merUnknown && m != 0)
            {
                bool peekUnbalanced;
                bool peekTilted;
                int peekDepth = kMersTable.GetDepthSum(PeekMer(packedMer, m, read, readProps.merCount+merSize-1, 2), minDepth, out peekUnbalanced, out peekTilted);
#if STATS
                Interlocked.Increment(ref gds4Calls_Peek);
#endif

                if (peekDepth < minDepth)
                {
                    if (!merHasAlternativesKnown)
                    {
                        merHasAlternatives = MerHasAlternatives(packedMer, minDepthForAlternative, unbalancedMer, m == 0 ? VariantWanted.varyAnyOne : VariantWanted.varyLast);
                        merHasAlternativesKnown = true;
                    }
                    if (merHasAlternatives)
                    {
                        checkReason = RatchetMerState(checkReason, kMerState.merCheck);
#if STATS
                        MIHCCauses["check peek alt"]++;
#endif
                    }
                }
            }

            if (checkReason == kMerState.merUnknown)
                checkReason = kMerState.merOK;

            return checkReason > kMerState.merOK;
        }

        private static kMerState RatchetMerState(kMerState currentState, kMerState newState)
        {
            if (newState > currentState)
                currentState = newState;
            return currentState;
        }

        private static ulong PeekMer(ulong mer, int m, Sequence read, int readLength, int peekBases)
        {
            ulong shiftedMer = mer;
            int basesLeft = readLength - m - merSize;
            if (peekBases > basesLeft)
                peekBases = basesLeft;

            for (int i = 0; i < peekBases; i++)
            {
                long nextBase = kMers.BaseCharToInt(read.Bases[m + merSize + i]);
                shiftedMer = (shiftedMer << 2) | (ulong)(nextBase << (64 - merSize * 2));
            }

            return shiftedMer;
        }

        // only start of read is bad
        private static bool StartOfReadIsBad(ReadProperties readProps)
        {
            int minDepth = readProps.minDepth;
            int firstAcceptableDepth = -1;
            int lastBadDepth = -1;

            for (int i = 0; i < readProps.merCount; i++)
            {
                int depth = readProps.depths[i];
                if (depth < minDepth)
                    lastBadDepth = i;
                else
                {
                    if (firstAcceptableDepth < 0)
                        firstAcceptableDepth = i;
                }
            }

            return lastBadDepth < firstAcceptableDepth;
        }

        private static bool MerHasAlternatives(ulong mer, int minDepth, bool unbalancedStartingMer, VariantWanted variantsWanted)
        {
            int acceptable = 0;
            int start = 0;
            if (variantsWanted == VariantWanted.varyLast)
                start = merSize - 1;

            ulong baseMask = 0xc000000000000000;

            for (int m = start; m < merSize; m++)
            {
                ulong merWithHole = mer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    int targetMinDepth = minDepth;

                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong alternateMer = merWithHole | newBase;
                    if (alternateMer == mer)
                        continue;

                    bool unbalancedAlternative;
                    bool tiltedAlternative;
                    int depth = kMersTable.GetDepthSum(alternateMer, minDepth, out unbalancedAlternative, out tiltedAlternative);
#if STATS
                    Interlocked.Increment(ref gds4Calls_Alts);
#endif

                    if (!unbalancedStartingMer && unbalancedAlternative)
                        continue;

                    if (unbalancedStartingMer && !unbalancedAlternative)
                        targetMinDepth = minDepth / 2;

                    if (depth >= targetMinDepth)
                    {
                        acceptable++;
                        break;
                    }
                }
            }

            return acceptable > 0;
        }

        // The only pair direction supported is backwards, as downstream errors could be blocking a pair while the current kMer is good.
        // The kMers involved in the pairs must both be at least minDepth deep. If it was possible to generate a backwards pair,
        // its depth is returned. If we're being called from TryHealingRead, the kMers and depths will be valid and will be used.
        // If we're being called with a modified read, we'll calculate the kMers and depths afresh.
        private static int GetPairDepth(Sequence read, ulong kMer, ulong[] mers, int[] depths, int m, int minDepth)
        {
            int pairDepth = -1;

            if (pairsTable == null)
                return -1;

#if STATS
            Interlocked.Increment(ref gpdCalls);
#endif

            ulong merPair;                      // a mer pair ending at 'm' (or possibly starting at 'm' if we're close to the start of the read)
            int firstMerDepth = 0;              // depth of kMer at start of the pair

            // can we have a pair finishing with the end of this mer (only backwards allowed now)
            if ((m + merSize) < merPairLength)
            {
#if STATS
                Interlocked.Increment(ref gpdNoRoom);
#endif
                return -1;
            }

            // if we're given a read, we'll generate the pair and its related kMers afresh
            if (read != null)
            {
                ulong firstMer;                     // a full k-mer corresponding to the first fragment in the pair string
                ulong secondMer;                    // and the kMer at the end of the pair string

                // generate the k-mer pair (separated by a 'gap' gap)
                if (!GenerateMerPairFromRead(read, m, PairDirection.backwards, out merPair, out firstMer, out secondMer))
                {
                    // found non-ACGT bases when constructing the pair
                    pairDepth = -1;
                    return pairDepth;
                }

                firstMerDepth = kMersTable.GetDepthSum(firstMer);
#if STATS
                Interlocked.Increment(ref gpdRead);
                Interlocked.Increment(ref gds1Calls);
#endif
            }
            else
            // already have kMers and depths available so take the fast path
            { 
                int firstMerIdx = m + merSize - merPairLength;
                ulong firstMer = mers[firstMerIdx];
                firstMerDepth = depths[firstMerIdx];
                merPair = (firstMer & kMerPairs.firstFragmentMask) | ((kMer >> (64 - merSize * 2)) & kMerPairs.lastFragmentMask);
#if STATS
                Interlocked.Increment(ref gpdFast);
#endif
            }

            if (firstMerDepth < minDepth)
            {
                pairDepth = -1;
                return pairDepth;
            }
            
            // look up the pair (in canonical form) in the pairs table, and return not found if it isn't there or if there aren't enough repetitions
            pairDepth = pairsTable.GetPairDepth(merPair);

            return pairDepth;
        }

        private static bool GenerateMerPairFromRead(Sequence read, int m, PairDirection backwardsOrForwards, out ulong merPair, out ulong startMer, out ulong endMer)
        {
            // Generate an individual canonical packed k-mer pair from a (possibly tentative) read. 
            // The second 16-mer is normally the last part of the k-mer starting at 'm', the first 16-mer starts (16+gap) bases before this.
            // If 'm' is too close to the start of the read, we'll try for a pair starting with the k-mer starting at 'm' and going into the uncorrected part of the read.
            // The pair string always has to include the last base in the 'm' kMer as it's the one being varied.  
            //  
            // backwards (into more correct/checked start region)         
            // <----                                m       
            // ..........1234567890123456....................1234567890123456......... (for k=25 and gap=20) 
            //           ---- start mer ----------  -------- end mer --------    
            //
            // forwards (less reliable as it includes the more error-prone (and uncorrected) tail of the read)
            // ---->     m                                         
            // ..........---------1234567890123456....................1234567890123456...... 
            //           ---- start mer ----------           -------- end mer --------

            merPair = 0;
            startMer = 0;                               // a full k-mer based on the starting 16-mer - used to detect errors in first part of pair
            endMer = 0;                                 // and the full kMer at the end of the pair fragment
            bool mersOK;                                // no problems constructing the kMers

            int startOfFirstFragment;
            if (backwardsOrForwards == PairDirection.backwards)
                startOfFirstFragment = m + merSize - merPairLength;
            else
                startOfFirstFragment = m + merSize - merPairFragmentLength;

            mersOK = kMerPairs.ConstructPair(read, startOfFirstFragment, merSize, merPairGap, out merPair, out startMer, out endMer);

            return mersOK;
        }


        private static bool MerIsHomopolymerEnd(ulong packedMer, char followingBaseChar)
        {
            // test whether we're at the end of a homopolymer (....XXX[y] or ....XXXy) (only used for 454 reads)

            // all possible homopolymer ends in packed form
            const ulong a3 = 0x0;               // AAA
            const ulong c3 = 0x15;              // CCC
            const ulong g3 = 0x2a;              // GGG
            const ulong t3 = 0x3f;              // TTT

            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shift mer
            ulong low3Bases = shiftedMer & 0x3f;                    // last 3 bases XXX

            // does current mer end in an HP (and followed by a non-HP base)
            if (low3Bases == a3 || low3Bases == c3 || low3Bases == g3 || low3Bases == t3)
            {
                // mer ends in an HP string
                ulong lastBase = low3Bases & 0x3;
                ulong followingBase = (ulong)kMers.BaseCharToInt(followingBaseChar);
                //ulong followingBase = (ulong)("ACGT".IndexOf(followingBaseChar));
                // so return true if the next base in the read isn't an extension of the HP string (ie we're at the end of the HP)
                return (lastBase != followingBase);
            }
            else
            {
                // mer doesn't end in an HP, are we just after the end of an HP (...XXXy)
                shiftedMer = shiftedMer >> 2;       // drop off the last base (we know it isn't the same as the bases before it)
                low3Bases = shiftedMer & 0x3f;      // get the XXX bases (from a possible initial XXXy mer)
                return (low3Bases == a3 || low3Bases == c3 || low3Bases == g3 || low3Bases == t3);
            }
        }

        private static bool MerIsHomopolymerEnd(ulong packedMer)
        {
            // test whether the end of the kMer is the end of a homopolymer (...XXXy). 
            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shifted mer
            ulong lowBases = shiftedMer & 0xff;                     // bases XXXy
            ulong lastBase = lowBases & 0x3;                        // last base
            ulong possibleHPBases = lowBases >> 2;                  // XXXX
            ulong lowHPBase = possibleHPBases & 0x3;                // X

            if (lastBase == lowHPBase)                              // can't be at the end of an HP if the last base is the same as the previous one
                return false;

            const ulong a3 = 0x0;
            const ulong c3 = 0x15;
            const ulong g3 = 0x2a;
            const ulong t3 = 0x3f;

            // say it's an HP end if there are 3 bases the same followed by a different one...
            return (possibleHPBases == a3) || (possibleHPBases == c3) || (possibleHPBases == g3) || (possibleHPBases == t3);
        }

        private static bool MerIsHomopolymer(ulong packedMer)
        {
            // test whether the end of the kMer is a homopolymer (...XXX). 
            ulong shiftedMer = packedMer >> (32 - merSize) * 2;     // right shifted mer
            ulong low3Bases = shiftedMer & 0x3f;                    // bases XXX
            const ulong a3 = 0x0;
            const ulong c3 = 0x15;
            const ulong g3 = 0x2a;
            const ulong t3 = 0x3f;

            return (low3Bases == a3) || (low3Bases == c3) || (low3Bases == g3) || (low3Bases == t3);
        }

        private static bool MerIsBad(int depth, int minDepth, int OKDepth, int pairDepth, int minPairDepth, bool merUnbalanced, int previousDepth, int previousPairDepth, int m)
        {
            bool merDepthBad = depth < minDepth || (depth < OKDepth && (pairDepth != -1 ? pairDepth < minPairDepth : false));  

            // but if this depth is much the same as the previous one (which we must have decided was OK), we'll let it be called 'poor' rather than 'bad'
            if (merDepthBad && (previousDepth > 0 && depth > previousDepth * 3 / 4) && (pairDepth == -1 || pairDepth > previousPairDepth * 3 / 4))
            {
                merDepthBad = false;    
#if STATS
                Interlocked.Increment(ref depthBadButRedeemed);
#endif
                //if (tracing == Tracing.traceFollowers)
                //    TraceClumpAdd("redeeming previous: d=" + depth + "/" + minDepth + "(" + previousDepth + ") pd=" + pairDepth + "/" + minPairDepth + "(" + previousPairDepth + ") @" + m);
            }

            return merDepthBad;
        }

        // Does this mer look like a low-rep error variant of a high-rep mer?
        private static bool MerIsLowRepVariant(ulong packedMer, int merDepth)
        {
            int merFill = 64 - merSize * 2;
            ulong rshiftedMer = packedMer >> merFill & 0xFFFFFFFFFFFFFFFC;
            int sumVariantCounts = 0;
            // generate all variants of this mer, look them up in the table and sum their counts
            for (ulong b = 0; b <= 3; b++)
            {
                ulong variantMer = (rshiftedMer | b) << merFill;
                int variantCountPlus = 0;
                int variantCountRC = 0;
                kMersTable.GetDepthCounts(variantMer, out variantCountPlus, out variantCountRC);
                sumVariantCounts += variantCountPlus + variantCountRC;
            }
            // and reject the current mer if it's counts are tiny relative to the sum for the variants (10%)
            return (sumVariantCounts > 0) && (((merDepth * 100) / sumVariantCounts) <= 10);
        }

        // Try healing the current mer (in packedMer). Returns the preferred variant kMer - possibly an unchanged kMer if this is the best choice
        // This function can be recursively called through CountFollowers, controlled by the repairsLeft parameter 
        static bool TryHealingMer(kMerState checkReason, bool subFixesOnly, bool recursiveCall, ulong kMer, int merDepth, int pairDepth, int previousMerDepth, int previousPairDepth, bool unbalancedMer,
                                            Sequence read, ReadProperties readProps, int m, int currentMerCount, int repairsLeft, FixTypes previousFixType, int mersToFirstChoice, out kMerProperties replacementMer,
                                            bool outer, int indentDepth, ref int thmCalls)
        {
            replacementMer = null;

            if (traceClump.Count > 1000)
                WriteTraceClump(traceClump);

            Stopwatch healMerTimer = null;
            if (tracing > Tracing.traceOff)
            {
                healMerTimer = new Stopwatch();
                healMerTimer.Start();
            }

            thmCalls++;
            if (thmCalls > maxTHMAllowed)
                return false;

            List<kMerProperties> variantSet = AllocateVariantSet();

            int bestChoiceIdx = 0;                                              // return this variant as the choice
            bool bestChoiceFound = false;                                       // set whenever we've decided on the best choice
            string bestChoiceReason = null;                                     // and why we chose it (for tracing)
            bool bestChoicePerfect = false;                                     // did we make a 'perfect' choice?

            string tracedStartingMer = null;
            if (tracing > Tracing.traceOff)
            {
                bestChoiceReason = "";
                tracedStartingMer = kMers.ExpandMer(kMer, merSize);
            }
            //if (startingMer == "GATTCCCATCGTCCACCCCATATTC" && tracing == traceFollowers)
            //    Debugger.Break();
            if (tracing == Tracing.traceFollowers)
                TraceClumpAdd(Indent(indentDepth) + "THM: heal " + tracedStartingMer + " @" + m + " " + checkReason + " d=" + merDepth + "/" + readProps.minDepth + "/" + readProps.OKDepth +
                                               ", pd=" + pairDepth + "/" + readProps.minPairDepth + "/" + readProps.OKPairDepth + "; prev=" + previousMerDepth + (unbalancedMer ? " unbalanced" : " balanced"));

            bool cpvOK = CollectPlausibleVariants(checkReason, subFixesOnly, repairsLeft, variantSet, readProps, m, currentMerCount, kMer, merDepth, pairDepth, unbalancedMer, read, 
                                                        previousMerDepth, previousPairDepth, previousFixType, mersToFirstChoice, indentDepth, ref thmCalls);

            if (!cpvOK)
            {
                ReturnVariantSet(variantSet);
                return false;
            }

            //// how many fixes are viable for the first mer (stops us just rewriting the entire read when there are too many errors at the start) 
            //int allowableFixes = (currentMerCount - m) / 2;
            // used to give a preference to Dels when they are extending HP runs (454-like data only, -hp option set)
            bool HPExtendingDel = false;

            //// initial pass through the variants
            //foreach (kMerProperties variant in variantSet)
            //{
            //    // adjust fixes count for HP-extending dels if we're trying harder to fix indels
            //    if (errorModel == ErrorModel.indelsCommon)
            //    {
            //        // treat HP extending dels specially... by reducing their fix count
            //        if (variant.fixType == FixTypes.fixDel && MerIsHomopolymer(variant.variant))
            //        {
            //            if (variant.fixes > 0)
            //                variant.fixes--;
            //            if (m != 0)
            //                HPExtendingDel = true;
            //        }
            //    }
            //}

            if ((tracing >= Tracing.traceChoices && outer) || tracing == Tracing.traceFollowers)
            {
                TraceVariants(m, indentDepth, variantSet, currentMerCount, readProps.merCount);
            }

            int basesLeftInRead = currentMerCount - m;

            //if (!bestChoiceFound && !recursiveCall && longestFollowers < Math.Min(20, basesLeftInRead))
            //{
            //    repairsLeft = maxFollowerRepairs * 2;
            //    int previousLongestFollowers = longestFollowers;

            //    if (repairsLeft < basesLeftInRead / 2)
            //    {
            //        longestFollowers = RefreshFollowerCounts(variantSet, read, m, currentMerCount, readProps, tryHarder, subFixesOnly, repairsLeft, indentDepth, ref thmCalls);

            //        if ((tracing >= Tracing.traceChoices && outer) || tracing == Tracing.traceFollowers)
            //        {
            //            traceClump.Add(Indent(indentDepth) + "recalculated followers: max allFollowers " + previousLongestFollowers + "-->" + longestFollowers);
            //            TraceVariants(m, indentDepth, variantSet, currentMerCount, readProps.merCount);
            //        }
            //    }
            //}

            // only one variant looks possible don't need to choose the best one
            if (variantSet.Count == 1)
            {
                bestChoiceFound = true;
                bestChoiceIdx = 0;
                bestChoiceReason = "only";
            }

            // multiple possible variants so choose the best one
            if (variantSet.Count > 1)
            {
                // try find the 'best' fix for this mer - or decide that it can't be 'fixed'         
                bestChoiceFound = ChooseBestVariant(variantSet, maxFollowerRepairs, merDepth, readProps.minDepth, readProps.OKDepth, pairDepth, readProps.minPairDepth, readProps.OKPairDepth, readProps.unbalancedRead, checkReason, HPExtendingDel,
                                                    out bestChoiceIdx, out bestChoiceReason, out bestChoicePerfect);
            }

            if (bestChoiceFound)
            {
                replacementMer = variantSet[bestChoiceIdx];
            }

            if ((tracing >= Tracing.traceChoices && outer) || tracing >= Tracing.traceRead)
            {
                healMerTimer.Stop();
                float timeTakenToHeal = (float)(healMerTimer.ElapsedMilliseconds) / 1000;

                string choiceMsg;
                if (replacementMer != null)
                {
                    choiceMsg = Indent(indentDepth) + "@" + m + " Choose " + fixNames[(int)replacementMer.fixType] + "(" + replacementMer.lengthDelta + ")" +
                                   " " + kMers.ExpandMer(replacementMer.variant, merSize) + " " +
                                   bestChoiceReason + " sum=" + replacementMer.sum + " fo=" + replacementMer.goodFollowers + "/" +
                                   replacementMer.allFollowers + "/" + replacementMer.maxFollowers + "/" + replacementMer.merCount + " fc@" + replacementMer.mersToFirstChoice + " " +
                                   fixNames[(int)replacementMer.nextFixType] + "@" + replacementMer.mersToNextFix + " fx=" + replacementMer.fixes + (replacementMer.perfectFix ? " (perfect) " : "") +
                                   " with " + repairsLeft + " fixes left. " + checkReason.ToString() + " thm=" + thmCalls;
                }
                else
                    choiceMsg = Indent(indentDepth) + "@" + m + " no viable kMer";

                if (outer || tracing == Tracing.traceFollowers)
                {
                    TraceClumpAdd(choiceMsg);
                    if (outer && replacementMer != null && replacementMer.breadcrumbs != null)
                        TraceClumpAdd(replacementMer.breadcrumbs);
                }
                else
                {
                    if (replacementMer != null)
                        replacementMer.breadcrumbs = replacementMer.breadcrumbs + Environment.NewLine + choiceMsg;
                }
            }

            // now return all the merProperty and other cached objects used here 

            for (int i = 0; i < variantSet.Count; i++)
                if (i != bestChoiceIdx)
                {
                    ReturnMerProperty(variantSet[i]);
                }
            ReturnVariantSet(variantSet);

            return true;
        }

        private static bool ChooseBestVariant(List<kMerProperties> variantSet, int maxFixes, int merDepth, int minDepth, int OKDepth, int pairDepth, int minPairDepth, int OKPairDepth, bool readUnbalanced, kMerState checkReason, bool HPExtendingDel, 
                                              out int bestChoiceIdx, out string bestChoiceReason, out bool bestChoicePerfect)
        {
            int variantsInSet = variantSet.Count;
            int viableVariantsInSet = variantsInSet;
            int chosenIdx = 0;
            bestChoicePerfect = false;

            int highestGoodFollowers;
            int highestAllFollowers;
            int lowestGoodFixes;
            int lowestAllFixes;
            int maxMaxFollowers = 0;
            GetHighestFollowers(variantSet, out highestGoodFollowers, out highestAllFollowers, out lowestGoodFixes, out lowestAllFixes);

            // look for a 'perfect' fix - one that takes us to the max followers - prefer 'good' followers
            int countPerfect = 0;
            int minFixesForPerfect = int.MaxValue;

            bool someVariantBalanced = true;
            bool pairsVariantFound = false;

            // if 'none' gives a perfect result, always choose it over any fix
            bool noneIsPerfect = false;
            // remember if we're looking only for 'good perfect'
            bool foundGoodPerfect = false;

            // see if any variant takes us straight to the end (allowing for downstream fixes)
            for (int v = 0; v < variantsInSet; v++)
            {
                kMerProperties variant = variantSet[v];
                if (!variant.validVariant)
                    continue;

                // be wary of ggg-style Illumina artefacts
                if (variant.fixType == FixTypes.fixNone && variant.unbalanced && MerIsHomopolymer(variant.variant))
                    continue;

                // remember if there were any balanced variants on this pass as well
                someVariantBalanced = someVariantBalanced & !variant.unbalanced;
                // and if any of them have OK pair depths
                pairsVariantFound |= variant.pairDepth >= minPairDepth; 
                // and the largest maxFollowers value
                if (variant.maxFollowers > maxMaxFollowers)
                    maxMaxFollowers = variant.maxFollowers;

                // is this variant a 'good perfect'?
                if (variant.maxFollowers == variant.goodFollowers)
                {
                    // if we haven't found any 'good perfect' yet, reset the counts etc in case we've come across some 'perfect all's
                    if (!foundGoodPerfect && countPerfect > 0)
                    {
                        countPerfect = 0;
                        minFixesForPerfect = int.MaxValue;
                        noneIsPerfect = false;

                        for (int rv = 0; rv < v; rv++)
                            if (variantSet[v].fixType != FixTypes.fixNone)
                                variantSet[rv].perfectFix = false;
                    }

                    foundGoodPerfect = true;
                    countPerfect++;
                    chosenIdx = v;
                    variant.perfectFix = true;
                    if (variant.fixes < minFixesForPerfect)        
                        minFixesForPerfect = variant.fixes;
                    if (variant.fixType == FixTypes.fixNone)
                        noneIsPerfect = true;
                }

                if (!foundGoodPerfect && variant.maxFollowers == variant.allFollowers)
                {
                    countPerfect++;
                    variant.perfectFix = true;
                    chosenIdx = v;

                    if (variant.fixType == FixTypes.fixNone)
                        noneIsPerfect = true;
                    if (variant.fixes < minFixesForPerfect)         
                        minFixesForPerfect = variant.fixes;
                }
            }

            // None can get us to the end of the read - always prefer to leave a kMer unchanged... 
            if (noneIsPerfect)
            {
                bestChoiceIdx = 0;
                bestChoiceReason = "perfectNone";
                bestChoicePerfect = true;
                return true;
            }

            // just one perfect choice - chosenIdx will be correct in this case
            if (countPerfect == 1)
            {
                bestChoiceIdx = chosenIdx;
                bestChoiceReason = "perfect";
                bestChoicePerfect = true;
                return true;
            }

            // multiple perfect so see if only one has the fewest fixes
            if (countPerfect > 1)
            {
                countPerfect = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];

                    if (!variant.perfectFix)
                        continue;

                    if (variant.fixes == minFixesForPerfect)
                    {
                        countPerfect++;
                        chosenIdx = v;
                    }
                }
            }

            // exactly one 'perfect' result with the fewest fixes
            if (countPerfect == 1)
            {
                bestChoiceIdx = chosenIdx;
                bestChoiceReason = "perfectFixes";
                bestChoicePerfect = true;
                return true;
            }
            
            // remove any non-paired variants if there are any paired ones
            if (pairsVariantFound)
            {
                foreach (kMerProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;

                    if (variant.pairDepth != -1 && variant.pairDepth < minPairDepth)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                }
            }

            // no single 'perfect' variant so clean up the variants - balanced depths with an accepable pair depth is preferred

            // get the best balanced depth from all of the variants
            int bestVariantDepth = 0;
            foreach (kMerProperties variant in variantSet)
            {
                if (!variant.validVariant)
                    continue;

                if (variant.depth > bestVariantDepth && !(variant.unbalanced && !readUnbalanced))
                    bestVariantDepth = variant.depth;
            }

            // if there is at least one depth-viable balanced read, cull any unbalanced variants (including the starting k-mer)
            if (bestVariantDepth >= minDepth)
            {
                foreach (kMerProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;

                    //if (variant.fixType == FixTypes.fixNone)
                    //    continue;

                    if (variant.unbalanced && !readUnbalanced)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                }
            }

            // are any of the variants 'perfect'? If so, delete all the non-perfect ones, and if there's only one left, choose it
            bool perfectFound = false;
            if (countPerfect > 0)
                foreach (kMerProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;

                    perfectFound = perfectFound | variant.perfectFix;
                }

            if (perfectFound)
            {
                countPerfect = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];

                    if (!variant.validVariant)
                        continue;

                    if (!variant.perfectFix)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                    else
                    {
                        countPerfect++;
                        chosenIdx = v;
                    }
                }

                if (countPerfect == 1)
                {
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "perfNext";
                    bestChoicePerfect = false;
                    return true;
                }
            }

            // calculate the 'pretty much the same' margin for follower counts
            int followersMargin = highestAllFollowers / 10;
            if (followersMargin == 0)
                followersMargin = 1;
            // and now the threshold for about-equally-as-good
            int goodEnoughFollowers = highestAllFollowers - followersMargin;
            if (goodEnoughFollowers < 1)
                goodEnoughFollowers = 1;

            // can we find the highest followers only once? (within the margin)
            bool foundFollowersTwice = false;
            bool foundFollowers = false;
            for (int v = 0; v < variantsInSet; v++)
            {
                kMerProperties variant = variantSet[v];
                if (!variant.validVariant)
                    continue;

                int followers = variant.allFollowers;
                if (followers >= goodEnoughFollowers)
                {
                    if (foundFollowers)
                        foundFollowersTwice = true;
                    foundFollowers = true;
                    chosenIdx = v;
                }
                else
                    variant.markedVariant = true;
            }

            if (foundFollowers && !foundFollowersTwice)
            {
                bestChoiceIdx = chosenIdx;
                bestChoiceReason = "followers";
                return true;
            }

            // found multiple 'good enough followers' variants, so remove any others from the list
            foreach (kMerProperties variant in variantSet)
                if (foundFollowersTwice)
                {
                    if (variant.markedVariant && variant.validVariant)
                    {
                        variant.validVariant = false;
                        viableVariantsInSet--;
                    }
                }
                else
                    variant.markedVariant = false;

            // multiple 'any' followers, can we pick a winner on the 'good' followers alone
            if (foundFollowersTwice)
            {
                bool foundInGoodList = false;
                bool foundinGoodListTwice = false;
                goodEnoughFollowers = highestGoodFollowers - followersMargin;
                if (goodEnoughFollowers < 1)
                    goodEnoughFollowers = 1;

                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;
                    if (variant.goodFollowers >= goodEnoughFollowers)
                    {
                        if (foundInGoodList)
                            foundinGoodListTwice = true;
                        foundInGoodList = true;
                        chosenIdx = v;
                    }
                    else
                        variant.markedVariant = true;               // mark for possible culling
                }

                if (foundInGoodList && !foundinGoodListTwice)
                {
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "good";
                    return true;
                }

                // found multiple goood-enough 'good' variants, so remove any others from the list
                foreach (kMerProperties variant in variantSet)
                    if (foundinGoodListTwice)
                    {
                        if (variant.markedVariant && variant.validVariant)
                        {
                            variant.validVariant = false;
                            viableVariantsInSet--;
                        }
                    }
                    // and reset mark if we're not culling
                    else
                        variant.markedVariant = false;
            }

            // all still-viable variants have the similar followers so choose None if it's still viable
            if (viableVariantsInSet > 0 && variantSet[0].validVariant && variantSet[0].fixType == FixTypes.fixNone)
            {
                bestChoiceIdx = 0;
                bestChoiceReason = "None";
                return true;
            }

            // no best choice by looking at followers... try for the best gap-fixes metric amongst the remaining candidates
            if (viableVariantsInSet > 0)
            {
                int bestGapsFixes = 0;

                foreach (kMerProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;
                    // replace mersToNextFix with combined metric
                    variant.mersToNextFix = variant.mersToNextFix + (maxFixes - variant.fixes);
                    if (variant.mersToNextFix > bestGapsFixes)
                        bestGapsFixes = variant.mersToNextFix;
                }

                // how many valid variants have this best gap/fixes?
                int foundTarget = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.mersToNextFix == bestGapsFixes)
                    {
                        chosenIdx = v;
                        foundTarget++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }
                }

                // and if there's only one variant with the best gap/fixes, choose that one
                if (foundTarget == 1)
                {
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "gaps";
                    return true;
                }
                else
                    foreach (kMerProperties variant in variantSet)
                        if (foundTarget > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // multiple choices with best followers and best gaps and fewest fixes, try for best sum if one is clearly better than the others
            if (viableVariantsInSet > 0)
            {
                int summedSums = 0;
                foreach (kMerProperties variant in variantSet)
                {
                    if (!variant.validVariant)
                        continue;
                    summedSums += variant.sum;
                }

                int foundTarget = 0;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.sum >= summedSums * 7 / 10)
                    {
                        chosenIdx = v;
                        foundTarget++;
                    }
                    else
                    {
                        variant.markedVariant = true;
                    }
                }

                // and if there's only one variant with the high-enough sum, choose that one
                if (foundTarget == 1)
                {
                    bestChoiceIdx = chosenIdx;
                    bestChoiceReason = "sums";
                    return true;
                }
                else
                    foreach (kMerProperties variant in variantSet)
                        if (foundTarget > 1)
                        {
                            if (variant.markedVariant && variant.validVariant)
                            {
                                variant.validVariant = false;
                                viableVariantsInSet--;
                            }
                        }
                        else
                            variant.markedVariant = false;
            }

            // multiple choices with best followers, best gaps, best fixes and near-equal depth-sums
            // try for status quo if this is a possible choice
            if (variantSet[0].validVariant && variantSet[0].fixType == FixTypes.fixNone && !HPExtendingDel)
            {
                // unchanged mer will always be the first in the lists if it was viable 
                if (merDepth > minDepth && (pairDepth == -1 || pairDepth >= minPairDepth))  
                {
                    bestChoiceIdx = 0;
                    bestChoiceReason = "none";
                    return true;
                }
            }

            // make an almost arbitrary decision - preferring del/ins over sub and del over ins to avoid rewriting too much
            // but preferring sub if we're at the end of a read or we are correcting Illumina data
            if (viableVariantsInSet > 0 && errorModel == ErrorModel.mostlySubs)
            {
                int bestSubSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    // find the Sub variant with the highest sum (if there is one)
                    if (variant.fixType == FixTypes.fixSub && variant.sum > bestSubSum)
                    {
                        bestSubSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceReason = "pSub";
                    return true;
                }
            }

            if (viableVariantsInSet > 0)
            {
                int bestDelSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixDel && variant.sum > bestDelSum)
                    {
                        bestDelSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceReason = "pDel";
                    return true;
                }
            }

            if (viableVariantsInSet > 0)
            {
                int bestInsSum = 0;
                bestChoiceIdx = -1;
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixIns && variant.sum > bestInsSum)
                    {
                        bestInsSum = variant.sum;
                        bestChoiceIdx = v;
                    }
                }
                if (bestChoiceIdx >= 0)
                {
                    bestChoiceReason = "pIns";
                    return true;
                }
            }

            if (viableVariantsInSet > 0)
            {
                // no del or ins in list but list isn't empty 
                for (int v = 0; v < variantsInSet; v++)
                {
                    kMerProperties variant = variantSet[v];
                    if (!variant.validVariant)
                        continue;

                    if (variant.fixType == FixTypes.fixSub)
                    {
                        bestChoiceIdx = v;
                        bestChoiceReason = "fSub";
                        return true;
                    }
                }
            }

            // couldn't find a suitable replacement
            bestChoiceIdx = -1;
            bestChoiceReason = "null";
            return false;
        }

        private static void GetHighestFollowers(List<kMerProperties> variants, out int highestGoodFollowers, out int highestAllFollowers,
                                                out int lowestGoodFixes, out int lowestAllFixes)
        {
            highestGoodFollowers = 0;
            highestAllFollowers = 0;
            lowestGoodFixes = int.MaxValue;
            lowestAllFixes = int.MaxValue;

            for (int v = 0; v < variants.Count; v++)
            {
                kMerProperties variant = variants[v];
                if (!variant.validVariant)
                    continue;
                int goodFollowers = variant.goodFollowers;
                if (goodFollowers > highestGoodFollowers)
                {
                    highestGoodFollowers = goodFollowers;
                    lowestGoodFixes = variant.fixes;
                }
                else
                {
                    if (goodFollowers == highestGoodFollowers && variant.fixes < lowestGoodFixes)
                        lowestGoodFixes = variant.fixes;
                }

                int allFollowers = variant.allFollowers;
                if (allFollowers > highestAllFollowers)
                {
                    highestAllFollowers = allFollowers;
                    lowestAllFixes = variant.fixes;
                }
                else
                {
                    if (allFollowers == highestAllFollowers && variant.fixes < lowestAllFixes)
                        lowestAllFixes = variant.fixes;
                }

            }
        }

        private static void TraceVariants(int m, int indentDepth, List<kMerProperties> variantSet, int merCount, int readMerCount)
        {
            for (int v = 0; v < variantSet.Count; v++)
                if (variantSet[v].validVariant)
                    TraceClumpAdd(Indent(indentDepth) + "@" + m + " " +
                                "#" + v + " " + fixNames[(int)variantSet[v].fixType] + "(" + variantSet[v].lengthDelta + ")" +
                                " s=" + variantSet[v].sum + " fo=" + variantSet[v].goodFollowers + "/" + variantSet[v].allFollowers + "/" + variantSet[v].maxFollowers + "/" + variantSet[v].merCount +
                               " fx=" + variantSet[v].fixes + " nf=" + fixNames[(int)variantSet[v].nextFixType] + "+" + variantSet[v].mersToNextFix + " fc+" + variantSet[v].mersToFirstChoice +
                               " mc=" + merCount + " rmc=" + readMerCount + " " + (variantSet[v].perfectFix ? "(perfect)" : "") + " " +
                               kMers.ExpandMer(variantSet[v].variant, merSize));
        }

        // Counts the followers for a kMer - possibly temporarily making a few downstream repairs if needed to overcome nearly adjacent errors (recursively calls TryHealingMer).
        // 
        private static bool CountFollowers(bool subFixesOnly, kMerProperties variantMer, int m, ReadProperties readProps, int currentMerCount, 
                                           int repairsLeft, int indentDepth, ref int thmCalls)
        {
            Stopwatch countFollowersTimer = null;
            if (tracing > Tracing.traceOff)
            {
                countFollowersTimer = new Stopwatch();
                countFollowersTimer.Start();
            }

            FixTypes firstFixType = FixTypes.fixNone;
            int mersToFirstFix = 0;
            int mersToFirstChoice = 0;
            int goodFollowerCount = 0;
            int allFollowerCount = 0;
            int healedMerCount = Math.Min(currentMerCount, readProps.merCount);               // use this count if the follower loop terminates here
            int fixes = 0;
            int followerSum = 0;
            bool perfectFix = false;
            string breadcrumbs = null;
            Sequence kMerContext = variantMer.fixContext;
            FixTypes previousFixType = variantMer.fixType;

            kMerProperties replacementMer = null;
            bool thmOK = true;

            int previousMerDepth = variantMer.depth;
            int previousPairDepth = variantMer.pairDepth;

            ulong nextMer = variantMer.variant;
            int nextM = m + 1;
            char followingBase;

            while (GenerateNextMer(nextMer, kMerContext, nextM, currentMerCount, readProps.OKDepth, readProps.OKPairDepth,  out nextMer, out followingBase))   // generate the next mer 
            {
                string watchMer = null;
                if (tracing >= Tracing.traceRead)
                    watchMer = kMers.ExpandMer(nextMer, merSize);

                kMerState checkReason = kMerState.merOK;
                bool unbalancedMer;
                bool tiltedMer;
                int nextMerDepth = kMersTable.GetDepthSum(nextMer, readProps.minDepth, out unbalancedMer, out tiltedMer);
#if STATS
                Interlocked.Increment(ref gds4Calls_Next);
#endif
                int nextPairDepth = GetPairDepth(kMerContext, nextMer, readProps.mers, readProps.depths, nextM, readProps.minDepth);

                // if the next mer (current mer + next base from read) looks OK, keep going 
                if (!MerIsHealingCandidate(kMerContext, readProps, nextM, nextMer, unbalancedMer, tiltedMer, followingBase, nextMerDepth, previousMerDepth, nextPairDepth, previousPairDepth, false, out checkReason))   
                {
                    // checkReason must be merOK if we get here
                    goodFollowerCount++;
                    allFollowerCount++;
                    followerSum += nextMerDepth;
                    mersToFirstFix++;
                    mersToFirstChoice++;
                    nextM++;
                    previousMerDepth = nextMerDepth;
                    continue;
                }

                //// reduce the number of repairs allowed if we're not just checking for an improvement
                //if (checkReason == kMerState.merBad)
                //    repairsLeft--;

                // ran out of good followers (or we've hit an HP), perhaps the mer can/should be repaired and we can keep going.
                // Only a few repairs are allowed to stop us just rewriting the read from the consensus (and to limit the tree exploration depth)
                if (repairsLeft > 0)
                {
                    thmOK = TryHealingMer(checkReason, subFixesOnly, true, nextMer, nextMerDepth, nextPairDepth, previousMerDepth, previousPairDepth, unbalancedMer,
                                              kMerContext, readProps, nextM, currentMerCount, repairsLeft, previousFixType, mersToFirstChoice, out replacementMer, false, indentDepth + 1, ref thmCalls);

                    if (replacementMer == null || !thmOK)
                        break;

                    bool merWasChanged = replacementMer.fixType != FixTypes.fixNone;

                    // did our attempt to repair a mer end up with the same mer as we started with?
                    if (replacementMer.variant == variantMer.variant && merWasChanged)
                    {
                        if (tracing >= Tracing.traceRead)
                            TraceClumpAdd("Killing follower path at potential loop @" + nextM + " " +
                                            kMers.ExpandMer(nextMer, merSize));
                        nextMerDepth = 0;
                        replacementMer.allFollowers = 0;
                        replacementMer.goodFollowers = 0;
                        replacementMer.maxFollowers = 0;
                        replacementMer.fixes = 99;
                        replacementMer.mersToNextFix = 0;
                        replacementMer.mersToFirstChoice = 0;
                        replacementMer.perfectFix = false;
                        replacementMer.fixType = FixTypes.fixAbandon;
                    }

                    // and add in the downstream followers to the follower counts we got from the loop prior to calling TryHealingMer
                    goodFollowerCount += replacementMer.goodFollowers;
                    allFollowerCount += replacementMer.allFollowers;
                    fixes = replacementMer.fixes;
                    if (merWasChanged)
                    {
                        firstFixType = replacementMer.fixType;
                    }
                    else
                    {
                        mersToFirstFix += replacementMer.mersToNextFix + 1;
                        firstFixType = replacementMer.nextFixType;
                    }
                    followerSum += replacementMer.sum;
                    healedMerCount = replacementMer.merCount;
                    breadcrumbs = replacementMer.breadcrumbs;

                    // no need to continue with this loop after trying a fix - recursion will have done it for us
                    break;
                }

                // not a good mer and no more fixes allowed
                firstFixType = FixTypes.fixAbandon;
                if (tracing >= Tracing.traceRead)
                    breadcrumbs = Indent(indentDepth+1) + "abandoned @" + nextM + " too many fixes";
                break;
            } // GenerateNextMer loop

            if (replacementMer != null)
            {
                ReturnMerProperty(replacementMer);
            }

            // populate the 'downstream' field in the mer object passed in as a parameter
            variantMer.goodFollowers = goodFollowerCount;
            variantMer.allFollowers = allFollowerCount;
            variantMer.merCount = healedMerCount;
            variantMer.fixes = fixes;
            variantMer.sum = followerSum;
            variantMer.mersToNextFix = mersToFirstFix;
            variantMer.mersToFirstChoice = mersToFirstChoice;
            variantMer.nextFixType = firstFixType;
            variantMer.perfectFix = perfectFix;
            variantMer.breadcrumbs = breadcrumbs; 

            if (tracing == Tracing.traceFollowers)
                TraceFollowerResults(indentDepth, m, variantMer, countFollowersTimer);

            return thmOK;
        }

        private static void TraceFollowerResults(int indentDepth, int m, kMerProperties variantMer, Stopwatch countFollowersTimer)
        {
            countFollowersTimer.Stop();
            float timeTakenToCount = (float)(countFollowersTimer.ElapsedMilliseconds) / 1000;

            TraceClumpAdd(Indent(indentDepth) + "counted followers " + variantMer.goodFollowers + "/" + variantMer.allFollowers + 
                " " + variantMer.fixes + " fixes" +
                " for " + fixNames[(int)variantMer.fixType] + " @" + m + " in " + timeTakenToCount.ToString("F3"));
        }

        private static int UpdateRead(FixTypes fixType, Sequence healingRead, int m, ulong replacementMer, int lengthDelta, int merCount, int initialMerCount)
        {
            if (fixType == FixTypes.fixDel)
            {
                // the replacement 'del' mer has had a base inserted somewhere. This means that we have to insert 
                // a placeholder base in the read string, prior to replacing the healed mer  
                // This can increase the length of a read pass its initial size but we don't truncate here in case
                // we need to redeem these extra bases after a later Ins has deleted bases.
                healingRead.Insert(m + merSize - 1, 'X');
                merCount++;
                if (merCount > initialMerCount)
                    merCount = initialMerCount;
                healingRead.Replace(m, replacementMer, merSize);
            }

            if (fixType == FixTypes.fixSub)
            {
                // for substitutions, the structure of the healed read doesn't change, 
                // so just overwrite the part of the read corresponding to the mer being healed
                healingRead.Replace(m, replacementMer, merSize);
            }

            if (fixType == FixTypes.fixIns)
            {
                // the replacement 'ins' mer has had a base deleted somewhere. This means that we have to delete 
                // a base from the read string, prior to replacing the healed mer (which includes the following base from the read)
                merCount += lengthDelta;
                for (int i = lengthDelta; i < 0; i++)
                    healingRead.Remove(m);
                healingRead.Replace(m, replacementMer, merSize);
            }

            return merCount;
        }

        private static void UpdateQuals(FixTypes fixType, Sequence quals, int m, int lengthDelta)
        {
            // fixing first mer in the read, so set all of its quals before making any other possible changes
            if (m == 0)
            {
                for (int i = 0; i < merSize; i++)
                    quals.Bases[i] = replacementQual;
            }

            // no change but mismatch between merDepth and qual score
            if (fixType == FixTypes.fixNone)
            {
                quals.Bases[m + merSize - 1] = replacementQual;
            }

            // fixing a Sub so just replace the qual
            if (fixType == FixTypes.fixSub)
            {
                quals.Bases[m + merSize - 1] = replacementQual;
            }
            // repair added a base so we need to shift remaining qual values right by one place
            if (fixType == FixTypes.fixDel)
            {
                quals.Insert(m + merSize - 1, replacementQual);
            }

            // repair deleted bases so we need to shift remaining qual values left by as many places as needed 
            if (fixType == FixTypes.fixIns)
            {
                for (int i = lengthDelta; i < 0; i++)
                    quals.Remove(m + merSize - 1);
            }
        }

        private static void UpdateMersDepths(FixTypes fixType, Sequence read, int m, ReadProperties readProps, ulong replacementMer, int replacementDepth, bool replacementBalanced)
        {
            // no change but mismatch between merDepth and qual score
            if (fixType == FixTypes.fixNone)
                return;

            if (fixType == FixTypes.fixSub)
            {
                // if we're at the start of the read, it's possible that any base in the first kMer was changed,
                // and this could affect the following kMers that contain the change. In this case we'll rescan the 
                // first MerSize kMers
                if (m == 0)
                {
                    ResetSetMerDepths(read, merSize, readProps);
                }
                else
                // fixing a sub in the middle of the read - just replace this kMer and depth
                {
                    //if (tracingThisRead)
                    //    TraceClumpAdd("replacing " + kMers.ExpandMer(mers[m], merSize) + " with " + kMers.ExpandMer(replacementMer, merSize) + " @" + m);
                    readProps.mers[m] = replacementMer;
                    readProps.depths[m] = replacementDepth;
                    if (replacementBalanced && !readProps.balancedDepths[m])
                        readProps.balancedMerCount++;
                    readProps.balancedDepths[m] = replacementBalanced;
                }
            }

            // repair added or deleted bases so we need to recalculate the kMers and depths for the tail of the read (post m)
            if (fixType == FixTypes.fixDel || fixType == FixTypes.fixIns)
            {
                int mersInRead = read.Length - merSize + 1;

                if (mersInRead > readProps.mers.Length)
                {
                    Array.Resize<ulong>(ref readProps.mers, mersInRead + 100);
                    Array.Resize<int>(ref readProps.depths, mersInRead + 100);
                    Array.Resize<int>(ref readProps.pairDepths, mersInRead + 100);      // need extra space for debug use
                    Array.Resize<int>(ref readProps.merChangeCost, mersInRead + 100);
                    Array.Resize<bool>(ref readProps.zeroStrand, mersInRead + 100);
                    Array.Resize<bool>(ref readProps.balancedDepths, mersInRead + 100);
                }

                ResetSetMerDepths(read, readProps.merCount, readProps);
            }
        }

        private static int CheckForRewriting(int m, int[] costs, bool inTail)
        {
            int startOfFixScan = m - rewriteRegion;
            if (startOfFixScan < 0)
                startOfFixScan = 0;

            int summedCost = costs[startOfFixScan];
            int basesSinceLastChange = summedCost == 0 ? 1 : 0;

            for (int i = startOfFixScan + 1; i <= m; i++)
            {
                int nextCost = costs[i];

                if (nextCost == 0)
                {
                    basesSinceLastChange++;
                    continue;
                }

                // reduce effective cost if there's been a gap between changes
                if (basesSinceLastChange > goodMerRunLength)
                    nextCost--;
                else
                {
                    // and penalise consecutive changes unless we're rewriting a tail
                    if (!inTail)
                       nextCost += gapCosts[basesSinceLastChange];
                }

                basesSinceLastChange = 0;

                summedCost += nextCost;
            }

            return summedCost;
        }

        private static int RetreatToGoodMers(int m, int[] merChangeCost)
        {
            int startOfFixScan = m - rewriteRegion;
            if (startOfFixScan < 0)
                startOfFixScan = 0;

            // find last long run of good kMers (at least goodMerRunLength long)
            int lastGoodRunStart = startOfFixScan;
            int lastGoodRunLength = goodMerRunLength;
            int lastGoodEnd = startOfFixScan;
            int lastRunLength = -1;
            int lastRunStart = -1;
            for (int i = startOfFixScan; i < m; i++)
            {
                if (merChangeCost[i] == 0)
                {
                    lastRunLength++;
                    if (lastRunStart == -1)
                        lastRunStart = i;
                }
                else
                {
                    if (lastRunLength >= goodMerRunLength)
                    {
                        lastGoodRunStart = lastRunStart;
                        lastGoodRunLength = lastRunLength;
                        lastGoodEnd = i - 1;
                    }
                    lastRunLength = 0;
                    lastRunStart = -1;
                }
            }
            // adjust if we're in a good run when we hit the end
            if (lastRunLength >= goodMerRunLength)
            {
                lastGoodRunStart = lastRunStart;
                lastGoodRunLength = lastRunLength;
                lastGoodEnd = m - 1;
            }

            return lastGoodEnd;
        }


        private static void TrimAbandonedRead(Sequence read, Sequence quals, ReadProperties readProps)
        {
            int m = readProps.abandonedAtM;
            bool healingAbandonedAtLC = false;

            // could have run into trouble as a result of a stutter or similar batch of unbalanced kMers
            // so trim back until we get good looking data
            while (m > 0)
            {
                ulong currentMer;
                bool depthUnbalanced = false;
                bool foundLCMer = false;

                bool merValid = Sequence.CondenseMer(read, m, merSize, out currentMer);
                if (merValid)
                {
                    bool tilted;
                    int depthSum = kMersTable.GetDepthSum(currentMer, readProps.minDepth, out depthUnbalanced, out tilted);
#if STATS
                    Interlocked.Increment(ref gds4Calls_TrimAbs);
#endif
                    foundLCMer = depthUnbalanced && lowComplexityTrap.Contains(currentMer & lowComplexityMask);
                    healingAbandonedAtLC |= foundLCMer;
                }

                if ((!readProps.unbalancedRead && depthUnbalanced) || foundLCMer || !merValid)
                    m--;
                else
                    break;
            }

            read.Length = m + merSize - 1;                      // truncate the read at this point
            if (quals.Length > 0)
                quals.Length = m + merSize - 1;

            readProps.merCount = m;
            if (pairsTable != null)
                readProps.pairsCount = read.Length - pairsTable.pairFullLength + 1;
            if (readProps.pairsCount < 0)
                readProps.pairsCount = 0;

            CheckReadForProblems(read, readProps);
        }

        private static bool CollectPlausibleVariants(kMerState checkReason, bool subFixesOnly, int repairsLeft, List<kMerProperties> variantSet, 
                                                    ReadProperties readProps, int m, int currentMerCount, ulong kMer, int merDepth, int pairDepth, bool unbalancedMer, Sequence fixContext, int previousMerDepth, int previousPairDepth,
                                                    FixTypes previousFixType, int mersToFirstChoice, int indentDepth, ref int thmCalls )
        {
            // always add the unchanged mer as the first variant if it's viable
            if (checkReason != kMerState.merBad)
            {
                kMerProperties startingMer = AllocateMerProperty();
                startingMer.variant = kMer;
                startingMer.fixType = FixTypes.fixNone;
                startingMer.depth = merDepth;
                startingMer.unbalanced = unbalancedMer;
                startingMer.pairDepth = pairDepth;
                startingMer.fixContext = fixContext;

                bool cfOK = CountFollowers(subFixesOnly, startingMer, m, readProps, currentMerCount, repairsLeft, indentDepth, ref thmCalls);

                if (!cfOK)
                {
                    ReturnMerProperty(startingMer);
                    return false;
                }

                // if the starting mer has a viable score, record it as producing a follower as well
                if (startingMer.depth >= readProps.OKDepth && (startingMer.pairDepth == -1 || startingMer.pairDepth >= readProps.OKPairDepth)) 
                    startingMer.goodFollowers++;
                if (startingMer.depth >= readProps.minDepth || startingMer.pairDepth >= readProps.minPairDepth) 
                    startingMer.allFollowers++;
                startingMer.sum += startingMer.depth;
                // calculate highest possible follower count (inclusive of current base ['m' is index --> adds 1 to count])
                startingMer.maxFollowers = startingMer.merCount - m;

                variantSet.Add(startingMer);

                if (tracing >= Tracing.traceRead)
                {
                    string traceMsg = Indent(indentDepth) + "None @" + m + " " + kMers.ExpandMer(kMer, merSize) +
                                    " sum=" + startingMer.sum + " fo=" + startingMer.goodFollowers + "/" + startingMer.allFollowers + "/" + startingMer.maxFollowers + "/" + startingMer.merCount +
                                    " fx=" + startingMer.fixes +
                                    " with " + repairsLeft + " fixes left";
                    if (tracing == Tracing.traceFollowers)
                        TraceClumpAdd(traceMsg);
                    startingMer.breadcrumbs = traceMsg + Environment.NewLine + startingMer.breadcrumbs;
                }
            }

            if (tracing == Tracing.traceFollowers)
                TraceClumpAdd(Indent(indentDepth) + "Finding variants for " + kMers.ExpandMer(kMer, merSize) + " @" + m + " with " + repairsLeft + " fixes left");

            // add best (perhaps tied) variants for each repair type to the set of candidate mers/fix types  
            // and stop if FindPlausibleVariants fails (meaning THM has hit the tree limit)
            bool fpvOK = true;

            // del+sub == sub+del; ins+sub == sub+ins - so don't generate these variants
            if (subFixesOnly || (mersToFirstChoice == 0 && !(previousFixType == FixTypes.fixDel || previousFixType == FixTypes.fixIns)))
            {
                int firstSubIdx = variantSet.Count; 
                int subVars; 
                fpvOK = FindPlausibleVariants(checkReason, subFixesOnly, FixTypes.fixSub, VariantWanted.varyAnyOne, kMer, unbalancedMer, previousMerDepth, previousPairDepth, fixContext, readProps, m, currentMerCount,
                                                    repairsLeft - 1, variantSet, out subVars, indentDepth, ref thmCalls);

                // if we have a particularly bad start (no plausible variants or original mer) try a bit harder to catch multiple (Illumina) sub errors in first mer
                // (we could try to catch these on the reverse pass, but there's a chance with shorter Illumina data that we'll not find any good mer to start from...)
                // also try harder if we're at the start of the read and all the single variants look to be poor
                if (fpvOK && m == 0)
                {
                    // if there were any single subVars, check them to make sure none of them are 'poor'
                    bool allSubVarsOK = true;
                    if (subVars > 0)
                    {
                        for (int i = firstSubIdx; i < variantSet.Count; i++)
                        {
                            // say a subVar is poor if we have to immediately try to fix the next kMer as well
                            if (variantSet[i].mersToNextFix == 0 || variantSet[i].unbalanced)
                            {
                                allSubVarsOK = false;
                                break;
                            }
                        }
                        // remove the single-sub variants (they'll be included in the double-sub)
                        if (!allSubVarsOK)
                        {
                            // leave the original kMer (if present) and delete the subs
                            for (int v = variantSet.Count-1; v >= firstSubIdx; v--)
                            {
                                ReturnMerProperty(variantSet[v]);
                                variantSet.RemoveAt(v);
                            }
                        }
                    }
                    if (subVars == 0 || !allSubVarsOK)
                        fpvOK = FindPlausibleVariants(checkReason, subFixesOnly, FixTypes.fixSub, VariantWanted.varyAnyTwo, kMer, unbalancedMer, previousMerDepth, previousPairDepth, fixContext, readProps, m, currentMerCount,
                                                        repairsLeft - 1, variantSet, out subVars, indentDepth, ref thmCalls);
                }
            }

            // if we already have a 'perfect' fix from only considering 'subs' then don't bother generating and checking the Del and Ins variants as these will not be chosen anyway
            foreach (kMerProperties variant in variantSet)
            {
                if (variant.maxFollowers == variant.goodFollowers)
                    subFixesOnly = true;
            }

            if (!subFixesOnly)
            {
                if (previousFixType != FixTypes.fixIns && mersToFirstChoice == 0 && fpvOK)
                {
                    int delVars;
                    fpvOK = FindPlausibleVariants(checkReason, subFixesOnly, FixTypes.fixDel, VariantWanted.varyAnyOne, kMer, unbalancedMer, previousMerDepth, previousPairDepth, fixContext, readProps, m, currentMerCount,
                                    repairsLeft - 1, variantSet, out delVars, indentDepth, ref thmCalls);
                }
                if (!(previousFixType == FixTypes.fixDel || previousFixType == FixTypes.fixIns) && mersToFirstChoice == 0 && fpvOK)
                {
                    int insVars;
                    fpvOK = FindPlausibleVariants(checkReason, subFixesOnly, FixTypes.fixIns, VariantWanted.varyAnyOne, kMer, unbalancedMer, previousMerDepth, previousPairDepth, fixContext, readProps, m, currentMerCount,
                                    repairsLeft - 1, variantSet, out insVars, indentDepth, ref thmCalls);
                }
            }

            int longestFollowers = 0;
            foreach (kMerProperties variant in variantSet)
            {
                // longest follower 
                if (variant.allFollowers > longestFollowers)
                    longestFollowers = variant.allFollowers;
            }

            return fpvOK;
        }

        private static bool FindPlausibleVariants(kMerState checkReason, bool subFixesOnly, FixTypes variantType, VariantWanted variantsWanted, 
                                                 ulong startingMer, bool unbalancedMer, int previousMerDepth, int previousPairDepth, Sequence currentRead, ReadProperties readProps, int m, int currentMerCount,
                                                 int repairsLeft, List<kMerProperties> variantMers, out int variantsAdded, int indentDepth, ref int thmCalls)
        {
            variantsAdded = 0;                      // count of variants added here (return value)

            if (m > 0)                              // only vary last base once past the first mer
                variantsWanted = VariantWanted.varyLast;

            int noMerVariants = 0;
            MerVariant[] merVariants = AllocateMerVariants(merSize * 4);         
            string[] merVariantsTrace = null;

            // generate all the variants of the current mer (including current mer itself)
            if (variantType == FixTypes.fixSub)
            {
                noMerVariants = GenerateMerSubVariants(startingMer, merVariants, 0, variantsWanted);
                if (variantsWanted == VariantWanted.varyAnyTwo)
                {
                    // collected variants of variants (unsorted, with duplicates) 
                    MerVariant[] doubleSubVariants = AllocateMerVariants(noMerVariants * merSize * 4);
                    int allIdx = 0;
                    int totalDoubleVariants = 0;

                    for (int v = 0; v < noMerVariants; v++)
                    {
                        // generate variants of next variant
                        int noNewSubVariants = GenerateMerSubVariants(merVariants[v].variantMer, doubleSubVariants, allIdx, variantsWanted);
                        totalDoubleVariants += noNewSubVariants;
                        allIdx += noNewSubVariants;
                    }

                    noMerVariants = totalDoubleVariants;
                    ReturnMerVariants(merVariants);
                    merVariants = doubleSubVariants;
                }
            }

            if (variantType == FixTypes.fixDel)
            {
                noMerVariants = GenerateMerDelVariants(startingMer, merVariants, variantsWanted);
            }

            if (variantType == FixTypes.fixIns)
            {
                int maxInslength = currentMerCount - m - 2;     // -2 because we don't want to go to the end of the read
                if (maxInslength > maxGap)
                    maxInslength = maxGap;
                // only generate Ins variants if there is at least one trailing base in the read
                if (maxInslength > 0)
                {
                    noMerVariants = GenerateMerInsVariants(startingMer, currentRead, m + merSize, maxInslength, merVariants, variantsWanted);
                }          
            }

            if (tracing >= Tracing.traceChoices)
            {
                merVariantsTrace = new string[noMerVariants];
                for (int v = 0; v < noMerVariants; v++)
                    merVariantsTrace[v] = kMers.ExpandMer(merVariants[v].variantMer, merSize);
            }

            // look up the scores (and followers) of each variant and add them to the return list if they meet the requirements
            Depth[] merVariantDepths = AllocateDepths(noMerVariants);
            int deepestVariantDepth = 0;

            // get the depths of each variant (and the deepest depth)
            for (int v = 0; v < noMerVariants; v++)
            {
                ulong variantMer = merVariants[v].variantMer;
                bool unbalanced;
                bool tilted;
                int sumDepth = kMersTable.GetDepthSum(variantMer, readProps.minDepth, out unbalanced, out tilted);
#if STATS
                Interlocked.Increment(ref gds4Calls_FPV);
#endif
                merVariantDepths[v].depth = sumDepth;
                merVariantDepths[v].unbalanced = unbalanced;
                if (sumDepth > deepestVariantDepth)
                    deepestVariantDepth = sumDepth;
            }

            // condense the list of variants by copying any likely-looking alternatives towards the start of the merVariants array,
            // discarding 'poor' alternatives (and any ...XXXXXXXXXXXX (long HP) alternatives as these attractors cause performance problems and incorrect healing)
            // also discard the startingMer as this will always be in the list somewhere
            int viableMerVariants = 0;
            for (int v = 0; v < noMerVariants; v++)
            {
                ulong variantMer = merVariants[v].variantMer;
                Depth variantDepth = merVariantDepths[v];
                int depth = variantDepth.depth;
                bool lowRepVariant = false;
                if (deepestVariantDepth > 0)
                    lowRepVariant = ((depth * 1000) / deepestVariantDepth) <= 100; // --> 10%
                if (lowRepVariant && depth >= readProps.OKDepth && !variantDepth.unbalanced)
                    lowRepVariant = false;
                if (variantMer != startingMer && !MerIsBad(depth, readProps.minDepth, readProps.OKDepth, -1, 0, variantDepth.unbalanced, previousMerDepth, previousPairDepth, m) && !lowRepVariant)
                {
                    merVariants[viableMerVariants] = merVariants[v];
                    merVariantDepths[viableMerVariants] = variantDepth;
                    viableMerVariants++;
                }
            }

            if ((deepestVariantDepth >= averageLoadedMerDepth * highDepthFactor) && errorModel == ErrorModel.mostlySubs)
                subFixesOnly = true;

            // multiple viable variants, so sort them so that the next loop can ignore duplicates (can't get duplicates if we're only varying the last base)
            if (viableMerVariants > 1 && variantsWanted != VariantWanted.varyLast)
                Array.Sort<MerVariant, Depth>(merVariants, merVariantDepths, 0, viableMerVariants, merVariantComparer);

            ulong previousMer = 0xffffffffffffffff;

            // stop looking variants once CountFollowers fails
            bool cfOK = true;

            for (int v = 0; v < viableMerVariants; v++)
            {
                ulong currentVariant = merVariants[v].variantMer;
                int lengthDelta = merVariants[v].lengthDelta;
                int variantRepairsLeft = repairsLeft;
                int variantMerCount = currentMerCount;

                // it's possible to have duplicates in this list, so just ignore them in this pass
                if (previousMer == currentVariant)
                    continue;
                previousMer = currentVariant;

                // do the next variants look likely (have suitable scores)
                Depth variantDepth = merVariantDepths[v];
                Sequence variantRead = AllocateSequence();
                currentRead.CopyTo(variantRead);

                variantMerCount = UpdateRead(variantType, variantRead, m, currentVariant, lengthDelta, variantMerCount, readProps.merCount);

                // don't continue with a variant without a pair
                int pairDepth = GetPairDepth(variantRead, currentVariant, readProps.mers, readProps.depths, m, readProps.minDepth);

                // deferred bad-pair check above (pair depth not known) so do it now - kMer is known not to be bad on just the merDepth
                if (pairDepth != -1 && MerIsBad(variantDepth.depth, readProps.minDepth, readProps.OKDepth, pairDepth, readProps.minPairDepth, variantDepth.unbalanced, previousMerDepth, previousPairDepth, m))
                {
                    ReturnSequence(variantRead);
                    continue;
                }

                // don't allow ggg-type artefacts through if there are alternatives
                if (viableMerVariants > 1 && variantDepth.unbalanced && MerIsHomopolymer(currentVariant))
                {
                    ReturnSequence(variantRead);
                    continue;
                }

                kMerProperties variantMer = AllocateMerProperty();
                //merProperties variantMer = new merProperties();
                variantMer.fixType = variantType;
                variantMer.variant = currentVariant;
                variantMer.lengthDelta = lengthDelta;
                variantMer.depth = variantDepth.depth;
                variantMer.unbalanced = variantDepth.unbalanced;
                variantMer.pairDepth = pairDepth;
                variantMer.fixContext = variantRead;

                //if (tracing == traceFollowers)
                //    TraceClumpAdd(Indent(indentDepth) +
                //                     "var " + v + " " + fixTypes[variantType] + " " + tempHealingRead + " @" + m);

                if (tracing == Tracing.traceFollowers)
                    TraceFollowers(checkReason, variantType, lengthDelta, m, indentDepth, kMers.ExpandMer(currentVariant, merSize), v, variantDepth.depth, variantRepairsLeft, variantRead.ToString(), variantMerCount, readProps);

                cfOK = CountFollowers(subFixesOnly, variantMer, m, readProps, variantMerCount, variantRepairsLeft, indentDepth, ref thmCalls);

                // if the variant mer has a viable score, record *this* repair as producing a follower
                if (variantMer.depth >= readProps.OKDepth && (variantMer.pairDepth == -1 || variantMer.pairDepth >= readProps.OKPairDepth)) 
                    variantMer.goodFollowers++;
                if (variantMer.depth >= readProps.minDepth || variantMer.pairDepth >= readProps.minPairDepth) 
                    variantMer.allFollowers++;

                variantMer.sum += variantDepth.depth;
                variantMer.maxFollowers = variantMer.merCount - m;
                variantMer.fixes++;

                variantMers.Add(variantMer);
                variantsAdded++;

                if (tracing >= Tracing.traceRead)
                {
                    string traceMsg = Indent(indentDepth) +
                                        "[" + v + "]" + fixNames[(int)variantType] + "(ld=" + variantMer.lengthDelta + ")" + "@" + m + " " + kMers.ExpandMer(currentVariant, merSize) +
                                        " sum=" + variantMer.sum + " fx " + variantMer.fixes +
                                        " fo=" + variantMer.goodFollowers + "/" + variantMer.allFollowers + "/" + variantMer.maxFollowers + "/" + variantMerCount + " nf@" + variantMer.mersToNextFix + " fc@" + variantMer.mersToFirstChoice + " cfOK " + cfOK;
                    if (tracing == Tracing.traceFollowers)
                        TraceClumpAdd(traceMsg);
                    variantMer.breadcrumbs = traceMsg + Environment.NewLine + variantMer.breadcrumbs;
                }

                if (!cfOK)
                    break;

            } // foreach viable variant

            ReturnMerVariants(merVariants);
            ReturnDepths(merVariantDepths);

            return cfOK;
        }

        //private static int RefreshFollowerCounts(List<kMerProperties> variants, Sequence read, int m, int currentMerCount, ReadProperties readProps, bool tryHarder, bool subFixesOnly, int repairsLeft, int indentDepth, ref int thmCalls)
        //{
        //    int longestFollowers = 0;

        //    foreach (kMerProperties variantMer in variants)
        //    {
        //        variantMer.allFollowers = 0;
        //        variantMer.goodFollowers = 0;
        //        variantMer.sum = 0;
        //        variantMer.fixes = 0;

        //        CountFollowers(subFixesOnly, variantMer, m, readProps, variantMer.merCount, repairsLeft, indentDepth, ref thmCalls);

        //        // if the variant mer has a viable score, record *this* repair as producing a follower
        //        if (variantMer.depth >= readProps.OKDepth && (variantMer.pairDepth == -1 || variantMer.pairDepth >= readProps.OKPairDepth)) 
        //            variantMer.goodFollowers++;
        //        if (variantMer.depth >= readProps.minDepth || variantMer.pairDepth >= readProps.minPairDepth) 
        //            variantMer.allFollowers++;

        //        variantMer.sum += variantMer.depth;
        //        if (variantMer.fixType != FixTypes.fixNone)
        //            variantMer.fixes++;

        //        // calculate highest possible follower count (inclusive of current base - count-index --> adds 1)
        //        variantMer.maxFollowers = variantMer.merCount - m;

        //        // longest follower 
        //        if (variantMer.allFollowers > longestFollowers)
        //            longestFollowers = variantMer.allFollowers;
        //    }

        //    return longestFollowers;
        //}

        private static void TraceFollowers(kMerState checkReason, FixTypes variantType, int lengthDelta, int m, int indentDepth, string merVariant, int v, 
                                           int depth, int repairsLeft, String readContext, int merCount, ReadProperties readProps)
        {
            int startOfRest = m + merSize;
            int contextLength = readContext.Length - startOfRest;
            if (contextLength > merSize)
                contextLength = merSize;
            string context = "(end)";
            if (contextLength > 0)
                context = readContext.Substring(startOfRest, contextLength);

            TraceClumpAdd(Indent(indentDepth) + "@" + m + " counting followers for " +
                            "[" + v + "] " + fixNames[(int)variantType] + "(ld=" + lengthDelta + ")" + " " + merVariant + "+" + context +
                            " with " + repairsLeft + " fixes left; sum=" + depth + "; mc=" + merCount + " " + 
                            "d=" + readProps.minDepth + "/" + readProps.OKDepth + " pd=" + readProps.minPairDepth + "/" + readProps.OKPairDepth + 
                            " " + checkReason.ToString());
        }

        // Generates a 'following' mer by taking the base following the end of a mer variant and 
        // concatenating it to the end of the variant. Returns true if such a mer could be generated.
        // Only used when counting followers. 
        private static bool GenerateNextMer(ulong packedMer, Sequence currentRead, int m, int currentMerCount, int OKDepth, int OKPairDepth,
                                            out ulong nextPackedMer, out char followingBase)
        {
            // already at end of read
            if (m == currentMerCount)                                
            {
                nextPackedMer = 0;
                followingBase = '\0';
                return false;
            }
            else
            {
                // get next base from read and shift into current kMer to generate the next kMer
                char nextBase = currentRead.Bases[m + merSize - 1];
                // and the base following the new kMer if possible
                if (m+1 < currentMerCount)
                    followingBase = currentRead.Bases[m + merSize];
                else
                    followingBase = '\0';   // avoid triggering end-of-homopolymer at end of read
                long nextBasePacked = kMers.BaseCharToInt(nextBase);
                //int nextBasePacked = "ACGT".IndexOf(nextBase);
                if (nextBasePacked < 0)
                {
                    return FindBestReplacementForNs(currentRead, m, OKDepth, out nextPackedMer);
                }
                else
                    nextPackedMer = packedMer << 2 | (ulong)nextBasePacked << (64 - merSize * 2);
                //string merIn = kMers.ExpandMer(packedMer);
                //string merOut = kMers.ExpandMer(nextPackedMer);
                return true;
            }
        }

        private static int GenerateMerSubVariants(ulong mer, MerVariant[] merVariants, int startingIdx, VariantWanted variantsWanted)
        {
            int start = 0;
            if (variantsWanted == VariantWanted.varyLast)
                start = merSize - 1;
            int variantsAdded = 0;

            ulong baseMask = 0xc000000000000000;

            for (int m = start; m < merSize; m++)
            {
                ulong merWithHole = mer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merWithHole | newBase;

                    merVariants[startingIdx + variantsAdded].variantType = FixTypes.fixSub;
                    merVariants[startingIdx + variantsAdded].variantMer = merVariant;
                    merVariants[startingIdx + variantsAdded].lengthDelta = 0;
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        private static int GenerateMerDelVariants(ulong mer, MerVariant[] merVariants, VariantWanted variantsWanted)
        {
            int start = 1;                                          // never try inserting at mer[0] as this is effectively shifting the read 1 base to the left
            if (variantsWanted == VariantWanted.varyLast)
                start = merSize - 1;
            int variantsAdded = 0;

            ulong maskMer = 0xffffffffffffffff << (64 - merSize * 2);   // just the RHS merSize bits

            for (int m = start; m < merSize; m++)
            {
                ulong maskLHS = 0xffffffffffffffff << (64 - (m * 2));   // the retained LHS part of the mer
                ulong merLHS = mer & maskLHS;
                ulong maskRHS = ~maskLHS;                               // the retained RHS part of the mer
                ulong merRHS = ((mer & maskRHS) >> 2) & maskMer;
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merLHS | newBase | merRHS;
                    if (variantsWanted == VariantWanted.varyAnyOne)                           // being careful about dels that effectively shift the read left one base in the genome
                    {
                        ulong topBase = merVariant & 0xc000000000000000;        // just the first base of the variant
                        ulong shiftedMer = ((mer >> 2) | topBase) & maskMer;    // starting mer with first base from variant inserted at start

                        if (merVariant == shiftedMer)                           // this variant didn't really change anything
                            continue;
                    }
                    merVariants[variantsAdded].variantType = FixTypes.fixDel;
                    merVariants[variantsAdded].variantMer = merVariant;
                    merVariants[variantsAdded].lengthDelta = +1;
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        private static int GenerateMerInsVariants(ulong mer, Sequence read, int fbi, int maxInsLength, MerVariant[] merVariants, VariantWanted variantsWanted)
        {
            int variantsAdded = 0;
            int start = merSize - 1;                                    // assume we're just varying the last base

            if (variantsWanted == VariantWanted.varyAnyOne)             // if we can vary any base, we want to adjust the starting point 
            {                                                           // to get past any initial homopolymer run as these will all be equivalent
                ulong maskTopBase = 0xc000000000000000;
                start = 0;
                ulong firstBase = mer & maskTopBase;
                ulong shiftedMer = mer;
                char followingBase = read.Bases[fbi];
                while (firstBase == (shiftedMer & maskTopBase) && start < merSize)
                {
                    start++;
                    shiftedMer = shiftedMer << 2;
                }

                // if it was a fully homopolymer kMer, just delete the last base
                if (start == merSize)
                    start = merSize - 1;

                ulong newBase = (ulong)kMers.BaseCharToInt(followingBase) << (64 - merSize * 2);

                for (int m = start; m < merSize; m++)
                {
                    ulong maskLHS = 0xffffffffffffffff << (64 - (m * 2));   // the retained LHS part of the mer
                    ulong merLHS = mer & maskLHS;
                    ulong maskRHS = ~maskLHS;                               // the RHS part of the mer (including the base to be excised)
                    ulong merRHS = ((mer << 2) & maskRHS);                  // and then without the right-most base
                    ulong merVariant = merLHS | newBase | merRHS;
                    if (merVariant == mer)                                  // don't include starting mer as one of its own variants
                        continue;
                    merVariants[variantsAdded].variantType = FixTypes.fixIns;
                    merVariants[variantsAdded].variantMer = merVariant;
                    merVariants[variantsAdded].lengthDelta = -1;
                    variantsAdded++;
                }
            }
            else
            {
                // vary only the last base but allow for multiple consecutive deletions
                int consecutiveIns = 0;

                for (int g = 0; g < maxInsLength; g++)
                {
                    char followingBase = read.Bases[fbi];
                    ulong newBase = (ulong)kMers.BaseCharToInt(followingBase) << (64 - merSize * 2);
                    ulong merWithoutLastBaseMask = (ulong)0xfffffffffffffff << ((64 - (merSize) * 2) + 2);
                    ulong merWithoutLastBase = mer & merWithoutLastBaseMask;
                    ulong merVariant = merWithoutLastBase | newBase;
                    consecutiveIns++;
                    merVariants[variantsAdded].variantType = FixTypes.fixIns;
                    merVariants[variantsAdded].variantMer = merVariant;
                    merVariants[variantsAdded].lengthDelta = -consecutiveIns;
                    variantsAdded++;
                    fbi++;
                }
            }
            
            return variantsAdded;
        }

        static void RateReporter()
        {
            int monitorCounter = 0;
            DateTime lastTimeAwake = DateTime.Now;
            long lastReadsCount = 0;

            while (true)
            {
                signalReporter.WaitOne(reportInterval);
                if (stopReporter)
                    break;

                monitorCounter += monitorInterval;
                if (monitorCounter >= reportInterval)
                {
                    monitorCounter = 0;
                    DateTime timeNow = DateTime.Now;
                    double timeTaken = (timeNow - lastTimeAwake).TotalSeconds;
                    lastTimeAwake = timeNow;
                    long currentReadsCount = progressReads;
                    long readsInTime = currentReadsCount - lastReadsCount;
                    lastReadsCount = currentReadsCount;
                    int readsRate = (int)((double)readsInTime / timeTaken);

                    if (programPhase == ProgramPhase.loading)
                        Console.WriteLine("loading kMers and pairs");
                    if (programPhase == ProgramPhase.healing)
                        Console.WriteLine(progressHealedReads + " healed" + ", " + progressOKReads + " OK from " + currentReadsCount + " reads at " + readsRate + " rps");
                }
            }
        }
    }

    public class HealingThreadParams
    {
        public int threadNumber;
        public BufferedReader[] bufferedReadsFiles;
        public BufferedWriter[] bufferedHealedReads;
        public BufferedWriter[] bufferedProblemReads;
        public BufferedWriter[] bufferedSingleReads;
        public BufferedWriter[] bufferedHealedQuals;
        public BufferedWriter[] bufferedProblemQuals;
        public BufferedWriter[] bufferedSingleQuals;
    }

    public class kMerProperties
    {
        public bool validVariant = true;                    // this variant can be considered (or has been deleted)
        public bool markedVariant = false;                  // this variant marked for possible deletion
        public FixTypes fixType = FixTypes.fixNone;
        public FixTypes nextFixType = FixTypes.fixNone;
        public int mersToNextFix = 0;
        public int mersToFirstChoice = 0;
        public bool perfectFix = false;
        public ulong variant = 0;
        public int lengthDelta = 0;
        public int depth = 0;
        public bool unbalanced = false;
        public int pairDepth = -1;
        public int sum = 0;
        public int goodFollowers = 0;
        public int allFollowers = 0;
        public int maxFollowers = 0;
        public int merCount = 0;
        public int fixes = 0;
        public string breadcrumbs = null;
        public Sequence fixContext = null;
    }

    public struct Depth
    {
        public int depth;
        public bool unbalanced;
    }

    public struct MerVariant
    {
        public FixTypes variantType;                    // none, sub, ins, del, ...
        public ulong variantMer;                        // kMer
        public int lengthDelta;                         // number of bases subsumed/added to read by this variant
    }

    public class MerComparerClass : IComparer <MerVariant>
    {
        public int Compare(MerVariant x, MerVariant y)
        {
            return x.variantMer.CompareTo(y.variantMer);
        }
    }

    public class ReadProperties
    {
        // the mers[] and depths[] are maintained during correction so they should always be accurate, and merCount says how many kMers there are 
        public ReadState initialReadState;              // state of read when first examined
        public ReadState finalReadState;                // state of read after correction
        public int merCount;                            // currently active kMers in the read 
        public ulong[] mers;                            // tiled kMers for the read (maintained during correction)
        public int[] depths;                            // summed kMer depths for the read (maintained during correction)
        public bool[] zeroStrand;                       // at least one strand has zero count
        public bool[] balancedDepths;                   // depths on both strands are about equal
        public int balancedMerCount;                    // how many of these 'balanced' kMers did we find?
        public int pairsCount;                          // currently active pairs in the read
        public int[] pairDepths;                        // the depths for the forward pairs starting with each of the corresponding kMers (not maintained during correction)
        public int[] merChangeCost;                     // 0/1 depending on whether the corresponding kMer was changed. Used to catch rewriting
        public int minDepth;                            // minimum depth for an acceptable mer 
        public int OKDepth;                             // target depth for a 'good' mer in this read
        public int initialOKDepth;                      // first OK depth calculated for the read - used to determine if recalculated depths have changed enough to force a re-scan
        public int minPairDepth;                        // minimum depth for an acceptable pair
        public int OKPairDepth;                         // target pair depth (good)
        public int firstGoodMer;                        // index of first kMer considered 'good' - used to limit reverse correction
        public int startOfNoisyTail;                    // start of poor quality region at the tail of the read (kMer#)
        public bool hmZeroPresent;                      // homopolymer with zero depth on one strand. Highly likely to be an artefact
        public bool deepUnbalancedPresent;              // were deep unbalanced kMers found in the read (adapters)
        public bool unbalancedRead;                     // is the read as a whole unbalanced (not necessarily deep, just one side is consistently low)
        public bool depthsRecalculated;                 // depth thresholds were changed as a result of a correction - need to retry correcting the whole read
        public bool readHasBeenChanged;                 // has the read been changed at all
        public bool healingAbandoned;                   // was the last healing pass abandoned - rewrites etc
        public int abandonedAtM;                        // kMer number where read was abandoned (for later trimming if needed)
        public AbandonReason abandonReason;             // and why was it abandoned
        public int remainingBadMers;                    // how many uncorrected 'bad' kmers remain in the 'corrected' read
        public int changedMers;                         // how many kMers were corrected

        public ReadProperties(int defaultLength)
        {
            mers = new ulong[defaultLength];
            depths = new int[defaultLength];
            pairDepths = new int[defaultLength];
            merChangeCost = new int[defaultLength];
            zeroStrand = new bool[defaultLength];
            balancedDepths = new bool[defaultLength];
            Initialise();
        }

        public void Initialise()
        {
            initialReadState = ReadState.readUnknown;
            finalReadState = ReadState.readUnknown;
            merCount = 0;
            pairsCount = 0;
            minDepth = 0;
            OKDepth = 0;
            minPairDepth = 0;
            OKPairDepth = 0;
            firstGoodMer = -1;
            startOfNoisyTail = 0;
            hmZeroPresent = false;
            deepUnbalancedPresent = false;
            unbalancedRead = false;
            depthsRecalculated = false;
            initialOKDepth = -1;
            readHasBeenChanged = false;
            healingAbandoned = false;
            abandonedAtM = -1;
            abandonReason = AbandonReason.notAbandoned;
            remainingBadMers = 0;
            changedMers = 0;
            Array.Clear(merChangeCost, 0, merChangeCost.Length);
        }

        public void Reset()
        {
            initialReadState = ReadState.readUnknown;
            finalReadState = ReadState.readUnknown;
            firstGoodMer = -1;
            readHasBeenChanged = false;
            healingAbandoned = false;
            abandonedAtM = -1;
            abandonReason = AbandonReason.notAbandoned;
            remainingBadMers = 0;
            changedMers = 0;
            Array.Clear(merChangeCost, 0, merChangeCost.Length);
        }
    }

    public class HealingStats
    {
        public long readsRead = 0;                                   // how many reads were read from the input files
        public long OKReadsWritten = 0;                              // how many OK reads were written
        public long correctedReadsWritten = 0;                       // how many reads were healed and written
        public long discardedBroken = 0;                             // how many 'broken' reads were discarded
        public long discardedOK = 0;                                 // how many good reads were dropped because their pair was faulty
        public long brokenReadsFound = 0;                            // how many reads were faulty and were still faulty after correction
        public long shortReadsFound = 0;                             // how many reads were too short to correct
        public long OKReads = 0;                                     // how many reads did not need healing (on initial look)
        public long OKReadsChecked = 0;                              // and how many reads looked like they might have a problem but proved to be OK
        public long healedReads = 0;                                 // how many reads were healed 
        public long readsNotHealedAtAll = 0;                         // how many reads needed healing but could not be healed at all 
        public long healedFirstPass = 0;                             // reads that were healed on the first pass
        public long healedAltPass = 0;                               // reads that were healed by a second try-harder (look for alts, indels) pass
        public long healedRCPass = 0;                                // reads that were healed by trying again in the reverse direction
        public long hdubReads = 0;                                   // how many reads had high-depth, unbalanced regions (adapters) trimmed?
        public long tooDeepReads = 0;                                // how many reads were skipped because they were deeper than the request maximum depth
        public long tooShortReads = 0;                               // how many reads were shorter/poorer than requested after healing/trimming
        public long abandonedReads = 0;                              // how many reads were abandoned (poor quality, rewriting, explosion or too many Ns)
        public long abandonedRewriting = 0;                          // ... decided we were just rewriting the read
        public long abandonedTree = 0;                               // ... tree got too big (cache or thm calls)
        public long abandonedNs = 0;                                 // ... read had too many Ns
        public long abandonedFutile = 0;                             // ... uncorrectable error and nothing downstream
        public long mers = 0;                                        // how many mers were in these reads
        public long replacedMers = 0;                                // how many of these mers were healed
        public long trimmed = 0;                                     // how many reads needed trimmming after correction
        public long extended = 0;                                    // how many shortened Illumina reads were restored their fixed length
        public long[] fixesByType = new long[(int)FixTypes.maxFixTypes];   // recording types of fixes used
    }

}


