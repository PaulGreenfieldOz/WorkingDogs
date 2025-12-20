using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.ExceptionServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using WorkingDogsCore;

namespace Kelpie
{
    // Kelpie is a targeted PCR-like assembler that will extract and assemble a genomic region specified by a pair of primer sequences.
    // A common use is extracting full-length marker genes from WGS metagenomic data, and the Kelpie-assembled sequences can be used
    // by conventional amplicon pipelines to build community profiles.
    //
    // A full description of Kelpie and its effectiveness can be found at https://peerj.com/articles/6174/
    //
    // Kelpie has been developed extensively since the initial release described in the PeerJ paper. The focus of this work has been
    // to improve accuracy for low abundance organisms (working further into the tail) and to support starting with unfiltered
    // sequence data files. These changes have been so extensive that this version is being called Kelpie V2. Kelpie (V1) and 
    // Kelpie V2 are fundamentally identical, but Kelpie V2 makes better use of kMer context to improve accuracy and performance, 
    // and to tighten up the initial filter construction/filtering stage to make working with unfiltered data feasible. Major parts 
    // of Kelpie V2 are now multi-threaded as well to improve overall performance.
    //
    // Kelpie runs in two phases: 
    //      - extracting and trimming the reads that cover the inter-primer region
    //      - extending the starting reads just found, using kMers derived from the covering reads
    //
    // usage: Kelpie -f forwardPrimer -r reversePrimer readsToFilterFNP extendedReadsFN (and other options)
    // 
    // Extracting/trimming 
    // -------------------
    //
    // Takes either reads filtered for a specific genomic region such as 16S (via FilterReads or an equivalent HMM) 
    // or unfiltered WGS reads (in v2) and first finds just those reads that include one of the specified primers
    //  at their start or end. These reads are then tiled for kMers to build an initial filter set. The remaining reads
    // are passed over this filter to find reads that overlap with the starting/ending reads. These additional reads are tiled and their
    // kMers are added to the filter. This process is repeated until the rate of new kMer discovery has dropped sufficiently. If pairs of files
    // were provided, the kMer sets are merged and only kMers found in both sets can be kept (loose/strict option).
    //
    // A full primer region looks like (double-stranded)
    // 5' ....F.....................(R)..... 3'  (+)  ------->
    // 3' ....(F)....................R...... 5'  (-)  <-------
    //
    // The sequencing reaction always goes in the same direction 5' --> 3' so (-) strand reads will be reversed
    // 
    // Each of these two strand-reads can be read either in a forward direction (R1) or reverse-complemented (R2)
    // 
    // (+) ....F.....................(R).....
    // R1  ....F.....................(R).....
    // R2  ....R.....................(F).....
    //
    // (-) ....R.....................(F)..... (after reversal)
    // R1  ....R.....................(F).....
    // R2  ....F.....................(R).....
    //
    // Scanning the reads and building up the filter always goes from left to right, so the starting primers are F and R,
    // and the ending primers are (F) and (R). Reads containing one of these starting primers are found, trimmed to the start of the
    // respective primer, and tiled for the initial set of filter kMers. These selected reads are added to the list of reads to be written.
    // The remaining reads are then scanned for any that start with a kMer in the starting set, and these reads are removed from the to-be-scanned list,
    // tiled and their kMers are added to the filtering set. Reads containing an ending primer are only tiled up to this primer. This process
    // continues until there are either no new matching reads (and so no additions to the filter), or the rate of new reads looks like it is indicating
    // that we've reached the end and may be about to shoot past one or more unrecognised ending primers.
    // 
    // Once the filter set is built, all the reads are passed by it again and matching reads are added to the to-be-written list. These
    // reads are then written. Those reads that start with a forward primer are marked with FP, and those that end with (F) are marked with 
    // FP'. The FP' reads will be reverse-complemented before being extended, ensuring all such 'starting' reads start with the forward primer 
    // (and can all be extended left-to-right). These starting reads are marked so they can be extended in the next phase.
    //
    // Extending
    // ---------
    //
    // The inter-primer reads found in the previous stage are then turned into a set of 32-mers, and this is then denoised to remove error-like low coverage 
    // kMers. The hashed longer kMers (contexts) are then generated from the reads, after their 32-mer ends are checked against the denoised 32-mer table. 
    // Such kMer contexts are generated for each length from 32 to the read length, in steps of 4 bases. The kMer pairs are saved as hashes 
    // of the 32-mers found at the start and end. This allows the use of fast binary hash tables.
    //
    // Each of the starting reads is checked, trimmed and then extended, one base at a time. At each iteration of the extension loop, the ending 32-mer for 
    // each possible extension is generated and then checked for viability against the 32-mer hash table. If more than one extension base is viable, all of
    // the longer kMers are generated and checked. If only a single extension is viable at any point, this base is chosen and the extension loop moves on.
    // The extension loop terminates when either the extended read ends with a terminal primer, or when none of the possible extensions are viable.
    //
    // If there are multiple viable extensions, each of the viable extended reads is recursively explored to see of any of them end up finishing with a terminal
    // primer. If only one of the extensions can reach the end of the inter-primer region, it is chosen and the extension of the reads finishes. 
    // If more than one extension can reach a terminal primer, one of them is chosen at random, in proportion to the repetition 
    // depth of the 32-mers at the end of the extended read. If none of the extension reach a terminal primer, the extension that results in the longest extended 
    // read is chosen, although this will often result in an incomplete extension that will be discarded anyway.
    //
    // Finally the set of extended reads are trimmed of their (partial) primers and written, and any reads that didn't reach a terminal primer are written to a 'discards' file. 
    // 
    //
    // -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT S1-2_270_reads_anonymous_R?_16S_20_fz_25_qt.fa CAMI_medium_16S_extended_1.2.fa -stats -log
    // -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA PlatesGH_all_arthmito_R?_tq30m100.fa All_arthmito_pairsTQ_amplicons.fa -min 310 -stats
    // -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT Fav95_?_16S_20_fz_25_qt.fa Fav_16S_extended_1.2.fa -stats -log
    // -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT c:\SeqTemp\MECSM\Fav95_?.fastq Fav_UF_extended_2.0.fa -unfiltered -tmp c:\SeqTemp\K2 -mm 1 -primers
    // -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT C:\SeqTemp\CAMI_medium\S1-2_270_reads_anonymous_R?.fq CAMI_medium_16S_extended_2.0S1.fa -unfiltered -tmp c:\SeqTemp\K2 -mm 1 -strict
    // -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA Sample_PRO1747_PlateG_*_R?_val_?_COI.fa PlatesGH_2003_BF3BR2_strictFF.fa -filtered -strict -tmp C:\SeqTemp\K2 -mm 0
    // -f TCNACNAAYCAYAARRAYATYGG -r TANACYTCNGGRTGNCCRAARAAYCA PlatesGH_2003\Sample_PRO1747_PlateG_A* _R? _val_?.fq PlatesGH_2003_FolmD_A_strict.fa -unfiltered -strict -tmp C:\SeqTemp\K2 -mm 0
    // -f CTAAATATAAATCTAGAGTTAATAA -r TTAATAATGCAAGTTTATCAAACTC F:\AmoebaTemp\IRE_CC45CANXX_TGACCACT_L00?_R?.fastq IRE_node4_plus.fa -mm 2 -min 1950 -max 2000 -unfiltered -strict -primers -matches -tmp c:\SeqTemp\K2 -kept I:\KelpieTests\KeptIRE -save IRE_node4_plus -nolcf
    // -t 16 -filtered -f AGAGTTTGATCMTGGCTCAG,TTCYGKTTGATCCYGSCRGA -r TACRGYTACCTTGTTACGACTT -mm 2 -length 1500 -tmp c:\KelpieTemp SeqData\JASP12LO_*_16S.fastq JASP12LO_16SABF_16SF.fa -log

    class Program
    {
        const string version = "V2.3.4";

        const bool cacheThings = true;

        enum PairedReadsWanted
        {
            unspecified,
            paired,
            unpaired
        }

        enum VariantsWanted
        {
            allSingle,
            allDouble,
            lastSingle
        }

        enum ScanOpts
        {
            atStart,
            atEnd,
            asIs,
            RC
        }

        const int FPpt = 0;
        const int RPpt = 1;
        const int rcFPpt = 2;
        const int rcRPpt = 3;

        const int fwd = 0;
        const int rvs = 1;
        const int asRead = 0;
        const int rced = 1;
        static readonly int[] directions = { fwd, rvs };
        const uint startingFwd = 0x00000000;
        const uint startingRvs = 0x40000000;
        const uint endingFwd = 0x80000000;
        const uint endingRvs = 0xc0000000;

        const int kMerSize = 32;                // size of kMers used in filters
        const int shortestContextSize = 40;     // shortest long kMer/kMer context length
        const int contextStride = 4;            // stride for (extension) context lengths
        const int maxRecursion = 10;            // limit depth of recursion during downstream tree exploration
        const int shortestCoreLength = 15;      // minimum size of a primer core (only matters for short primers)
        const int degenerateHCL = 2;            // rhs part of primer that is not allowed to have variants - primer has to latch somewhere - used for degenerate primers only
        const int errorRate = 100;              // estimated error rate - 1/N - default is 1/100 == 1%. 
        const int pairedReadKML = 100;          // length of kMer used for paired-read coverage test

        const int shortestContextLength = 44;        // shortest allowed context length (and min length allowed for final starting reads)
        const int shortestOtherContextLength = 44;   // shortest allowed context from other reads
        const int filterContextStride = 1;      // filter context lengths increase by this length

        const int readsPerPartition = 5000000;  // WGS reads in each partition
        const int readsInBatch = 1000;          // how many reads are fetched in each ReadReads call

        const int statsNextIdx = 0;             // how many times did we look for the next kMer
        const int statsSingleIdx = 1;           // only a single next kMer
        const int statsDownstreamIdx = 2;       // needed to look downstream to resolve ambiguity
        const int statsSingleDSTIdx = 3;        // only a single downstream path led to a terminal primer
        const int statsRolledIdx = 4;           // chose kMer in proportion to kMer counts
        const int statsLongestIdx = 5;          // chose kMer that gave the longest downstream extension
        const int statsCachedResultIdx = 6;     // starting read already extended
        const int statsCachedBranchIdx = 7;     // cached branch encountered in tree exploration
        const int statsCostIdx = 8;             // alternative chosen on lowest cost
        const int statsPairedReadsRIdx = 9;     // fork resolved using paired reads (<--)
        const int statsPairedReadsFIdx = 10;    // fork resolved using paired reads (-->)
        const int statsPairedReadsFullIdx = 11; // fork resolved using paired reads (<-->)
        const int statsReadsCoveredIdx = 12;    // fork fully covered by reads
        const int statsContextsIdx = 13;        // first context-length index in stats (for contexts used to decide between alternatives)

        static char[] bases = new char[] { 'A', 'C', 'G', 'T' };

        const int head = 0;
        const int core = 1;
        static readonly int[] primerParts = { head, core };

        static bool chooseClosest = false;
        static int traceOnReadNo = -1;

        static Dictionary<ulong, List<string>> startsOfReads = new Dictionary<ulong, List<string>>(20000);

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine("usage: Kelpie [-h] [-t #threads] [-filtered|-unfiltered] -f forwardPrimer(s) -r reversePrimer(s) [-length nn[-mm]] WGSReadsFNP extendedReadsFN (" + version + ")");
                return;
            }

            // options (and their defaults)
            List<string> forwardPrimers = new List<string>();
            List<string> reversePrimers = new List<string>();
            string filteredReadsTag = null;
            int maxAllowedReadLength = 0;
            bool writeFullHelp = false;
            bool saveLog = false;
            int minAmpliconLength = 0;
            int maxAmpliconLength = 0;
            int minExtendedLength = 0;
            int maxExtendedLength = 0;              // obsolete - kept for script compatibility
            int minDepth = 2;
            bool dropLowComplexity = true;
            string extendedReadsFN = null;
            bool generateStats = false;
            bool savePrimers = false;
            bool strict = true;
            int mismatchesFP = 1;
            int mismatchesRP = 1;
            bool prefilteredReads = true;
            string tempDir = null;
            int minQual = 30;
            string keptDir = null;
            bool leaveTempFilesBehind = false;
            bool saveMatches = false;
            PairedReadsWanted pairedReadsArg = PairedReadsWanted.unspecified;
            bool pairedReads = false;
            int threads = 0;
            ParallelOptions threadsWanted = new ParallelOptions();
            bool derepAmplicons = false;

            List<string> FNParams = new List<string>();     // the set of file names or patterns to filter

            string kelpieCommand = "Kelpie (" + version + ") ";
            for (int p = 0; p < args.Length; p++)
                kelpieCommand += args[p] + " ";

            Console.WriteLine(kelpieCommand);

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    string lcArg = args[p].ToLower();

                    if (lcArg == "-h" || lcArg == "-help")
                    {
                        writeFullHelp = true;
                        continue;
                    }

                    // gather & write choice stats 
                    if (lcArg == "-stats")
                    {
                        generateStats = true;
                        continue;
                    }

                    // save actual primers found in the reads
                    if (lcArg == "-primers")
                    {
                        savePrimers = true;
                        continue;
                    }

                    // save the reads that are matched as the primer-region filter is being constructed (debug use)
                    if (lcArg == "-matches")
                    {
                        saveMatches = true;
                        continue;
                    }

                    // save extended seq files derepicated (rather than expanding each replicate set)
                    if (lcArg == "-derep")
                    {
                        derepAmplicons = true;
                        continue;
                    }

                    // force ExtendRead to choose closest depth choice (experimental)
                    if (lcArg == "-closest")
                    {
                        chooseClosest = true;
                        continue;
                    }

                    // remove kMers from filter that aren't found in both the R1 and R2 reads
                    if (lcArg == "-strict")
                    {
                        strict = true;
                        continue;
                    }

                    // keep kMers that only appear in R1 or R2
                    if (lcArg == "-loose")
                    {
                        strict = false;
                        continue;
                    }

                    // reads have been pre-filtered, and this small set will be kept in memory rather than being re-read
                    if (lcArg == "-filtered")
                    {
                        prefilteredReads = true;
                        continue;
                    }

                    // dealing with raw WGS reads, The WGS files are re-read at every iteration of the filter construction process
                    if (lcArg == "-unfiltered")
                    {
                        prefilteredReads = false;
                        continue;
                    }

                    // default is to not use low complexity kMers when building the filter set
                    if (lcArg == "-nolcf")
                    {
                        dropLowComplexity = false;
                        continue;
                    }

                    // forward primer
                    if (lcArg == "-f" || lcArg == "-forward")
                    {
                        if (!CheckForParamValue(p, args.Length, "forward primer expected after -forward"))
                            return;
                        string forwardPrimerArg = args[p + 1].ToUpper();
                        if (forwardPrimerArg.Contains(','))
                        {
                            string[] splitPrimerArg = forwardPrimerArg.Split(',');
                            foreach (string fp in splitPrimerArg)
                                forwardPrimers.Add(fp);
                        }
                        else
                            forwardPrimers.Add(forwardPrimerArg);
                        p++;
                        continue;
                    }

                    // reverse primer
                    if (lcArg == "-r" || lcArg == "-reverse")
                    {
                        if (!CheckForParamValue(p, args.Length, "reverse primer expected after -reverse"))
                            return;
                        string reversePrimerArg = args[p + 1].ToUpper();
                        if (reversePrimerArg.Contains(','))
                        {
                            string[] splitPrimerArg = reversePrimerArg.Split(',');
                            foreach (string rp in splitPrimerArg)
                                reversePrimers.Add(rp);
                        }
                        else
                            reversePrimers.Add(reversePrimerArg);
                        p++;
                        continue;
                    }

                    // save the set of reads matched by the primer-region filter
                    if (lcArg == "-save")
                    {
                        if (!CheckForParamValue(p, args.Length, "primer name expected after -save"))
                            return;
                        filteredReadsTag = args[p + 1];
                        p++;
                        continue;
                    }

                    // treat pairs of reads files as being paired (R1 & R2)
                    if (lcArg == "-paired")
                    {
                        pairedReadsArg = PairedReadsWanted.paired;
                        continue;
                    }

                    // reads files are NOT paired
                    if (lcArg == "-unpaired")
                    {
                        pairedReadsArg = PairedReadsWanted.unpaired;
                        continue;
                    }

                    // option added for EBI incorrectly merged pairs...
                    if (lcArg == "-maxreadlength")
                    {
                        if (!CheckForParamValue(p, args.Length, "max read length expected after -maxreadlength"))
                            return;
                        try
                        {
                            maxAllowedReadLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -maxreadlength parameter: " + args[p + 1]);
                            return;
                        }
                        if (maxAllowedReadLength <= 0)
                        {
                            Console.WriteLine("-maxreadlength parameter must be > 0: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    // trim read tails back until this qual level is reached (sliding window)
                    if (lcArg == "-qt" || lcArg == "-qualtrim" || lcArg == "-tq" || lcArg == "-trimqual")
                    {
                        if (!CheckForParamValue(p, args.Length, "min qual value expected after -qualTrim"))
                            return;
                        try
                        {
                            minQual = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -qualTrim parameter: " + args[p + 1]);
                            return;
                        }
                        if (minQual < 0)
                        {
                            Console.WriteLine("-qualTrim parameter must be >= 0: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    // write debugging log from extension phase
                    if (lcArg == "-log")
                    {
                        saveLog = true;
                        continue;
                    }

                    // save any extended reads that are at least this long, even if a TP wasn't found
                    if (lcArg == "-min" || lcArg == "-minlength")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a minimum length after -min");
                            return;
                        }
                        try
                        {
                            minExtendedLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a minimum length after -min: " + args[p + 1]);
                            return;
                        }
                        if (minExtendedLength <= 0)
                        {
                            Console.WriteLine("-min parameter must be > 0: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    // set maximum number of cores to use when multi-threading
                    if (lcArg == "-t" || lcArg == "-threads")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected #threads after -threads");
                            return;
                        }

                        string threadsParam = args[p + 1];
                        if (threadsParam == "max" || threadsParam == "all")
                        {
                            threads = Environment.ProcessorCount;
                        }
                        else
                        {
                            try
                            {
                                threads = Convert.ToInt32(args[p + 1]);
                            }
                            catch
                            {
                                Console.WriteLine("expected #threads after -threads: " + args[p + 1]);
                                return;
                            }
                            if (threads <= 0)
                            {
                                Console.WriteLine("-threads parameter must be > 0: " + args[p + 1]);
                                return;
                            }
                        }
                        p++;
                        continue;
                    }

                    // trace extension in log for this read
                    if (lcArg == "-trn")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected readNo after -trn");
                            return;
                        }

                        string trnParam = args[p + 1];
                        try
                        {
                            traceOnReadNo = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected readNo after -trn: " + args[p + 1]);
                            return;
                        }
                        if (traceOnReadNo <= 0)
                        {
                            Console.WriteLine("-trn parameter must be > 0: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }


                    // legacy 'max' option - replaced with 'length'
                    if (lcArg == "-max" || lcArg == "-maxlength")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a maximum length after -max");
                            return;
                        }
                        try
                        {
                            maxExtendedLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a maximum length after -max: " + args[p + 1]);
                            return;
                        }
                        if (maxExtendedLength < shortestContextSize)
                        {
                            Console.WriteLine("max length must be >= " + shortestContextSize);
                            return;
                        }
                        p++;
                        continue;
                    }

                    // expected amplicon/inter-primer region length. Used to control filtering extension iterations, and stop amplicon extension getting out of control
                    if (lcArg == "-length")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a length/range after -length");
                            return;
                        }
                        string lengthArg = args[p + 1];
                        string lowerLength = null;
                        string upperLength = null;
                        int dashIdx = lengthArg.IndexOf('-');
                        if (dashIdx == -1)
                            lowerLength = lengthArg;
                        else
                        {
                            lowerLength = lengthArg.Substring(0, dashIdx);
                            upperLength = lengthArg.Substring(dashIdx + 1);
                        }
                        int lowerLengthInt = 0;
                        int upperLengthInt = 0;
                        try
                        {
                            lowerLengthInt = Convert.ToInt32(lowerLength);
                        }
                        catch
                        {
                            Console.WriteLine("expected a length after -length: " + args[p + 1]);
                            return;
                        }
                        if (upperLength != null)
                            try
                            {
                                upperLengthInt = Convert.ToInt32(upperLength);
                            }
                            catch
                            {
                                Console.WriteLine("expected a length range after -length: " + args[p + 1]);
                                return;
                            }
                        if (lowerLengthInt < shortestContextSize)
                        {
                            Console.WriteLine("length must be >= " + shortestContextSize);
                            return;
                        }

                        minAmpliconLength = upperLengthInt > 0 ? lowerLengthInt : lowerLengthInt - (lowerLengthInt / 10);
                        maxAmpliconLength = upperLengthInt > 0 ? upperLengthInt : lowerLengthInt + (lowerLengthInt / 10);
                        p++;
                        continue;
                    }


                    // default minDepth for kMers/reads. Only useful when trying to assemble very high depth regions
                    if (lcArg == "-md" || lcArg == "-mindepth")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a minimum depth after -md");
                            return;
                        }
                        try
                        {
                            minDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a minimum depth after -md: " + args[p + 1]);
                            return;
                        }
                        if (minDepth < 1)
                        {
                            Console.WriteLine("min depth must be >= 1");
                            return;
                        }
                        p++;
                        continue;
                    }

                    // how many mismatches are allowed when matching primers (be careful as it's easy to get off-target matches)
                    // added to better handle COI primers (and used for CAMI tests as well)
                    if (args[p] == "-mm" || args[p] == "-mismatches")
                    {
                        if (!CheckForParamValue(p, args.Length, "no. of bases expected after -mismatches"))
                            return;

                        string mmParam = args[p + 1];
                        string mmpFP = mmParam;
                        string mmpRP = null;

                        int plusIdx = Math.Max(mmParam.IndexOf('+'), mmParam.IndexOf(','));
                        if (plusIdx > 0)
                        {
                            mmpFP = mmParam.Substring(0, plusIdx);
                            mmpRP = mmParam.Substring(plusIdx + 1);
                        }

                        try
                        {
                            mismatchesFP = Convert.ToInt32(mmpFP);
                            if (mmpRP != null)
                                mismatchesRP = Convert.ToInt32(mmpRP);
                            else
                                mismatchesRP = mismatchesFP;
                        }
                        catch
                        {
                            Console.WriteLine("Expected number(s) for the -mismatches parameter: " + args[p + 1]);
                            return;
                        }
                        if (mismatchesFP < 0 || mismatchesRP < 0)
                        {
                            Console.WriteLine("-mismatches parameters must be >= 0: " + mmParam);
                            return;
                        }
                        p++;
                        continue;
                    }

                    // where to store the temp files used when processing WGS files
                    if (args[p] == "-tmp" || args[p] == "-temp")
                    {
                        if (!CheckForParamValue(p, args.Length, "temp directory name expected after -tmp"))
                            return;
                        tempDir = args[p + 1];
                        p++;
                        continue;
                    }

                    // where to find already pre-processed WGS (temp) files
                    if (args[p] == "-kept")
                    {
                        if (!CheckForParamValue(p, args.Length, "kept directory name expected after -kept"))
                            return;
                        keptDir = args[p + 1];
                        p++;
                        continue;
                    }

                    Console.WriteLine("unrecognised option: " + args[p]);
                    return;

                } // options

                // not an option - must be a file name
                FNParams.Add(args[p]);
            }

            if (writeFullHelp)
            {
                WriteFullUsage();
                return;
            }

            if (forwardPrimers == null)
            {
                Console.WriteLine("no forward primer found");
                return;
            }

            if (reversePrimers == null)
            {
                Console.WriteLine("no reverse primer found");
                return;
            }
            int forwardPrimerLength = 0;
            foreach (string fp in forwardPrimers)
            {
                if (fp.Length > forwardPrimerLength)
                    forwardPrimerLength = fp.Length;
                if (fp.Length > 32)
                {
                    Console.WriteLine("forward primers have to be <= 32bp - " + fp);
                    return;
                }
            }
            foreach (string fp in forwardPrimers)
                if (fp.Length != forwardPrimerLength)
                {
                    Console.WriteLine("all forward primers have to be the same length (" + forwardPrimerLength + ") - " + fp);
                    return;
                }

            int reversePrimerLength = 0;
            foreach (string rp in reversePrimers)
            {
                if (rp.Length > reversePrimerLength)
                    reversePrimerLength = rp.Length;
                if (rp.Length > 32)
                {
                    Console.WriteLine("reverse primers have to be <= 32bp");
                    return;
                }
            }
            foreach (string rp in reversePrimers)
                if (rp.Length != reversePrimerLength)
                {
                    Console.WriteLine("all reverse primers have to be the same length (" + reversePrimerLength + ") - " + rp);
                    return;
                }

            int shortestAllowed = shortestContextSize + forwardPrimerLength + reversePrimerLength;
            if (minExtendedLength > 0 && minExtendedLength < shortestAllowed)
            {
                Console.WriteLine("min extended length must be >= " + shortestAllowed);
                return;
            }

            if (minAmpliconLength > 0 && minAmpliconLength < shortestAllowed)
            {
                Console.WriteLine("min amplicon length must be >= " + shortestAllowed);
                return;
            }

            if (minAmpliconLength > 0 && maxAmpliconLength > 0 && minAmpliconLength >= maxAmpliconLength)
            {
                Console.WriteLine("min amplicon length must be > max amplicon length");
                return;
            }

            if (minExtendedLength > 0 && maxAmpliconLength > 0 && minExtendedLength > maxAmpliconLength)
            {
                Console.WriteLine("min extended length must be <= max amplicon length");
                return;
            }

            if (minExtendedLength > 0 && maxExtendedLength > 0 && minExtendedLength > maxExtendedLength)
            {
                Console.WriteLine("min extended length must be < max extended length");
                return;
            }

            if (minExtendedLength > 0 && minAmpliconLength == 0)
                minAmpliconLength = minExtendedLength;

            if (maxExtendedLength > 0 && maxAmpliconLength == 0)
                maxAmpliconLength = maxExtendedLength;

            if (maxExtendedLength > maxAmpliconLength)
                maxAmpliconLength = maxExtendedLength;

            if (threads == 0)
                threads = Environment.ProcessorCount / 2;
            threadsWanted.MaxDegreeOfParallelism = threads;

            if (FNParams.Count == 0)
            {
                Console.WriteLine("no file names found");
                return;
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("expected both readsFNP & extendedReadsFN. Got only 1 file name");
                return;
            }

            // expand file name patterns to get a list of files to be processed
            List<string> FNs = new List<string>();
            extendedReadsFN = FNParams[FNParams.Count - 1];
            string extendedReadsPrefix = extendedReadsFN;
            int lastERDotIdx = extendedReadsFN.LastIndexOf('.');
            if (lastERDotIdx > 0)
                extendedReadsPrefix = extendedReadsFN.Substring(0, lastERDotIdx);

            for (int f = 0; f < FNParams.Count - 1; f++)
            {
                string FNP = FNParams[f];
                string fnDir;
                string fnFN;
                SeqFiles.SplitFileName(FNP, out fnDir, out fnFN);
                string[] FNsFromFNP = Directory.GetFiles(fnDir, fnFN);
                if (FNsFromFNP.Length == 0)
                {
                    Console.WriteLine("file name pattern '" + FNP + "' did not match any files");
                    return;
                }
                foreach (string FN in FNsFromFNP)
                    FNs.Add(FN);
            }

            if (pairedReadsArg == PairedReadsWanted.unpaired)
                pairedReads = false;
            if (pairedReadsArg == PairedReadsWanted.paired)
            {
                if (FNs.Count % 2 == 0)
                    pairedReads = true;
                else
                {
                    Console.WriteLine("-paired option set but found uneven number of reads files (" + FNs.Count + ")");
                    return;
                }
            }
            // even number of files - assume R1/R2 pairs - filter kMers must be found in both files
            if (FNs.Count % 2 == 0 && pairedReadsArg == PairedReadsWanted.unspecified)
                pairedReads = true;

            string fnSeparator = Path.DirectorySeparatorChar.ToString();  // \ for Windows; / for Unix/Linux

            StreamWriter log = null;
            if (saveLog)
            {
                log = new StreamWriter("KelpieLog.txt");
                log.WriteLine(kelpieCommand);
            }

            // set the temp dir - and test if we've been given a valid path
            if (tempDir == null)
                tempDir = Directory.GetCurrentDirectory() + fnSeparator;
            else
            {
                try
                {
                    // add a trailing \ if the temp directory name doesn't already have one
                    if (!tempDir.EndsWith(fnSeparator))
                        tempDir += fnSeparator;
                    if (!Directory.Exists(tempDir))
                        Directory.CreateDirectory(tempDir);
                    string testTempFN = tempDir + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                    StreamWriter testTemp = new StreamWriter(testTempFN);
                    testTemp.Close();
                    File.Delete(testTempFN);
                }
                catch
                {
                    Console.WriteLine("Temp directory: " + tempDir + " was invalid");
                    return;
                }
            }

            Stopwatch totalElapsedTime = new Stopwatch();
            totalElapsedTime.Start();
            Stopwatch phase1ElapsedTime = new Stopwatch();
            phase1ElapsedTime.Start();

            // Generate the sets of allowable primer variants (degenerate bases and mismatches)
            // --------------------------------------------------------------------------------
            // Keep both forward and reverse complement forms of all primer variants
            // primers are held as head+core to better match PCR primer matching and to avoid exponential explosion when 
            // generating mismatch variants of long degenerate primers (such as COI)

            // Primers are split into two parts - head & core - to better emulate PCR (and primer binding) and so to allow more mismatches
            // The 3' end (RHS) binds more tightly and can tolerate the fewest differences. The 5' end binds less tightly and can tolerate more mismatches.
            // The 3' end is always towards the centre of the amplicon, given reversals and strand directions. Any 'hard' matching bases will be at this 3' end.
            // F hhhhcccccccccX, F' Xcccccccccchhhh, R hhhhcccccccccX, R' Xcccccccccchhhh

            char[] degenerateBases;
            Dictionary<char, List<char>> degenerateBaseExpansions;
            kMers.InitialiseDegenerateBaseTables(out degenerateBases, out degenerateBaseExpansions);

            // starting primers are either the Forward primer or the Reverse primer 
            // filter extension stops when we hit a termination primer - rc(Reverse) or rc(Forward)
            HashSet<string>[] forwardPrimersFwd = new HashSet<string>[2];
            HashSet<string>[] forwardPrimersRC = new HashSet<string>[2];
            HashSet<string>[] reversePrimersFwd = new HashSet<string>[2];
            HashSet<string>[] reversePrimersRC = new HashSet<string>[2];
            HashSet<string>[] terminatingPrimers = new HashSet<string>[2]; // derived from reversePrimersRC
            foreach (int p in primerParts)
            {
                forwardPrimersFwd[p] = new HashSet<string>();
                forwardPrimersRC[p] = new HashSet<string>();
                reversePrimersFwd[p] = new HashSet<string>();
                reversePrimersRC[p] = new HashSet<string>();
                terminatingPrimers[p] = new HashSet<string>();
            }

            // primer core length - balance between being long enough for distinctiveness and short enough to allow some bases for the 'head' 
            int primerLength = Math.Min(forwardPrimerLength, reversePrimerLength);
            int primerCoreLength = primerLength * 3 / 4;
            if (primerCoreLength < shortestCoreLength)
                primerCoreLength = shortestCoreLength;
            if (primerLength < primerCoreLength)
                primerCoreLength = primerLength;

            int fwdPrimerHeadLength = forwardPrimerLength - primerCoreLength;
            int rvsPrimerHeadLength = reversePrimerLength - primerCoreLength;


            // we trim away the start of the starting reads so they are compatible with the salvaged starting reads
            int startingReadTrim = forwardPrimerLength / 2;

            // remember what primer forms can be found
            List<string>[] primers = new List<string>[4];
            primers[FPpt] = new List<string>();
            primers[rcFPpt] = new List<string>();
            primers[RPpt] = new List<string>();
            primers[rcRPpt] = new List<string>();

            // forward primers
            // ---------------
            bool degeneratePrimers = PrimersQuiteDegenerate(forwardPrimers);
            // demand that the few most critical primers bases (RHS) must match exactly (no mismatches)
            int hardCoreLength = degeneratePrimers ? degenerateHCL : 0;

            foreach (string forwardPrimer in forwardPrimers)
            {
                // build variant tables for forward primer
                string forwardPrimerHead = forwardPrimer.Substring(0, fwdPrimerHeadLength);
                string forwardPrimerCore = forwardPrimer.Substring(fwdPrimerHeadLength);
                int fphAllowedMismatches = 0; ;
                int coreAllowedMismatches = 0;
                if (mismatchesFP == 0)
                {
                    fphAllowedMismatches = 0;
                    coreAllowedMismatches = 0;
                }
                if (mismatchesFP == 1)
                {
                    // allow mismatches anywhere but primer seqs with >1 mismatches will be rejected later
                    fphAllowedMismatches = 1;
                    coreAllowedMismatches = 1;
                }
                if (mismatchesFP > 1)
                {
                    // multiple mismatches - allow one(+) in head and specified number in core. Matching 'primer' seqs with too many mismatches will be rejected later
                    fphAllowedMismatches = mismatchesFP * fwdPrimerHeadLength / forwardPrimer.Length;
                    if (fphAllowedMismatches == 0)
                        fphAllowedMismatches = 1;
                    coreAllowedMismatches = mismatchesFP;
                }

                // forward primer
                // F Xcccccccccchhhh
                kMers.GenerateSeqVariants(forwardPrimerHead, 0, forwardPrimerHead.Length, degenerateBases, degenerateBaseExpansions, forwardPrimersFwd[head], fphAllowedMismatches);
                kMers.GenerateSeqVariants(forwardPrimerCore, hardCoreLength, (forwardPrimerCore.Length - hardCoreLength), degenerateBases, degenerateBaseExpansions, forwardPrimersFwd[core], coreAllowedMismatches);

                primers[FPpt].Add(forwardPrimer);

                // forward primer RC
                // F' hhhhcccccccccX
                string forwardPrimerRC = kMers.ReverseComplement(forwardPrimer);
                string forwardPrimerRCHead = forwardPrimerRC.Substring(forwardPrimerCore.Length);
                string forwardPrimerRCCore = forwardPrimerRC.Substring(0, forwardPrimerCore.Length);
                kMers.GenerateSeqVariants(forwardPrimerRCHead, 0, forwardPrimerRCHead.Length, degenerateBases, degenerateBaseExpansions, forwardPrimersRC[head], fphAllowedMismatches);
                kMers.GenerateSeqVariants(forwardPrimerRCCore, 0, (forwardPrimerRCCore.Length - hardCoreLength), degenerateBases, degenerateBaseExpansions, forwardPrimersRC[core], coreAllowedMismatches);

                primers[rcFPpt].Add(forwardPrimerRC);

            }
            Console.WriteLine("built kMer sets from forward primers");

            // reverse primers
            // ---------------
            degeneratePrimers = PrimersQuiteDegenerate(reversePrimers);
            // demand that the few most critical primers bases (RHS) must match exactly (no mismatches)
            hardCoreLength = degeneratePrimers ? degenerateHCL : 0;
            foreach (string reversePrimer in reversePrimers)
            {
                // build variant tables for reverse primer
                string reversePrimerHead = reversePrimer.Substring(0, rvsPrimerHeadLength);
                string reversePrimerCore = reversePrimer.Substring(rvsPrimerHeadLength);
                int rphAllowedMismatches = 0;
                int coreAllowedMismatches = 0;

                if (mismatchesRP == 0)
                {
                    rphAllowedMismatches = 0;
                    coreAllowedMismatches = 0;
                }
                if (mismatchesRP == 1)
                {
                    // allow mismatches anywhere but primer seqs with >1 mismatches will be rejected later
                    rphAllowedMismatches = 1;
                    coreAllowedMismatches = 1;
                }
                if (mismatchesRP > 1)
                {
                    // multiple mismatches - allow one(+) in head and specified number in core. Matching 'primer' seqs with too many mismatches will be rejected later
                    rphAllowedMismatches = mismatchesRP * rvsPrimerHeadLength / reversePrimer.Length;
                    if (rphAllowedMismatches == 0)
                        rphAllowedMismatches = 1;
                    coreAllowedMismatches = mismatchesRP;
                }

                // reverse primer
                // R Xcccccccccchhhh
                kMers.GenerateSeqVariants(reversePrimerHead, 0, reversePrimerHead.Length, degenerateBases, degenerateBaseExpansions, reversePrimersFwd[head], rphAllowedMismatches);
                kMers.GenerateSeqVariants(reversePrimerCore, hardCoreLength, (reversePrimerCore.Length - hardCoreLength), degenerateBases, degenerateBaseExpansions, reversePrimersFwd[core], coreAllowedMismatches);

                primers[RPpt].Add(reversePrimer);

                // reverse primer RC
                // R' hhhhcccccccccX
                string reversePrimerRC = kMers.ReverseComplement(reversePrimer);
                string reversePrimerRCHead = reversePrimerRC.Substring(reversePrimerCore.Length);
                string reversePrimerRCCore = reversePrimerRC.Substring(0, reversePrimerCore.Length);
                kMers.GenerateSeqVariants(reversePrimerRCHead, 0, reversePrimerRCHead.Length, degenerateBases, degenerateBaseExpansions, reversePrimersRC[head], rphAllowedMismatches);
                kMers.GenerateSeqVariants(reversePrimerRCCore, 0, (reversePrimerRCCore.Length - hardCoreLength), degenerateBases, degenerateBaseExpansions, reversePrimersRC[core], coreAllowedMismatches);

                primers[rcRPpt].Add(reversePrimerRC);

                // don't include the 'hard core' in terminating primers - only really exists to cut down on false primer matches during initial WGS filtering
                foreach (string ph in reversePrimersRC[head])
                    terminatingPrimers[head].Add(ph);
                kMers.GenerateSeqVariants(reversePrimerRCCore, 0, reversePrimerRCCore.Length, degenerateBases, degenerateBaseExpansions, terminatingPrimers[core], coreAllowedMismatches);

            }
            Console.WriteLine("built kMer sets from reverse primers");

            Dictionary<char, HashSet<char>> allowedBases = new Dictionary<char, HashSet<char>>();
            foreach (char b in new char[] { 'A', 'C', 'G', 'T' })
            {
                allowedBases.Add(b, new HashSet<char>());
                allowedBases[b].Add(b);
            }
            foreach (char b in degenerateBases)
            {
                allowedBases.Add(b, new HashSet<char>());
                foreach (char be in degenerateBaseExpansions[b])
                    allowedBases[b].Add(be);
            }

            HashSet<ulong> allPrimerCores = new HashSet<ulong>();
            foreach (string primerCore in forwardPrimersFwd[core])
                allPrimerCores.Add(kMers.CondenseMer(primerCore, primerCoreLength));
            foreach (string primerCore in forwardPrimersRC[core])
                allPrimerCores.Add(kMers.CondenseMer(primerCore, primerCoreLength));
            foreach (string primerCore in reversePrimersFwd[core])
                allPrimerCores.Add(kMers.CondenseMer(primerCore, primerCoreLength));
            foreach (string primerCore in reversePrimersRC[core])
                allPrimerCores.Add(kMers.CondenseMer(primerCore, primerCoreLength));

            long totalReadsCount = 0;
            long totalBasesTrimmed = 0;

            Dictionary<string, int>[] primersFound = new Dictionary<string, int>[4];
            for (int i = 0; i < 4; i++)
                primersFound[i] = new Dictionary<string, int>();

            // Process the pairs of files 
            // --------------------------
            //      assume that we're dealing with a single DNA sample, spread (possibly) across multiple (pairs of) files.
            //      the usual case will be a single R1/R2 pair, but some HiSeq runs produce two pairs of pairs (L1 & L2) for each sample.
            //      the R1 reads are effectively merged, as are the R2 reads.
            //
            //      the starting and ending filters are built separately for the R1 and R2 files (if 'paired' is set)
            //      if 'paired' & 'strict', reduce the filters by excluding any kMers that were seen in only the R1 or R2 set.

            // if we're handling paired reads file, loop over the files two at a time
            int fileStride = pairedReads ? 2 : 1;
            // and remember the longest read we came across
            int longestReadLength = 0;

            // Kelpie can handle either pre-filtered or raw WGS files.
            // ------------------------------------------------------
            //      unfiltered WGS reads - could be large - and with few reads of interest - so don't keep all the reads in memory
            //                             raw reads (possibly fastq) are read and saved as partitioned fasta files
            //                             subsequent passes over these reads are multi-threaded. 
            //      pre-filtered reads   - read once and kept in-memory to make it faster to do multiple passes over them)            

            // in-memory copies of (pre-filtered) reads
            List<string>[] headersToFilter = null;
            List<string>[] readsToFilter = null;

            // track which reads have already been tiled for the region filter (for each partition)
            // outer index is file within pair, list is for partitions (pre-filtered are in [0])
            // this list is pre-populated with primer reads to avoid selecting them twice
            // and to pick up the trimmed reads rather than the original ones from the WGS files
            List<HashSet<int>>[] readsAlreadyTiled = new List<HashSet<int>>[fileStride];

            // partitioned temp files of WGS reads (fasta)
            List<string>[] WGSReadsFNs = new List<string>[fileStride];
            string WGSReadsFNPrefix = extendedReadsPrefix;
            int[] partitionCount = new int[fileStride];
            // and how many reads are in each partition
            List<int>[] readsInPartitions = new List<int>[fileStride];

            // read all the WGS reads and either leave them in memory (filtered) or write them out as single-line FASTA temp files (WGS)
            // -------------------------------------------------------------------------------------------------------------------------
            if (prefilteredReads)
            {
                // saved headers+reads from all pre-filtered WGS data files
                headersToFilter = new List<string>[fileStride];
                readsToFilter = new List<string>[fileStride];

                for (int fip = 0; fip < fileStride; fip++)
                {
                    headersToFilter[fip] = new List<string>(100000);
                    readsToFilter[fip] = new List<string>(100000);
                    // in-memory pre-filtered reads regarded as partition [0]
                    readsAlreadyTiled[fip] = new List<HashSet<int>>();
                    readsInPartitions[fip] = new List<int>();
                    partitionCount[fip] = 1;
                }
            }
            else
            {
                // (temp files) saved, converted raw (unfiltered) WGS reads in FASTA format, with a single line for each sequence
                for (int fip = 0; fip < fileStride; fip++)
                {
                    WGSReadsFNs[fip] = new List<string>();
                    readsAlreadyTiled[fip] = new List<HashSet<int>>();
                    readsInPartitions[fip] = new List<int>();
                    partitionCount[fip] = 0;
                }
            }

            // we've already got a saved set of temp WGS files that we can just copy across
            if (keptDir != null)
            {
                string[][] keptFNs = new string[fileStride][];
                int keptFilesCopied = 0;

                if (!Directory.Exists(keptDir))
                    Directory.CreateDirectory(keptDir);

                Dictionary<string, int> readsInTmp = new Dictionary<string, int>(500);
                string[] metadataFNs = Directory.GetFiles(keptDir, "*_kept_metadata.txt");

                if (metadataFNs.Length == 1)
                {
                    StreamReader metadata = new StreamReader(metadataFNs[0]);
                    bool metaEOF = false;
                    while (!metaEOF)
                    {
                        string metaLine = metadata.ReadLine();
                        if (metaLine == null)
                            break;
                        int tabIdx = metaLine.IndexOf('\t');
                        string fn = metaLine.Substring(0, tabIdx);
                        string val = metaLine.Substring(tabIdx + 1);
                        if (fn == "longest")
                            longestReadLength = Convert.ToInt32(val);
                        else
                            readsInTmp.Add(fn, Convert.ToInt32(val));
                    }
                }
                else
                    Console.WriteLine("found multiple metadata files in kept directory.");

                // metadata file found and read - so trying copying the data files
                if (readsInTmp.Count > 0)
                    for (int fip = 0; fip < fileStride; fip++)
                    {
                        keptFNs[fip] = Directory.GetFiles(keptDir, "*_" + fip + ".tmp");

                        for (int f = 0; f < keptFNs[fip].Length; f++)
                        {
                            string keptFN = keptFNs[fip][f];
                            keptFilesCopied++;

                            int lastSloshIdx = keptFN.LastIndexOf(Path.DirectorySeparatorChar);
                            string lastFNPart = keptFN.Substring(lastSloshIdx + 1);
                            string tempFN = tempDir + lastFNPart;

                            File.Copy(keptFN, tempFN, true);

                            WGSReadsFNs[fip].Add(tempFN);
                            int readsInFile = readsInTmp[lastFNPart];
                            readsInPartitions[fip].Add(readsInFile);
                            partitionCount[fip]++;
                            totalReadsCount += readsInFile;
                            readsAlreadyTiled[fip].Add(new HashSet<int>());
                            Console.WriteLine("recovered " + lastFNPart);
                        }
                    }

                // told temp files had been kept, but didn't find them
                if (keptFilesCopied == 0)
                {
                    Console.WriteLine("no kept file set found in " + keptDir + ". Will leave temp files behind for copying to kept directory");
                    leaveTempFilesBehind = true;
                    keptDir = null;
                }
            }

            // normal case: read all the reads from all the files (possibly in pairs) and save/write them as temp files
            if (keptDir == null)
            {
                for (int fip = 0; fip < fileStride; fip++)
                {
                    // read the headers and seqs for every file and save them in memory or in the WGS temp file

                    StreamWriter WGSReads = null;
                    int readsInWGSPartition = 0;
                    if (!prefilteredReads)
                    {
                        WGSReads = OpenNewWGSPartition(partitionCount[fip], fip, tempDir, WGSReadsFNPrefix, WGSReadsFNs[fip]);
                        partitionCount[fip]++;
                    }

                    Sequence[] readHeaders = new Sequence[readsInBatch];
                    Sequence[] readSeqs = new Sequence[readsInBatch];
                    Sequence[] qualHeaders = new Sequence[readsInBatch];
                    Sequence[] qualScores = new Sequence[readsInBatch];
                    for (int i = 0; i < readsInBatch; i++)
                    {
                        readHeaders[i] = new Sequence(60);
                        readSeqs[i] = new Sequence(200);
                        qualHeaders[i] = new Sequence(10);
                        qualScores[i] = new Sequence(200);
                    }

                    int fileFormat = SeqFiles.formatNone;
                    bool fullQualHeader;
                    int qualOffset = SeqFiles.ResolveFastqQualAmbiguity(FNs[fip], out fullQualHeader);

                    // for each R1 or R2 file
                    for (int f = 0; f < FNs.Count; f += fileStride)
                    {
                        string FN = FNs[f + fip];
                        if (fileFormat == SeqFiles.formatNone)
                            fileFormat = SeqFiles.DetermineFileFormat(FN);

                        StreamReader reads = SeqFiles.OpenSeqStream(FN);
                        BufferedReader bufferedReader = new BufferedReader(fileFormat, reads, null);

                        Console.WriteLine("reading " + FN);

                        bool EOF = false;
                        while (!EOF)
                        {
                            int readsRead = bufferedReader.ReadReads(readsInBatch, readHeaders, readSeqs, qualHeaders, qualScores);
                            if (readsRead != readsInBatch)
                                EOF = true;

                            for (int i = 0; i < readsRead; i++)
                            {
                                Sequence header = readHeaders[i];
                                Sequence read = readSeqs[i];
                                Sequence quals = qualScores[i];
                                totalReadsCount++;

                                // used to skip any reads that contained non-ACGT bases - but CAMI data has Ns! And checks for invalid Kmers after tiling will stop these being an issue (SeqContainsAmbiguousBases(read))

                                // Poor Illumina data will have a trailing run of Gs - the default base - these should have poor qual scores but not always it would seem.
                                // Cleanest just trim them away here
                                if (read.Length >= kMerSize)
                                {
                                    ulong kMerAtEnd;
                                    Sequence.CondenseMer(read, read.Length - kMerSize / 2, kMerSize / 2, out kMerAtEnd);
                                    if (kMerAtEnd == 0xAAAAAAAA00000000)
                                    {
                                        int lastBaseIdx = read.Length - 1;
                                        while (lastBaseIdx >= 0 && read.Bases[lastBaseIdx] == 'G')
                                        {
                                            read.Length--;
                                            if (fileFormat == SeqFiles.formatFASTQ)
                                                quals.Length--;
                                            lastBaseIdx--;
                                        }
                                    }
                                }

                                if (maxAllowedReadLength > 0 && read.Length > maxAllowedReadLength)
                                    read.Length = maxAllowedReadLength;

                                //if (header.ToString() == "@A00632:395:HLMHJDSX7:4:2604:19479:18599 2:N:0:GTCTCTAC+ATAGGATA")
                                //    Debugger.Break();

                                if (fileFormat == SeqFiles.formatFASTQ && minQual > 0)
                                {
                                    totalBasesTrimmed += SeqFiles.TrimTrailingPoorQuals(read, quals, minQual, qualOffset);
                                    int poorQualsRemaining = SeqFiles.CountPoorQuals(quals, minQual, qualOffset);
                                    if (poorQualsRemaining > 0)
                                        header.Append(";pq=" + poorQualsRemaining);
                                }

                                if (prefilteredReads)
                                {
                                    headersToFilter[fip].Add(header.ToString());
                                    readsToFilter[fip].Add(read.ToString());
                                }
                                else
                                {
                                    header.Bases[0] = '>'; // force into single-line FASTA format
                                    WGSReads.WriteLine(header.Bases, 0, header.Length);
                                    WGSReads.WriteLine(read.Bases, 0, read.Length);
                                    readsInWGSPartition++;
                                    if (readsInWGSPartition == readsPerPartition)
                                    {
                                        // close off the just-filled partition
                                        WGSReads.Close();
                                        readsInPartitions[fip].Add(readsInWGSPartition);
                                        readsAlreadyTiled[fip].Add(new HashSet<int>());

                                        // and prep the next one
                                        WGSReads = OpenNewWGSPartition(partitionCount[fip], fip, tempDir, WGSReadsFNPrefix, WGSReadsFNs[fip]);
                                        partitionCount[fip]++;
                                        readsInWGSPartition = 0;
                                    }
                                }
                                if (read.Length > longestReadLength)
                                    longestReadLength = read.Length;
                            }
                        } // for all the reads in a single file
                    } // for all files (R1 or R2 only)

                    if (WGSReads == null)
                    {
                        readsInPartitions[fip].Add(readsToFilter[fip].Count);
                        readsAlreadyTiled[fip].Add(new HashSet<int>());
                    }
                    else
                    {
                        WGSReads.Close();
                        readsInPartitions[fip].Add(readsInWGSPartition);
                        readsAlreadyTiled[fip].Add(new HashSet<int>());
                    }

                } // for paired files (fip)
            } // normal (non-kept) path

            if (leaveTempFilesBehind)
            {
                StreamWriter kmd = new StreamWriter(tempDir + WGSReadsFNPrefix + "_kept_metadata.txt");
                kmd.WriteLine("longest\t" + longestReadLength);
                for (int fip = 0; fip < fileStride; fip++)
                    for (int i = 0; i < WGSReadsFNs[fip].Count; i++)
                    {
                        string keptFN = WGSReadsFNs[fip][i];
                        int lastSloshIdx = keptFN.LastIndexOf(Path.DirectorySeparatorChar);
                        string lastFNPart = keptFN.Substring(lastSloshIdx + 1);
                        kmd.WriteLine(lastFNPart + "\t" + readsInPartitions[fip][i]);
                    }
                kmd.Close();
            }

            Console.WriteLine("read " + totalReadsCount + " reads from " + FNs.Count + " files. Avg tail trim = " + (totalBasesTrimmed / totalReadsCount));

            // remember the reads that were found to start/end with one of the primers.  
            // these are the actual trimmed read and tagged header strings. Saving just the readIdx 
            // doesn't work for the WGS reads as it's not possible to update the reads within the files
            List<string>[][] startingPrimerReadsHeaders = new List<string>[fileStride][];
            List<string>[][] startingPrimerReadsReads = new List<string>[fileStride][];
            List<string>[][] endingPrimerReadsHeaders = new List<string>[fileStride][];
            List<string>[][] endingPrimerReadsReads = new List<string>[fileStride][];
            // these reads will be 'selected' prior to the final filtering pass, just to make sure they're
            // not missed, and primerReadsIdx is used to avoid selecting them twice. 
            Dictionary<int, uint>[][] primerReadsIdx = new Dictionary<int, uint>[fileStride][];

            int[][] startingPrimerReadsCount = new int[fileStride][];
            int[][] endingPrimerReadsCount = new int[fileStride][];
            int initialReadListLength = 100000;
            if (prefilteredReads)
                initialReadListLength = readsToFilter[0].Count / 2;

            for (int fip = 0; fip < fileStride; fip++)
            {
                startingPrimerReadsHeaders[fip] = new List<string>[2];
                startingPrimerReadsReads[fip] = new List<string>[2];
                endingPrimerReadsHeaders[fip] = new List<string>[2];
                endingPrimerReadsReads[fip] = new List<string>[2];
                startingPrimerReadsCount[fip] = new int[2];
                endingPrimerReadsCount[fip] = new int[2];
                foreach (int d in directions)
                {
                    startingPrimerReadsHeaders[fip][d] = new List<string>(initialReadListLength);
                    startingPrimerReadsReads[fip][d] = new List<string>(initialReadListLength);
                    endingPrimerReadsHeaders[fip][d] = new List<string>(initialReadListLength);
                    endingPrimerReadsReads[fip][d] = new List<string>(initialReadListLength);
                }

                primerReadsIdx[fip] = new Dictionary<int, uint>[partitionCount[fip]];

                // Get an initial set of reads with primers near their starts. These will be trimmed and used to generate the initial region filter set.
                // The terminal primers are also found and used to populate the ending kMer set (along with the terminal primers themselves)
                if (prefilteredReads)
                    FindReadsWithPrimers(headersToFilter[fip], readsToFilter[fip], primers, mismatchesFP, mismatchesRP, allowedBases,
                                         fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength, kMerSize, kMerSize + startingReadTrim,
                                         allPrimerCores, forwardPrimersFwd, forwardPrimersRC, reversePrimersFwd, reversePrimersRC,
                                         startingPrimerReadsHeaders[fip], startingPrimerReadsReads[fip], endingPrimerReadsHeaders[fip], endingPrimerReadsReads[fip],
                                         primerReadsIdx[fip], readsAlreadyTiled[fip], primersFound);
                else
                    FindReadsWithPrimers(readsInPartitions[fip], WGSReadsFNs[fip], primers, mismatchesFP, mismatchesRP, allowedBases,
                                         fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength, kMerSize, kMerSize + startingReadTrim,
                                         allPrimerCores, forwardPrimersFwd, forwardPrimersRC, reversePrimersFwd, reversePrimersRC,
                                         startingPrimerReadsHeaders[fip], startingPrimerReadsReads[fip], endingPrimerReadsHeaders[fip], endingPrimerReadsReads[fip],
                                         primerReadsIdx[fip], readsAlreadyTiled[fip], (fileStride > 1 ? readsAlreadyTiled[fip == 0 ? 1 : 0] : null), primersFound, threadsWanted);

                startingPrimerReadsCount[fip][fwd] = startingPrimerReadsReads[fip][fwd].Count;
                endingPrimerReadsCount[fip][fwd] = endingPrimerReadsReads[fip][fwd].Count;
                startingPrimerReadsCount[fip][rvs] = startingPrimerReadsReads[fip][rvs].Count;
                endingPrimerReadsCount[fip][rvs] = endingPrimerReadsReads[fip][rvs].Count;
            }

            allPrimerCores = null;

            int totalStartingPrimerReads = 0;
            int totalEndingPrimerReads = 0;
            for (int fip = 0; fip < fileStride; fip++)
            {
                totalStartingPrimerReads += startingPrimerReadsCount[fip][fwd] + startingPrimerReadsCount[fip][rvs];
                totalEndingPrimerReads += endingPrimerReadsCount[fip][fwd] + endingPrimerReadsCount[fip][rvs];
            }

            if (fileStride == 1)
                Console.WriteLine("found " + totalStartingPrimerReads + " (" + startingPrimerReadsCount[0][fwd] + "F+" + startingPrimerReadsCount[0][rvs] + "R) reads matching starting primers, " +
                                             totalEndingPrimerReads + " (" + endingPrimerReadsCount[0][fwd] + "R'+" + endingPrimerReadsCount[0][rvs] + "F') matching ending primers");
            else
                Console.WriteLine("found " + totalStartingPrimerReads + " (" + startingPrimerReadsCount[0][fwd] + "F+" + startingPrimerReadsCount[0][rvs] + "R/" + startingPrimerReadsCount[1][fwd] + "F+" + startingPrimerReadsCount[1][rvs] + "R) reads matching starting primers, " +
                                             totalEndingPrimerReads + " (" + endingPrimerReadsCount[0][fwd] + "R'+" + endingPrimerReadsCount[0][rvs] + "F'/" + endingPrimerReadsCount[1][fwd] + "R'+" + endingPrimerReadsCount[1][rvs] + "F') matching ending primers");

            if (totalStartingPrimerReads == 0)
            {
                Console.WriteLine("No starting primer reads found - exiting.");
                return;
            }

            if (savePrimers)
            {
                StreamWriter primersFile = new StreamWriter(extendedReadsPrefix + "_primers.txt");
                string[] primerType = new string[] { "FP", "RP", "FP'", "RP'" };
                for (int pt = FPpt; pt <= rcRPpt; pt++)
                {
                    string[] pv = new string[primersFound[pt].Count];
                    int[] pc = new int[primersFound[pt].Count];
                    primersFound[pt].Keys.CopyTo(pv, 0);
                    primersFound[pt].Values.CopyTo(pc, 0);
                    Array.Sort<int, string>(pc, pv);
                    Array.Reverse(pv);
                    Array.Reverse(pc);

                    primersFile.WriteLine(primerType[pt]);
                    for (int j = 0; j < pv.Length; j++)
                        primersFile.WriteLine(pv[j] + "\t" + kMers.ReverseComplement(pv[j]) + "\t" + pc[j]);
                }
                primersFile.Close();
            }

            // regionFilter: An initial set of filtering kMers (starting from any detected starting primer and continuing at single-base increments)
            // The filter is built by scanning from the starts of both the R1 and R2 reads (if paired) and only those
            // kMers present in both read sets will be kept for the final filtering pass (-strict).
            // The kMers in these tables are packed kMerLength strings (longer than minimum primer length)
            // This filter is then used for a final pass over all the reads - similar to FilterReads
            HashSet<ulong>[][] regionFilter = new HashSet<ulong>[fileStride][]; // [fip][dir]
            // the set of 48-mer hashed kMers that precede (and include the terminal primers). 
            // augmented as new starting-derived reads are found. Used to stop building the region filter when the starting and ending kMer clouds run into each other
            HashSet<ulong>[][] endingFilter = new HashSet<ulong>[fileStride][];

            for (int fip = 0; fip < fileStride; fip++)
            {
                regionFilter[fip] = new HashSet<ulong>[2];
                endingFilter[fip] = new HashSet<ulong>[2];
            }

            // The detection of overlaps is done using a combination of kMer overlaps and contexts
            // kMers at the start (or end) of a read are first matched against the current (directional) region filter
            // and reads with matching kMers are then checked to see if the kMer appears in the correct context
            // Simple kMer checks, and simple single-length (40, 48) proved unable to handle COI primers on WGS data
            // Obvious data structure is kMer --> distinct contexts (as a kMer can have multiple contexts) and these contexts
            // are effectively the remainder of the read in which the kMer appeared. This is far too expensive to implement
            // so an approximation is used. Context lengths are quantised, and the presence
            // of a kMer at a context length says that this length context exists and should be checked.
            // These data structures are common to both files in a pair to better find low abundance reads.

            // Filter kMerContexts are distinct (non-overlapping), and maintained that way by ExtendFilterSet.
            // contextExists - are contexts known to exist for a kMer at a given length? [dir][lengthIdx]
            HashSet<ulong>[][] filterContextExists = new HashSet<ulong>[2][];
            // and the contexts themselves - hashed rest-of-reads, including the initial kMer
            HashSet<ulong>[][] filterContexts = new HashSet<ulong>[2][];
            // kMer counts table used to detect and avoid likely error kMers
            Dictionary<ulong, int>[] kMerCounts = new Dictionary<ulong, int>[2];

            // what context lengths are we interested in? 
            List<int> filterContextLengths = new List<int>(longestReadLength);
            for (int cl = shortestContextLength; cl <= longestReadLength; cl += filterContextStride)
                filterContextLengths.Add(cl);

            // adapterTrap is a small set of kMers thought be be from errant adapter seqs
            Dictionary<ulong, int> adapterTrap = new Dictionary<ulong, int>();

            foreach (int d in directions)
            {
                for (int fip = 0; fip < fileStride; fip++)
                {
                    regionFilter[fip][d] = new HashSet<ulong>();
                    endingFilter[fip][d] = new HashSet<ulong>();
                }
                kMerCounts[d] = new Dictionary<ulong, int>(100000);
                filterContextExists[d] = new HashSet<ulong>[filterContextLengths.Count];
                filterContexts[d] = new HashSet<ulong>[filterContextLengths.Count];
                for (int li = 0; li < filterContextLengths.Count; li++)
                {
                    filterContextExists[d][li] = new HashSet<ulong>();
                    filterContexts[d][li] = new HashSet<ulong>();
                }
            }

            // keep extending filter until both files have been exhausted
            bool[] fileFinished = new bool[fileStride];
            // how many files have finished
            int finishedFiles = 0;
            // stop filter extension loop once both threads have reported they're finished
            bool stopExtending = false;
            // how many extension iterations
            int filterIterations = 0;

            // we're told how long to expect the between-primers region to be so don't stop the filter iterations too early
            int minIterations = minAmpliconLength > 0 ? ((minAmpliconLength / 2) / (longestReadLength / 2)) : 0;
            // and don't go on too long if we know how long the amplicons should be
            int maxIterations = maxAmpliconLength > 0 ? (maxAmpliconLength / (longestReadLength / 2)) : int.MaxValue;

            StreamWriter[] matching = new StreamWriter[fileStride];
            if (saveMatches)
            {
                for (int fip = 0; fip < fileStride; fip++)
                    matching[fip] = new StreamWriter(extendedReadsPrefix + "_matching_" + (fip + 1) + ".txt");
            }

            int readsExpected = InitialiseIterativeFilters(startingPrimerReadsHeaders, startingPrimerReadsReads, endingPrimerReadsHeaders, endingPrimerReadsReads, primerReadsIdx, forwardPrimerLength, reversePrimerLength,
                                                           regionFilter[0], endingFilter[0], filterContextExists, filterContexts, kMerCounts, kMerSize, minDepth, adapterTrap, matching);
            int initialReadsExpected = readsExpected;

            // regionsFilter & endingFilter start off the same for both files in a pair
            if (fileStride > 1)
            {
                foreach (int d in directions)
                {
                    regionFilter[1][d] = new HashSet<ulong>(regionFilter[0][d]);
                    endingFilter[1][d] = new HashSet<ulong>(endingFilter[0][d]);
                }
            }

            // counters etc used to terminate filter extension loop
            // kept out here as filter extension is done pair-at-a-time
            bool[] foundEndingRegion = new bool[fileStride];
            bool[] bigDrop = new bool[fileStride];
            int[] totalEndingReadsFound = new int[fileStride];
            int[] totalMatchingReadsFound = new int[fileStride];
            int[] previousMatchingReadsFound = new int[fileStride];
            int[] rafNearZeroCount = new int[fileStride];

            List<string>[] matchingReads = new List<string>[2];
            int[] readsAddingToFilter = new int[2];
            List<string>[] endingReadsFound = new List<string>[2];
            int[] endingReadsCount = new int[2];
            HashSet<ulong>[] kMersFromSelectedReads = new HashSet<ulong>[2];

            foreach (int d in directions)
            {
                matchingReads[d] = new List<string>(10000);
                endingReadsFound[d] = new List<string>(10000);
                kMersFromSelectedReads[d] = new HashSet<ulong>();
            }

            // The region filter is populated by starting with the set of starting reads (fwd and rvs) and finding reads whose starts
            // overlap already scanned reads. The first iteration finds overlaps with the starting reads, and these newly found reads 
            // are then tiled in preparation for the next iteration. The loop finishes when the newly found reads run into reads found from the other direction. 
            // do this across both the files in a pair at each iteration to extract from both files at an equal pace (and sharing contexts etc across both files)

            while (!stopExtending)
            {
                filterIterations++;

                for (int fip = 0; fip < fileStride; fip++)
                {
                    if (fileFinished[fip])
                        continue;

                    int droppedLowComplexityReads;

                    // find a new set of reads that start with kMers already in the region filter (overlapping reads)
                    if (prefilteredReads)
                        FindMatchingReads(filterIterations, readsToFilter[fip], readsAlreadyTiled[fip][0], kMerSize, dropLowComplexity, regionFilter[fip], filterContextLengths, filterContextExists, filterContexts, matchingReads, out droppedLowComplexityReads, threadsWanted);
                    else
                        FindMatchingReads(filterIterations, readsInPartitions[fip], WGSReadsFNs[fip], readsAlreadyTiled[fip], kMerSize, dropLowComplexity, regionFilter[fip], filterContextLengths, filterContextExists, filterContexts, matchingReads, out droppedLowComplexityReads, threadsWanted);

                    int matchingReadsCount = matchingReads[fwd].Count + matchingReads[rvs].Count;
                    totalMatchingReadsFound[fip] += matchingReadsCount;

                    // add these newly-found reads to the growing filter
                    ExtendFilterSet(fip, matchingReads, adapterTrap, regionFilter[fip], kMerSize, dropLowComplexity, endingFilter[fip], kMerCounts, kMersFromSelectedReads, filterContextLengths, filterContextExists, filterContexts, readsAddingToFilter, endingReadsCount, endingReadsFound);
                    totalEndingReadsFound[fip] += endingReadsCount[fwd] + endingReadsCount[rvs];

                    if (saveMatches)
                    {
                        foreach (int d in directions)
                        {
                            matching[fip].WriteLine((d == 0 ? "fwd" : "rvs") + ": matching #" + filterIterations);
                            for (int r = 0; r < matchingReads[d].Count; r++)
                            {
                                string read = matchingReads[d][r];
                                matching[fip].WriteLine(read);
                            }
                        }
                        foreach (int d in directions)
                        {
                            matching[fip].WriteLine((d == 0 ? "fwd" : "rvs") + ": ending #" + filterIterations);
                            foreach (string read in endingReadsFound[d])
                                matching[fip].WriteLine(read);
                        }
                    }

                    // Stop when we either found no more reads to tile - or if we appear to be finding an increasing number of reads to tile
                    // The number of reads to tile should be approximately flat while traversing the inter-primer region and then drop rapidly
                    // as we get to the end of the region. An increase in this number indicates that we've overshot the terminating primer.

                    // kludge to handle the spike-in amplicon data which is largely missing its reverse primers. This results in the initial estimate for readsExpected being far too low. Only adjust if the jump is big though
                    int minMatchingCount = Math.Min(matchingReads[fwd].Count, matchingReads[rvs].Count);
                    //  not yet in ending region           much bigger than initial estimate            and 25% larger than any previous revision
                    if (totalEndingReadsFound[fip] == 0 && minMatchingCount > 5 * initialReadsExpected && minMatchingCount > 5 * readsExpected / 4)
                    {
                        Console.WriteLine("adjusted readsExpected " + readsExpected + " --> " + minMatchingCount);
                        readsExpected = minMatchingCount;
                    }

                    // say we're into the ending region when the number of ending reads is near the expected reads
                    foundEndingRegion[fip] |= totalEndingReadsFound[fip] > 3 * readsExpected / 4;

                    int minReadsAddingToFilter = Math.Min(readsAddingToFilter[fwd], readsAddingToFilter[rvs]);
                    if (minReadsAddingToFilter <= readsExpected / 100)
                        rafNearZeroCount[fip]++;
                    int minEndingReadsCount = Math.Min(endingReadsCount[fwd], endingReadsCount[rvs]);
                    bool tinyEndCount = foundEndingRegion[fip] && minEndingReadsCount <= totalEndingReadsFound[fip] / 100;
                    bool terReached = foundEndingRegion[fip] && totalEndingReadsFound[fip] > readsExpected && !VeryClose(totalEndingReadsFound[fip], readsExpected);
                    bool tooManyReads = foundEndingRegion[fip] && minMatchingCount > readsExpected * 2;
                    bigDrop[fip] = (bigDrop[fip] && CloseOrLower(matchingReadsCount, previousMatchingReadsFound[fip])) || (matchingReadsCount < previousMatchingReadsFound[fip] / 2);
                    fileFinished[fip] = (rafNearZeroCount[fip] > 1 || tinyEndCount || (terReached && bigDrop[fip])) && (filterIterations > minIterations || minReadsAddingToFilter == 0) ||
                                        tooManyReads ||
                                        filterIterations > maxIterations ||
                                        (finishedFiles > 0 && minEndingReadsCount == 0);
                    if (fileFinished[fip])
                        finishedFiles++;
                    string finishedMsg = "";
                    if (fileFinished[fip])
                        finishedMsg = " - finished.";
                    if (tooManyReads)
                        finishedMsg = " - finished (too many reads).";
                    if (filterIterations > maxIterations)
                        finishedMsg = " - finished (max length).";

                    Console.WriteLine((fip + 1) + "-" + filterIterations + ". matched reads: " + matchingReads[fwd].Count + "+" + matchingReads[rvs].Count + "=" + matchingReadsCount +
                                                                           " (" + readsAddingToFilter[fwd] + "+" + readsAddingToFilter[rvs] + " adding, " +
                                                                           endingReadsCount[fwd] + "+" + endingReadsCount[rvs] + "=>" + totalEndingReadsFound[fip] + " ending, " + readsExpected + " expected)" + finishedMsg);
                    //Console.WriteLine("stop=" + fileFinished[fip] + ", fer=" + foundEndingRegion[fip] + ", rafNearZero=" + rafNearZeroCount[fip] + ", tinyEnd=" + tinyEndCount + ", end=" + endingReadsCount[fwd] + "+" + endingReadsCount[rvs] + "=>" + totalEndingReadsFound[fip] +
                    //                  ", terReached=" + terReached + ", tooManyReads=" + tooManyReads + ", bigDrop=" + bigDrop[fip] + ", terf=" + totalEndingReadsFound[fip] + ", pmrf=" + previousMatchingReadsFound[fip] + ", maxIters=" + (filterIterations > maxIterations));
                    //if (stopExtending)
                    //  Debugger.Break();

                    // update the region filter now that we know whether or not this is the last iteration
                    foreach (int d in directions)
                    {
                        // only keep kMers from 'ending' reads if we terminated because things looked like they had got out of hand
                        if (tooManyReads)
                        {
                            //Console.WriteLine("tooManyReads: " + matchingReadsCount + " matched, " + readsExpected + " expected");
                            // replace the kMers from all the selected reads with just those from the ending reads
                            TileReadsForMers(endingReadsFound[d], kMersFromSelectedReads[d], kMerSize);
                        }

                        // add these new kMers to the region filter
                        regionFilter[fip][d].UnionWith(kMersFromSelectedReads[d]);
                    }

                    previousMatchingReadsFound[fip] = matchingReadsCount;
                } // for each file in the pair 

                stopExtending = true;
                for (int fip = 0; fip < fileStride; fip++)
                    stopExtending &= fileFinished[fip];

            } // while still growing filter sets

            // something is probably horribly wrong with the primers if we found no overlap between the primer clouds
            int totalTotalEndingReadsFound = 0;
            for (int fip = 0; fip < fileStride; fip++)
                totalTotalEndingReadsFound += totalEndingReadsFound[fip];
            if (totalTotalEndingReadsFound == 0)
                Console.WriteLine("WARNING - no overlap found in primer-derived reads. Results may be unreliable.");

            if (saveMatches)
                for (int fip = 0; fip < fileStride; fip++)
                    matching[fip].Close();

            matchingReads = null;
            endingReadsFound = null;
            kMersFromSelectedReads = null;
            endingFilter = null;
            filterContextExists = null;
            filterContexts = null;
            filterContextLengths = null;
            kMerCounts = null;

            // build final region filter used to find the reads that will go into phase 2
            HashSet<ulong> finalRegionFilter = new HashSet<ulong>();
            for (int fip = 0; fip < fileStride; fip++)
            {
                // merge fwd and rvs kMers clouds --> into fwd for convenience
                regionFilter[fip][fwd].UnionWith(regionFilter[fip][rvs]);
                // and add the fwd+rvs from each strand to the final region filter
                finalRegionFilter.UnionWith(regionFilter[fip][fwd]);
            }

            Console.WriteLine("initial region filter contains " + finalRegionFilter.Count + " kMers");
            //Console.WriteLine("initial ending filter contains " + finalEndingFilter.Count + " kMers");

            //foreach (ulong adk in adapterTrap.Keys)
            //    if (finalRegionFilter.Contains(adk))
            //        Debugger.Break();

            // strict pairs - kMers must also exist (in either form) in both files in the pair
            if (strict && fileStride > 1)
            {
                List<ulong> kMersToDelete = new List<ulong>();
                foreach (ulong kMer in finalRegionFilter)
                {
                    //if (kMer == 0x1EC02E5F631E9FDF)
                    //    Debugger.Break();
                    ulong kMerRC = kMers.ReverseComplement(kMer, kMerSize);
                    bool presentInBoth = regionFilter[0][fwd].Contains(kMer) && (regionFilter[1][fwd].Contains(kMer) || regionFilter[1][fwd].Contains(kMerRC)) ||
                                         regionFilter[1][fwd].Contains(kMer) && (regionFilter[0][fwd].Contains(kMer) || regionFilter[0][fwd].Contains(kMerRC));
                    if (!presentInBoth)
                        kMersToDelete.Add(kMer);
                }
                finalRegionFilter.ExceptWith(kMersToDelete);
                //StreamWriter deleted = new StreamWriter("deleted_" + filteredReadsTag + ".txt");
                //foreach (ulong dm in kMersToDelete)
                //    deleted.WriteLine(kMers.ExpandMer(dm, kMerSize));
                //deleted.Close();
            }

            // finally ensure RC forms are present for all kMers by adding them if necessary
            HashSet<ulong> missingRCForms = new HashSet<ulong>();
            foreach (ulong kMer in finalRegionFilter)
            {
                ulong kMerRC = kMers.ReverseComplement(kMer, kMerSize);
                if (!finalRegionFilter.Contains(kMerRC))
                    missingRCForms.Add(kMerRC);
            }
            finalRegionFilter.UnionWith(missingRCForms);

            regionFilter = null;
            missingRCForms = null;

            Console.WriteLine("final region filter contains " + finalRegionFilter.Count + " kMers");
            //Console.WriteLine("final ending filter contains " + finalEndingFilter.Count + " kMers");

            // now do a final pass though the starting reads, matching the starts of the reads against the region filter
            // and writing out the selected reads if requested. The selected reads are saved for the extension phase.   
            List<string> selectedHeaders = new List<string>(5000);
            List<string> selectedReads = new List<string>(5000);
            Dictionary<int, int> readPairs = new Dictionary<int, int>(5000);

            // selectedReads are merged across all files/pairs, so remember which ones have already been written
            int selectedReadsWritten = 0;
            StreamWriter filteredReads = null;
            if (filteredReadsTag != null)
            {
                string filteredReadsFN = "Kelpie_filtered_reads_" + filteredReadsTag + ".fa";
                filteredReads = new StreamWriter(filteredReadsFN);
            }

            // pre-select the primer-containing reads (these will have been trimmed to primer edges, unlike their WGS originals)
            for (int fip = 0; fip < fileStride; fip++)
            {
                SelectPrimerReads(startingFwd, primerReadsIdx[fip], startingPrimerReadsHeaders[fip][fwd], startingPrimerReadsReads[fip][fwd], selectedHeaders, selectedReads);
                SelectPrimerReads(endingFwd, primerReadsIdx[fip], endingPrimerReadsHeaders[fip][fwd], endingPrimerReadsReads[fip][fwd], selectedHeaders, selectedReads);
                SelectPrimerReads(startingRvs, primerReadsIdx[fip], startingPrimerReadsHeaders[fip][rvs], startingPrimerReadsReads[fip][rvs], selectedHeaders, selectedReads);
                SelectPrimerReads(endingRvs, primerReadsIdx[fip], endingPrimerReadsHeaders[fip][rvs], endingPrimerReadsReads[fip][rvs], selectedHeaders, selectedReads);
            }

            if (filteredReadsTag != null)
            {
                WriteFilteredReads(filteredReads, selectedReadsWritten, selectedHeaders, selectedReads, kMerSize, "primer");
                selectedReadsWritten = selectedReads.Count;
            }

            startingPrimerReadsHeaders = null;
            startingPrimerReadsReads = null;
            endingPrimerReadsHeaders = null;
            endingPrimerReadsReads = null;

            // go through all of the reads and select/filter those that match the inter-primer kMers (primer-containing reads have already been selected)
            if (prefilteredReads)
                FinalFilterReads(primerReadsIdx, readsInPartitions, headersToFilter, readsToFilter, kMerSize, finalRegionFilter, selectedHeaders, selectedReads, longestReadLength, pairedReads, readPairs);
            else
                FinalFilterReads(primerReadsIdx, readsInPartitions, WGSReadsFNs, kMerSize, finalRegionFilter, selectedHeaders, selectedReads, longestReadLength, pairedReads, readPairs, threadsWanted);

            // sanity check the paired reads - and turn off pair-checking if the reads turn out not to have been paired
            if (pairedReads)
            {
                bool readsArePaired = CheckReadPairings(readPairs, selectedHeaders);

                if (!readsArePaired)
                {
                    readPairs.Clear();
                    if (pairedReads)
                        Console.WriteLine("'paired' reads don't seem to be paired. Not doing paired-read path resolution");
                }
            }

            // and write them out if requested
            if (filteredReadsTag != null)
            {
                WriteFilteredReads(filteredReads, selectedReadsWritten, selectedHeaders, selectedReads, kMerSize, "filtered");
                filteredReads.Close();
            }
            selectedReadsWritten = selectedReads.Count;

            headersToFilter = null;
            readsToFilter = null;
            finalRegionFilter = null;
            if (!prefilteredReads && !leaveTempFilesBehind)
            {
                for (int fip = 0; fip < fileStride; fip++)
                    foreach (string FN in WGSReadsFNs[fip])
                        File.Delete(FN);
                string metadataFN = WGSReadsFNPrefix + "_kept_metadata.txt";
                if (File.Exists(metadataFN))
                    File.Delete(metadataFN);
            }

            //endingPairs.Close();
            phase1ElapsedTime.Stop();
            Console.WriteLine("filtered and kept " + selectedReads.Count + "/" + totalReadsCount + " primer-region reads from " + FNs.Count + " files in " + phase1ElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s");

            // ====================
            // Read extension phase
            // ====================

            Stopwatch phase2ElapsedTime = new Stopwatch();
            phase2ElapsedTime.Start();

            // build the (packed) terminating primers (derived from the RC reverse primer only for the extension phase as all starting reads will be in F --> (R) direction - after reversing if necessary)
            HashSet<ulong>[] extensionTerminatingPrimers = new HashSet<ulong>[2];
            extensionTerminatingPrimers[head] = new HashSet<ulong>();
            extensionTerminatingPrimers[core] = new HashSet<ulong>();
            foreach (string rpRCh in terminatingPrimers[head])
            {
                ulong rpRChp = kMers.CondenseMer(rpRCh, rvsPrimerHeadLength);
                extensionTerminatingPrimers[head].Add(rpRChp);
            }
            foreach (string rpRCc in terminatingPrimers[core])
            {
                ulong rpRCcp = kMers.CondenseMer(rpRCc, primerCoreLength);
                extensionTerminatingPrimers[core].Add(rpRCcp);
            }

            Console.WriteLine("built terminating primer set for extension phase");

            // just the reads with forward primers (FP or FP') - our starting sequences for later extension
            List<int> readsWithStartingPrimers = new List<int>();
            // and remember which of these starting seqs need to be RCed before extension (those ending in FP')
            // (RC is deferred as long as possible to allow Denoising etc to always go in the direction of the read - with more errors at the RHS)
            HashSet<int> readsToBeRCed = new HashSet<int>();
            // and all the other reads (to be checked later for partial primers)
            List<int> nonStartingReads = new List<int>();
            // the reads starting with a primer (F or R)
            HashSet<int> readsStartingWithPrimer = new HashSet<int>();
            // the reads ending with a primer (F' or R')
            HashSet<int> readsEndingWithPrimer = new HashSet<int>();
            // actual primers found at the starts of the starting reads
            Dictionary<int, string> FPReadPrimers = new Dictionary<int, string>();
            // FP' reads are RCed so they can be used as starting reads
            HashSet<int> startingReadRCed = new HashSet<int>();
            // reads with bases other than ACGT
            HashSet<int> nonACGTReads = new HashSet<int>();
            // longest reads found by ExtractPrimerReads
            int longestSelectedRead;

            // a canonical kMers table
            Dictionary<ulong, int> kMerTable = new Dictionary<ulong, int>(100000);
            // save kMer contexts for a quantised range of lengths. Shorter ones have a better chance of being found, but are less likely to be distinct.
            List<int> contextLengths = new List<int>();
            Dictionary<int, Dictionary<ulong, int>> kMerContexts = new Dictionary<int, Dictionary<ulong, int>>();
            Dictionary<int, int> kMerContextDepths = new Dictionary<int, int>();
            // remember the start of forward primer reads so we can salvage reads that don't have complete primers at their start
            Dictionary<string, int> startsOfFPReads = new Dictionary<string, int>();
            // kMers that are only found in one form (n,0)
            HashSet<ulong> kMersWithZero = new HashSet<ulong>(100000);

            // partitioned reads (for parallelism)
            List<int> selectedPartitionStarts = new List<int>();
            List<int> selectedPartitionEnds = new List<int>();
            PartitionReadList(selectedReads, selectedPartitionStarts, selectedPartitionEnds);

            // find the reads containing primers. Reads ending with FP' are reversed. Reads starting with FP are trimmed. The starts of the pre-trimmed FP reads are saved in startsOfFPReads. 
            ExtractPrimerReads(selectedHeaders, selectedReads, selectedPartitionStarts, selectedPartitionEnds, startingReadTrim, forwardPrimerLength,
                               readsWithStartingPrimers, startingReadRCed, startsOfFPReads, FPReadPrimers, nonStartingReads, readsStartingWithPrimer, readsEndingWithPrimer, out longestSelectedRead, threadsWanted);
            // tile the filtered reads (both R1 and R2) and build the initial (noisy) kMerTable
            GenerateInitialkMerTable(selectedReads, startingReadRCed, selectedPartitionStarts, selectedPartitionEnds, kMerTable, kMerSize, nonACGTReads, kMersWithZero, threadsWanted);
            //Console.WriteLine("readsWithStartingPrimers="+readsWithStartingPrimers.Count+", nonStartingReads="+nonStartingReads.Count+ ", nonACGTReads="+ nonACGTReads.Count);

            // remove kMers that appear to be sequencing errors and save the depth stats for the starting reads
            ReadStats[] statsForReads = new ReadStats[selectedReads.Count];
            int kMersCulled = DeNoiseMerTable(selectedHeaders, selectedReads, selectedPartitionStarts, selectedPartitionEnds, readsStartingWithPrimer, readsEndingWithPrimer, longestReadLength, Math.Min(forwardPrimerLength, reversePrimerLength), 
                                              kMerSize, kMerTable, minDepth, kMersWithZero, statsForReads, log, threadsWanted);
            Console.WriteLine("removed " + kMersCulled + " kMers from kMerTable");
            int kMersAfterCulling = kMerTable.Count - kMersCulled;

            // fix up per-read stats now that denoise has finished (DN will have altered the kMer counts used to calculate the stats)
            RefreshStatsForSelectedReads(selectedReads, selectedPartitionStarts, selectedPartitionEnds, kMerTable, kMerSize, minDepth, readPairs, statsForReads, threadsWanted);

            // make contexts for all possible lengths (in contextStride quanta). 
            for (int pl = shortestContextSize; pl < longestSelectedRead; pl += contextStride)
            {
                // initialise each context length
                if (!kMerContexts.ContainsKey(pl))
                {
                    contextLengths.Add(pl);
                    kMerContexts.Add(pl, new Dictionary<ulong, int>(kMersAfterCulling));
                }
            }

            Stopwatch sw = new Stopwatch();
            sw.Start();
            Parallel.ForEach(contextLengths, pl =>
            {
                int coverage = GenerateContextTable(pl, selectedHeaders, selectedReads, selectedPartitionStarts, selectedPartitionEnds, nonACGTReads, kMerTable, kMerSize, statsForReads, kMerContexts, threadsWanted);
                kMerContextDepths.Add(pl, coverage);
            });
            sw.Stop();
            Console.WriteLine("generated " + kMerContexts.Count + " contexts in " + sw.Elapsed.TotalSeconds.ToString("F1") + "s");

            int totalContexts = 0;
            int largestContexts = 0;
            foreach (int ps in contextLengths)
            {
                totalContexts += kMerContexts[ps].Count;
                if (kMerContexts[ps].Count > largestContexts)
                    largestContexts = kMerContexts[ps].Count;
            }

            //if (generateStats)
            //{
            //    foreach (int ps in contextLengths)
            //        Console.WriteLine(ps + ": " + kMerContexts[ps].Count);
            //}

            // trim contexts table to avoid low coverage at higher lengths. The longest contexts become increasingly noisy
            // can can lead to good extensions being blocked unnecessarily. Reducing the maximum context length makes this trimming
            // somewhat redundant, but one of the test datasets has spike-in amplicons which resist this natural culling.

            // dropped by 30% - time to stop
            int contextCountThreshold = largestContexts * 7 / 10;
            int lastContextIdxUsed = 0;
            for (int i = kMerContexts.Count - 1; i >= 0; i--)
            {
                if (kMerContexts[contextLengths[i]].Count > contextCountThreshold)
                {
                    lastContextIdxUsed = i;
                    break;
                }
            }
            // remove the not-to-be-used context tables
            int initialContextLengthsCount = contextLengths.Count;
            int contextTablesWanted = lastContextIdxUsed + 1;
            if (contextTablesWanted != initialContextLengthsCount)
            {
                for (int i = contextTablesWanted; i < initialContextLengthsCount; i++)
                    kMerContexts.Remove(contextLengths[i]);
                contextLengths.RemoveRange(contextTablesWanted, initialContextLengthsCount - contextTablesWanted);
            }

            // stats on context length usage during extension phase
            int[] contextLengthStats = new int[kMerContexts.Count + statsContextsIdx];

            Console.WriteLine("finished building kMer tables: " + (kMerTable.Count - kMersCulled) + " " + kMerSize + "-mers; " +
                               totalContexts + " kMer contexts over " + kMerContexts.Count + "/" + initialContextLengthsCount + " context lengths");

            Console.WriteLine("found " + readsWithStartingPrimers.Count + " reads starting with full F primer");

            int minReadLength = shortestContextSize;

            // find any other reads that start with a partial forward primer 
            List<int> additionalStartingReads = new List<int>();
            if (log != null)
                log.WriteLine("partial starting reads");
            foreach (int r in nonStartingReads)
            {
                string read = selectedReads[r];

                if (read == null)
                    continue;

                //if (read == "GTCGTCTCCAATTAGGGCACCAGGGTGTCCTAATTCAGCTCGAATTAGAATACTTAATGAAGTACCAACTATTCCTGCCCAAGCTCCAAATATAAAATATAAAGTTCCAATGTCTTTATGGTTTG")
                //    Debugger.Break();
                //if (r == 1916)
                //    Debugger.Break();

                if (read.Length < minReadLength)
                    continue;

                // see if the start of the read is in our partial primers set
                string startFragment = read.Substring(0, shortestContextSize);
                if (startsOfFPReads.ContainsKey(startFragment))
                {
                    selectedReads[r] = read.Substring(startingReadTrim - startsOfFPReads[startFragment]);
                    additionalStartingReads.Add(r);
                    if (log != null)
                        log.WriteLine(selectedReads[r]);
                    continue;
                }
                // try the end of the read in case it overlaps with FP'
                startFragment = kMers.ReverseComplement(read.Substring(read.Length - shortestContextSize, shortestContextSize));
                if (startsOfFPReads.ContainsKey(startFragment))
                {
                    read = kMers.ReverseComplement(read);
                    selectedReads[r] = read.Substring(startingReadTrim - startsOfFPReads[startFragment]);
                    additionalStartingReads.Add(r);
                    if (log != null)
                        log.WriteLine(selectedReads[r]);
                    continue;
                }
            }
            startsOfFPReads = null;

            Console.WriteLine("found " + additionalStartingReads.Count + "/" + nonStartingReads.Count + " additional starting reads with partial primers");

            for (int a = 0; a < additionalStartingReads.Count; a++)
            {
                int r = additionalStartingReads[a];
                readsWithStartingPrimers.Add(r);
                RefreshStatsForRead(selectedReads[r], kMerTable, kMerSize, statsForReads[r].minDepthAllowed, statsForReads[r]);
            }

            int extendedShortStartingReads = 0;
            // build final set of starting reads - after cleaning, trimming, extending, dropping
            List<int> finalStartingReads = new List<int>();

            GenerateFinalStartingReads(readsWithStartingPrimers, selectedHeaders, selectedReads, statsForReads, finalStartingReads, minReadLength, longestReadLength, kMerSize, minDepth, kMerTable, contextLengths, kMerContexts, kMersWithZero, readPairs,
                                       out extendedShortStartingReads, log);

            Console.WriteLine("extended " + extendedShortStartingReads + " shorter starting reads");
            Console.WriteLine("extending " + finalStartingReads.Count + " reads starting with F primer");

            readsWithStartingPrimers = null;           // finished with the initial set of starting reads
            readsStartingWithPrimer = null;
            readsEndingWithPrimer = null;

            // sort the starting reads to make the log deterministic
            finalStartingReads.Sort();

            //if (log != null)
            //{
            //    foreach (int ri in finalStartingReads)
            //    {
            //        log.WriteLine(ri + selectedHeaders[ri]);
            //        log.WriteLine(selectedReads[ri]);
            //    }
            //}

            for (int i = 0; i < selectedReads.Count; i++)
            {
                string read = selectedReads[i];
                if (read == null)
                    continue;

                ulong kMerAtStart = 0;
                if (kMers.CondenseMer(read, 0, kMerSize, out kMerAtStart))
                {
                    if (!startsOfReads.ContainsKey(kMerAtStart))
                        startsOfReads.Add(kMerAtStart, new List<string>(10));
                    startsOfReads[kMerAtStart].Add(read);
                    string readRC = kMers.ReverseComplement(read);
                    kMers.CondenseMer(readRC, 0, kMerSize, out kMerAtStart);
                    if (!startsOfReads.ContainsKey(kMerAtStart))
                        startsOfReads.Add(kMerAtStart, new List<string>(10));
                    startsOfReads[kMerAtStart].Add(readRC);
                }
            }

            // start from a set of starting reads and incrementally grow them from the available kMers until they can be grown no further
            // --------------------------------------------------------------------------------------------------------------------------
            ConcurrentDictionary<int, string> extendedStartingReads = ExtendFromStartingReads(selectedHeaders, selectedReads, finalStartingReads, startingReadRCed, minReadLength, maxAmpliconLength, longestReadLength, statsForReads, kMerSize, kMerTable, contextLengths, kMerContexts, contextLengthStats, extensionTerminatingPrimers,
                                                                                    fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength, readPairs, threadsWanted, log);

            //if (log != null)
            //{
            //    log.WriteLine("final starting reads");
            //    foreach (int r in finalStartingReads)
            //        log.WriteLine(r + ": " + selectedHeaders[r]);
            //    foreach (KeyValuePair<int, int> kvp in readPairs)
            //        log.WriteLine(kvp.Key + " " + kvp.Value);
            //}

            // now have a set of final extended reads - almost ready for writing
            // -----------------------------------------------------------------

            // but first trim away the starting and terminating primers (and discard any reads that are too short or too long)
            Dictionary<string, int> discards = new Dictionary<string, int>();
            Dictionary<int, string> trimmedExtendedReads = new Dictionary<int, string>(finalStartingReads.Count);
            Dictionary<string, int> discardsR = new Dictionary<string, int>();
            Dictionary<int, string> TPsFound = new Dictionary<int, string>(finalStartingReads.Count);
            Dictionary<string, int> distinctGoodReads = new Dictionary<string, int>(finalStartingReads.Count);

            int fullyExtendedCount = TrimExtendedReads(extendedStartingReads, fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength, reversePrimersRC, minExtendedLength, maxExtendedLength, trimmedExtendedReads, distinctGoodReads, TPsFound, discards, discardsR, log);

            // and write out all the extended reads
            // ------------------------------------
            StreamWriter extendedReads = new StreamWriter(extendedReadsFN);

            if (derepAmplicons)
            {
                int r = 0;
                foreach (KeyValuePair<string, int> kvp in distinctGoodReads)
                {
                    r++;
                    SeqFiles.WriteRead(extendedReads, ">R" + r + ";size=" + kvp.Value + ";", kvp.Key, SeqFiles.formatFNA, 100);
                }
            }
            else
            {
                foreach (KeyValuePair<int, string> kvp in trimmedExtendedReads)
                {
                    int r = kvp.Key;
                    string extendedRead = kvp.Value;

                    string primersInRead = "";
                    if (savePrimers)
                    {
                        string FP = FPReadPrimers.ContainsKey(r) ? FPReadPrimers[r] : "partialPrimer";
                        string TP = TPsFound[r];
                        primersInRead = ";FP=" + FP + ";TP=" + TP;
                    }

                    SeqFiles.WriteRead(extendedReads, (">R" + r + primersInRead), extendedRead, SeqFiles.formatFNA, 100); 
                }
            }

            int droppedReads = 0;
            int lastDotIdx = extendedReadsFN.LastIndexOf('.');
            string discardsFN = extendedReads + "_discards";
            if (lastDotIdx > 0)
            {
                string prefix = extendedReadsFN.Substring(0, lastDotIdx);
                string suffix = extendedReadsFN.Substring(lastDotIdx);
                discardsFN = prefix + "_discards" + suffix;
            }
            StreamWriter discardedReads = new StreamWriter(discardsFN);
            int discardedNo = 1;
            foreach (KeyValuePair<string, int> kvp in discards)
            {
                string discardedExtendedRead = kvp.Key;
                int discardedExtendedCount = kvp.Value;
                int discardR = discardsR[discardedExtendedRead];
                string FP;
                if (!FPReadPrimers.TryGetValue(discardR, out FP))
                    FP = "partial";
                SeqFiles.WriteRead(discardedReads, (">D" + discardedNo + ";size=" + discardedExtendedCount) + ";FP=" + FP, discardedExtendedRead, SeqFiles.formatFNA, 100);
                droppedReads += discardedExtendedCount;
                discardedNo++;
            }
            Console.WriteLine("dropped " + droppedReads + " reads (too short after extension or no terminal primer found)");

            extendedReads.Close();
            discardedReads.Close();
            if (log != null)
                log.Close();

            phase2ElapsedTime.Stop();
            totalElapsedTime.Stop();

            string ampliconRange = "full length";
            if (minExtendedLength > 0 || maxExtendedLength > 0)
            {
                string minLengthBit = minExtendedLength > 0 ? (minExtendedLength + "bp") : "";
                string middleBit = minLengthBit != "" ? " - " : "";
                string maxLengthBit = maxExtendedLength > 0 ? (maxExtendedLength + "bp") : "";
                ampliconRange = minLengthBit + middleBit + maxLengthBit;
            }

            Console.WriteLine("extracted, extended & wrote " + trimmedExtendedReads.Count + " F->R reads (" + ampliconRange + ") to " + extendedReadsFN +
                              " in " + totalElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s. (" + phase1ElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s+" + phase2ElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s)");

            if (generateStats)
            {
                Console.WriteLine("read extensions:\t" + contextLengthStats[statsNextIdx]);
                Console.WriteLine("read cached:\t" + contextLengthStats[statsCachedResultIdx]);
                Console.WriteLine("branch cached:\t" + contextLengthStats[statsCachedBranchIdx]);
                Console.WriteLine("single next kMer:\t" + contextLengthStats[statsSingleIdx]);
                Console.WriteLine("looked downstream:\t" + contextLengthStats[statsDownstreamIdx]);
                Console.WriteLine("single good downstream:\t" + contextLengthStats[statsSingleDSTIdx]);
                Console.WriteLine("chose by lowest cost:\t" + contextLengthStats[statsCostIdx]);
                Console.WriteLine("chose in proportion by depth:\t" + contextLengthStats[statsRolledIdx]);
                Console.WriteLine("chose longest downstream:\t" + contextLengthStats[statsLongestIdx]);
                Console.WriteLine("chose using paired reads (<--): \t" + contextLengthStats[statsPairedReadsRIdx]);
                Console.WriteLine("chose using paired reads (-->): \t" + contextLengthStats[statsPairedReadsFIdx]);
                Console.WriteLine("chose using paired reads (<-->): \t" + contextLengthStats[statsPairedReadsFullIdx]);
                Console.WriteLine("chose using read coverage: \t" + contextLengthStats[statsReadsCoveredIdx]);

                for (int pl = 0; pl < kMerContexts.Count; pl++)
                    Console.WriteLine("single choice at length " + contextLengths[pl] + ":\t" + contextLengthStats[statsContextsIdx + pl]);
            }

            //Console.ReadLine();
        }

        private static bool CheckReadPairings(Dictionary<int, int> readPairs, List<string> selectedHeaders)
        {
            List<int> readNos = new List<int>(readPairs.Keys);
            for (int i = 0; i < readNos.Count; i += 100)
            {
                string h1 = selectedHeaders[readNos[i]];
                string h2 = selectedHeaders[readPairs[readNos[i]]];

                int h1Idx = h1.IndexOf(";pq=");
                if (h1Idx == -1)
                    h1Idx = h1.IndexOf(";FP");
                if (h1Idx == -1)
                    h1Idx = h1.IndexOf(";RP");
                if (h1Idx != -1)
                    h1 = h1.Substring(0, h1Idx);

                int h2Idx = h2.IndexOf(";pq=");
                if (h2Idx == -1)
                    h2Idx = h2.IndexOf(";FP");
                if (h2Idx == -1)
                    h2Idx = h2.IndexOf(";RP");
                if (h2Idx != -1)
                    h2 = h2.Substring(0, h2Idx);

                // normalise headers to remove 1/2 differences
                h1 = h1.Replace('2', '1');
                h2 = h2.Replace('2', '1');

                if (h1 != h2)
                    return false;
            }

            return true;
        }

        private static int InitialiseIterativeFilters(List<string>[][] startingPrimerReadsHeaders, List<string>[][] startingPrimerReadsReads, List<string>[][] endingPrimerReadsHeaders, List<string>[][] endingPrimerReadsReads, Dictionary<int, uint>[][] primerReadsIdx,
                                                      int forwardPrimerLength, int reversePrimerLength, HashSet<ulong>[] regionFilter, HashSet<ulong>[] endingFilter, HashSet<ulong>[][] filterContextExists, HashSet<ulong>[][] filterContexts,
                                                      Dictionary<ulong, int>[] kMerCounts, int kMerSize, int minDepth, Dictionary<ulong, int> adapterTrap, StreamWriter[] matching)
        {
            int fileStride = startingPrimerReadsReads.Length;
            int readsExpected = 0;

            // prepared (checked (reversed)) starting & ending reads
            List<string>[] preppedStartingReads = new List<string>[2];
            List<string>[] preppedEndingReads = new List<string>[2];
            Dictionary<ulong, int>[] endingReadsMerCounts = new Dictionary<ulong, int>[2];

            foreach (int d in directions)
            {
                preppedStartingReads[d] = new List<string>(startingPrimerReadsReads[0][fwd].Count * 2);
                preppedEndingReads[d] = new List<string>(endingPrimerReadsReads[0][fwd].Count * 2);
                endingReadsMerCounts[d] = new Dictionary<ulong, int>();
            }

            // prepare the starting/ending reads and count their kMers
            for (int fip = 0; fip < fileStride; fip++)
            {
                int[] preppedReadsCount = new int[2];
                // count F' and R' into this temporary kMer set. These new kMers will be used to validate the starting read kMers.
                Dictionary<ulong, int>[][] kMerCountsForChecking = new Dictionary<ulong, int>[2][];
                kMerCountsForChecking[fwd] = new Dictionary<ulong, int>[2];
                kMerCountsForChecking[rvs] = new Dictionary<ulong, int>[2];

                // prepare starting & ending reads (check primers, reverse as needed) - prepped reads will always start with the primer sequence
                preppedReadsCount[fwd] += PrepStartingReads(startingPrimerReadsHeaders[fip][fwd], startingPrimerReadsReads[fip][fwd], ScanOpts.asIs, ScanOpts.atStart, kMerSize, forwardPrimerLength, preppedStartingReads[fwd]);   // F
                preppedReadsCount[fwd] += PrepStartingReads(endingPrimerReadsHeaders[fip][rvs], endingPrimerReadsReads[fip][rvs], ScanOpts.RC, ScanOpts.atEnd, kMerSize, forwardPrimerLength, preppedEndingReads[rvs]);          // (F')
                preppedReadsCount[rvs] += PrepStartingReads(startingPrimerReadsHeaders[fip][rvs], startingPrimerReadsReads[fip][rvs], ScanOpts.asIs, ScanOpts.atStart, kMerSize, reversePrimerLength, preppedStartingReads[rvs]);   // R
                preppedReadsCount[rvs] += PrepStartingReads(endingPrimerReadsHeaders[fip][fwd], endingPrimerReadsReads[fip][fwd], ScanOpts.RC, ScanOpts.atEnd, kMerSize, reversePrimerLength, preppedEndingReads[fwd]);          // (R')

                // preliminary counting the kMers in the starting reads - used to detect adapters as these will be unbalanced (not found in both as-is and RC forms) - and will be in both the F and R 
                // all as-read reads - could contain adapter kMers
                kMerCountsForChecking[fwd][asRead] = new Dictionary<ulong, int>(10000);
                CountMersInReads(preppedStartingReads[fwd], kMerSize, ScanOpts.asIs, kMerCountsForChecking[fwd][asRead]);   // F
                kMerCountsForChecking[rvs][asRead] = new Dictionary<ulong, int>(10000);
                CountMersInReads(preppedStartingReads[rvs], kMerSize, ScanOpts.asIs, kMerCountsForChecking[rvs][asRead]);   // R
                // all RC-ed reads - will not contain adapters as these are always unbalanced
                kMerCountsForChecking[fwd][rced] = new Dictionary<ulong, int>(10000);
                CountMersInReads(preppedEndingReads[rvs], kMerSize, ScanOpts.asIs, kMerCountsForChecking[fwd][rced]);    // (F')
                kMerCountsForChecking[rvs][rced] = new Dictionary<ulong, int>(10000);
                CountMersInReads(preppedEndingReads[fwd], kMerSize, ScanOpts.asIs, kMerCountsForChecking[rvs][rced]);    // (R')

                Dictionary<ulong, int> allForward = new Dictionary<ulong, int>(kMerCountsForChecking[fwd][asRead].Count + kMerCountsForChecking[fwd][rced].Count);
                Dictionary<ulong, int> allReverse = new Dictionary<ulong, int>(kMerCountsForChecking[fwd][asRead].Count + kMerCountsForChecking[fwd][rced].Count);
                foreach (KeyValuePair<ulong, int> kvp in kMerCountsForChecking[fwd][asRead])
                    if (allForward.ContainsKey(kvp.Key))
                        allForward[kvp.Key] += kvp.Value;
                    else
                        allForward.Add(kvp.Key, kvp.Value);
                foreach (KeyValuePair<ulong, int> kvp in kMerCountsForChecking[fwd][rced])
                    if (allForward.ContainsKey(kvp.Key))
                        allForward[kvp.Key] += kvp.Value;
                    else
                        allForward.Add(kvp.Key, kvp.Value);
                foreach (KeyValuePair<ulong, int> kvp in kMerCountsForChecking[rvs][asRead])
                    if (allReverse.ContainsKey(kvp.Key))
                        allReverse[kvp.Key] += kvp.Value;
                    else
                        allReverse.Add(kvp.Key, kvp.Value);
                foreach (KeyValuePair<ulong, int> kvp in kMerCountsForChecking[rvs][rced])
                    if (allReverse.ContainsKey(kvp.Key))
                        allReverse[kvp.Key] += kvp.Value;
                    else
                        allReverse.Add(kvp.Key, kvp.Value);

                // cull any reads that seem to have adapters (adapters will be also found in other as-read reads but not in RCed reads
                CullReadsWithAdapters(kMerCountsForChecking[fwd][asRead], kMerCountsForChecking[fwd][rced], allReverse, kMerSize, preppedStartingReads[fwd], adapterTrap, matching[fip], "F-->");
                CullReadsWithAdapters(kMerCountsForChecking[rvs][asRead], kMerCountsForChecking[rvs][rced], allForward, kMerSize, preppedStartingReads[rvs], adapterTrap, matching[fip], "R-->");
                CullReadsWithAdapters(kMerCountsForChecking[fwd][rced], kMerCountsForChecking[fwd][asRead], allForward, kMerSize, preppedEndingReads[rvs], adapterTrap, matching[fip], "(F-->)");
                CullReadsWithAdapters(kMerCountsForChecking[rvs][asRead], kMerCountsForChecking[rvs][rced], allReverse, kMerSize, preppedEndingReads[fwd], adapterTrap, matching[fip], "(R-->)");

                List<ulong> adapterToDrop = new List<ulong>();
                foreach (ulong kMer in adapterTrap.Keys)
                {
                    //string fate = " kept";
                    if (adapterTrap[kMer] < 10)
                    {
                        //fate = " dropped";
                        adapterToDrop.Add(kMer);
                    }

                    //Console.WriteLine(kMers.ExpandMer(kMer, kMerSize) + " " + adapterTrap[kMer] + " " + fate);
                }
                foreach (ulong kMer in adapterToDrop)
                    adapterTrap.Remove(kMer);

                CountMersInReads(preppedStartingReads[fwd], kMerSize, ScanOpts.asIs, kMerCounts[fwd]);      // F
                CountMersInReads(preppedEndingReads[rvs], kMerSize, ScanOpts.asIs, kMerCounts[fwd]);     // (F')
                CountMersInReads(preppedStartingReads[rvs], kMerSize, ScanOpts.asIs, kMerCounts[rvs]);      // R
                CountMersInReads(preppedEndingReads[fwd], kMerSize, ScanOpts.asIs, kMerCounts[rvs]);     // (R')

                readsExpected = Math.Max(Math.Min(preppedReadsCount[fwd], preppedReadsCount[rvs]), readsExpected);
            }

            // initialise the region filter, contexts and ending filter
            for (int fip = 0; fip < fileStride; fip++)
            {
                // tile the extracted/trimmed primer-containing reads to the region filter
                int[] kMersAddedToRegionFilter = new int[2];
                // F and (F') for the initial region filter (fwd)
                kMersAddedToRegionFilter[fwd] += InitialiseRegionFilter(preppedStartingReads[fwd], regionFilter[fwd], filterContextExists[fwd], filterContexts[fwd], kMerSize, forwardPrimerLength, kMerCounts[fwd], minDepth);      // F
                kMersAddedToRegionFilter[fwd] += InitialiseRegionFilter(preppedEndingReads[rvs], regionFilter[fwd], filterContextExists[fwd], filterContexts[fwd], kMerSize, forwardPrimerLength, kMerCounts[fwd], minDepth);        // (F')
                // R and (R') for the initial region filter (rvs)
                kMersAddedToRegionFilter[rvs] += InitialiseRegionFilter(preppedStartingReads[rvs], regionFilter[rvs], filterContextExists[rvs], filterContexts[rvs], kMerSize, reversePrimerLength, kMerCounts[rvs], minDepth);      // R
                kMersAddedToRegionFilter[rvs] += InitialiseRegionFilter(preppedEndingReads[fwd], regionFilter[rvs], filterContextExists[rvs], filterContexts[rvs], kMerSize, reversePrimerLength, kMerCounts[rvs], minDepth);        // (R')

                // count the kMers in the ending reads
                CountMersInReads(preppedEndingReads[fwd], kMerSize, ScanOpts.RC, endingReadsMerCounts[fwd]);              // (R') 
                CountMersInReads(preppedStartingReads[rvs], kMerSize, ScanOpts.RC, endingReadsMerCounts[fwd]);            // R
                CountMersInReads(preppedEndingReads[rvs], kMerSize, ScanOpts.RC, endingReadsMerCounts[rvs]);              // (F')
                CountMersInReads(preppedStartingReads[fwd], kMerSize, ScanOpts.RC, endingReadsMerCounts[rvs]);            // F

                // and tile the ending reads to the ending filter
                int[] kMersAddedToEndingFilter = new int[2];
                // R' & (R) to end F & (F')
                kMersAddedToEndingFilter[rvs] += InitialiseEndingFilter(preppedEndingReads[fwd], endingFilter[fwd], kMerSize, reversePrimerLength, endingReadsMerCounts[fwd], minDepth);               // (R') --> R'
                kMersAddedToEndingFilter[rvs] += InitialiseEndingFilter(preppedStartingReads[rvs], endingFilter[fwd], kMerSize, reversePrimerLength, endingReadsMerCounts[fwd], minDepth);             // R --> (R)
                // F' & (F) to end R & (R')
                kMersAddedToEndingFilter[fwd] += InitialiseEndingFilter(preppedEndingReads[rvs], endingFilter[rvs], kMerSize, forwardPrimerLength, endingReadsMerCounts[rvs], minDepth);               // (F') --> F'
                kMersAddedToEndingFilter[fwd] += InitialiseEndingFilter(preppedStartingReads[fwd], endingFilter[rvs], kMerSize, forwardPrimerLength, endingReadsMerCounts[rvs], minDepth);             // F --> (F)

                // write initial (primer) matching reads if requested
                if (matching[fip] != null)
                {
                    foreach (int d in directions)
                    {
                        matching[fip].WriteLine("starting primer reads: " + fip + (d == fwd ? " F" : " R"));
                        for (int i = 0; i < preppedStartingReads[d].Count; i++)
                            matching[fip].WriteLine(preppedStartingReads[d][i]);
                        matching[fip].WriteLine("ending primer reads: " + fip + (d == fwd ? " (R')" : " (F')"));
                        for (int i = 0; i < preppedEndingReads[d].Count; i++)
                            matching[fip].WriteLine(preppedEndingReads[d][i]);
                    }
                }
            }

            //Console.WriteLine("initial filters: " + "s=" + (initialRegionFilter[fwd].Count + initialRegionFilter[rvs].Count) + ", e=" +
            //                                                (initialEndingFilter[fwd].Count + initialEndingFilter[rvs].Count) + "; target=" + readsExpected + " ending reads");

            return readsExpected;
        }

        // Removes sequences that look to have adapters from a set of starting reads. The adapters are only ever found in one form, never in RC form. They will be present in both the as-read (forward direction)
        // reads but not in the corresponding RCed reads. The 'unwanted' depths should always be zero and these depths are scanned from the end of the read until a zero depth is reached.
        // The read will then be further trimmed back to a solid kMer (both wanted strands and not in unwanted). 
        private static void CullReadsWithAdapters(Dictionary<ulong, int> kMerDepthsForRead, Dictionary<ulong, int> kMerDepthsWanted, Dictionary<ulong, int> kMerDepthsUnwanted, int kMerSize, List<string> startingReads, Dictionary<ulong, int> adapterTrap, StreamWriter matches, string direction)
        {
            ulong[] kMersFromRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int[] kMerDepthsFromRead = new int[1000];
            int[] kMerDepthsFromWanted = new int[1000];
            int[] kMerDepthsFromUnwanted = new int[1000];

            List<int> readsToCull = new List<int>();
            int readsTrimmed = 0;

            for (int r = 0; r < startingReads.Count; r++)
            {
                bool targetFound = false;
                //string targetRead = "AGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTAGGCGGGGTAATGGCCCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACCTGTCTCTTATACACATCTCCGAGCCCACGAGACGGTCGGCGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAA";
                //if (startingReads[r] == targetRead)
                //{
                //    targetFound = true;
                //    Debugger.Break();
                //}

                int kMersFound = kMers.GenerateMersFromRead(startingReads[r], kMerSize, ref kMersFromRead, ref kMersValid);
                if (kMersFound > kMerDepthsFromRead.Length)
                {
                    Array.Resize<int>(ref kMerDepthsFromRead, kMersFound + 500);
                    Array.Resize<int>(ref kMerDepthsFromWanted, kMersFound + 500);
                    Array.Resize<int>(ref kMerDepthsFromUnwanted, kMersFound + 500);
                }

                int nonZeroUnwantedCount = 0;
                for (int i = 0; i < kMersFound; i++)
                {
                    int kMerDepthForRead = 0;
                    kMerDepthsForRead.TryGetValue(kMersFromRead[i], out kMerDepthForRead);
                    kMerDepthsFromRead[i] = kMerDepthForRead;

                    int kMerDepthWanted = 0;
                    kMerDepthsWanted.TryGetValue(kMersFromRead[i], out kMerDepthWanted);
                    kMerDepthsFromWanted[i] = kMerDepthWanted;

                    int kMerDepthUnwanted = 0;
                    kMerDepthsUnwanted.TryGetValue(kMersFromRead[i], out kMerDepthUnwanted);
                    kMerDepthsFromUnwanted[i] = kMerDepthUnwanted;
                    if (kMerDepthUnwanted != 0)
                        nonZeroUnwantedCount++;

                    if (targetFound)
                        Console.WriteLine(i + " " + kMers.ExpandMer(kMersFromRead[i], kMerSize) + " " + kMerDepthForRead + ", " + kMerDepthWanted + ", " + kMerDepthUnwanted);
                }

                // only interested in trimming the read if there are any kMers that are unwanted
                if (nonZeroUnwantedCount == 0)
                    continue;

                // scan back until we hit a non-zero unwanted. Could be zeros at the end (strings of As)
                int lastNonZeroIdx = 0;
                bool unwantedFound = false;
                for (int i = kMersFound - 1; i >= 0; i--)
                {
                    if (kMerDepthsFromUnwanted[i] == 0 && unwantedFound)
                    {
                        lastNonZeroIdx = i + 1;
                        break;
                    }
                    if (kMerDepthsFromUnwanted[i] != 0)
                        unwantedFound = true;
                }

                // final check for good --> bad before trimming away likely adaptor sequence
                if (lastNonZeroIdx > 0)
                {
                    if (matches != null)
                        for (int i = 0; i < kMersFound; i++)
                            matches.WriteLine(i + " " + kMers.ExpandMer(kMersFromRead[i], kMerSize) + " " + kMerDepthsFromRead[i] + ", " + kMerDepthsFromWanted[i] + ", " + kMerDepthsFromUnwanted[i]);

                    // move back in the read until we get a solid kMer as there's typically a small length of random seq preceding the adaptor
                    int trimFromIdx = lastNonZeroIdx;
                    for (int i = lastNonZeroIdx - 1; i >= 0; i--)
                    {
                        if (kMerDepthsFromRead[i] > 0 && kMerDepthsFromWanted[i] > 0 && kMerDepthsFromUnwanted[i] == 0)
                        {
                            lastNonZeroIdx = i;
                            break;
                        }
                    }

                    // remember what we've trimmed in the trap
                    for (int i = lastNonZeroIdx + 1; i < kMersFound; i++)
                        if (adapterTrap.ContainsKey(kMersFromRead[i]))
                            adapterTrap[kMersFromRead[i]]++;
                        else
                            adapterTrap.Add(kMersFromRead[i], 1);

                    // trim off adaptor
                    string trimmedRead = startingReads[r].Substring(0, lastNonZeroIdx + kMerSize);
                    if (matches != null)
                        matches.WriteLine("trimmed " + direction + " [" + r + "] " + startingReads[r] + " -> " + trimmedRead);
                    startingReads[r] = trimmedRead;

                    readsTrimmed++;
                }
            }

            if (readsTrimmed > 0)
                Console.WriteLine("trimmed " + readsTrimmed + " " + direction + " starting reads containing adapters");
        }

        // reduces a set of reads to distinct reads+counts and orders them lexically
        private static List<string> ReduceToDistinctReads(List<string> reads, out List<int> counts)
        {
            // sort the starting reads so that replicates will be together, and sub-reads will come before their extended relatives
            reads.Sort();
            List<string> distinctReads = new List<string>(reads.Count / 2);
            counts = new List<int>(distinctReads.Count);
            int currentClusterCount = 0;

            for (int i = 0; i < reads.Count - 1; i++)
            {
                string currentRead = reads[i];
                string nextRead = reads[i + 1];

                //if (currentRead == "ACGTTCATTGCGTAGTTCGGGTAGTTTGCGCCTCTGAGTTCGCCAACACAACCTTCGTCTGATGCGACTGAGAACACGTTGGTAGCGCCGCACTGATCCTGCAGATCGTATCCGAAGA")
                //    Debugger.Break();

                // skip a read that is either the same or a substring of the following read
                if (currentRead.Length <= nextRead.Length &&
                    nextRead.Substring(0, currentRead.Length) == currentRead)
                {
                    currentClusterCount++;
                    continue;
                }

                // found a suitably distinct read so save it (and its count)
                distinctReads.Add(currentRead);
                currentClusterCount++;
                counts.Add(currentClusterCount);
                // reset the reads-in-cluster count
                currentClusterCount = 0;
            }

            if (reads.Count > 0)
            {
                distinctReads.Add(reads[reads.Count - 1]);
                currentClusterCount++;
                counts.Add(currentClusterCount);
            }

            return distinctReads;
        }

        private static bool PrimersQuiteDegenerate(List<string> primers)
        {
            //  H   A or C or T     not - G, H follows G in the alphabet    T G A   D
            //  B   G or T or C     not - A, B follows A                    C A G   V
            //  V   G or C or A     not-T(not - U), V follows U             C G T   B
            //  D   G or A or T     not-C, D follows C                      C T A   H
            //  N   any             aNy                                     N       N
            //  I   any             Inosine (a match-anything base)         N       N  

            bool degenerate = false;

            int ATCG = 0;
            int veryDegenerateBases = 0;

            foreach (string primer in primers)
            {
                foreach (char b in primer)
                {
                    if (b == 'A' || b == 'C' || b == 'G' || b == 'T')
                        ATCG++;
                    if (b == 'H' || b == 'B' || b == 'V' || b == 'D' || b == 'N' || b == 'I')
                        veryDegenerateBases++;
                }
                degenerate |= (ATCG - veryDegenerateBases) * 100 / primer.Length < 80;
            }

            return degenerate;
        }

        private static void SelectPrimerReads(uint primerType, Dictionary<int, uint>[] primerReadsIdx, List<string> primerReadsHeaders, List<string> primerReadsReads, List<string> selectedHeaders, List<string> selectedReads)
        {
            // for each partition
            for (int p = 0; p < primerReadsIdx.Length; p++)
            {
                // primerReadsIdx[partition] holds mappings from read# (in partition) to location in primerReadsHeaders/primerReadsReads lists. 
                // The list indexes are per-list, and there are 4 such lists. The top two bits of each entry says which list is the target (primerType).
                // These indexes will be rewritten recno --> selected idx. 

                // which records in primerReadsReads come from this primerType?
                HashSet<int> idxForType = new HashSet<int>(primerReadsReads.Count);
                // and the corresponding recordNos so we can rewrite primerReadsIdx
                Dictionary<int, int> recNosForType = new Dictionary<int, int>(primerReadsReads.Count);

                foreach (KeyValuePair<int, uint> kvp in primerReadsIdx[p])
                {
                    uint lrIdx = kvp.Value;
                    uint idxType = lrIdx & 0xc0000000;
                    if (idxType == primerType)
                    {
                        int primerListIdx = (int)(lrIdx & 0x3fffffff);
                        recNosForType.Add(primerListIdx, kvp.Key);
                        idxForType.Add(primerListIdx);
                    }
                }

                // for each selected primer read of the calling type (SelectPrimerReads is called 8 times - file stride, direction, forward+reverse]
                foreach (int pr in idxForType)
                {
                    string header = primerReadsHeaders[pr];
                    string read = primerReadsReads[pr];

                    selectedHeaders.Add(header);
                    selectedReads.Add(read);

                    // rewrite index now that primer read has been copied 
                    primerReadsIdx[p][recNosForType[pr]] = (uint)(selectedReads.Count - 1);

                    //if (header.StartsWith("@A00632:395:HLMHJDSX7:4:1623:28881:8578"))
                    //    Debugger.Break();
                }
            }
        }

        private static StreamWriter OpenNewWGSPartition(int currentPartitionNo, int fip, string tempDir, string WGSReadsFNPrefix, List<string> WGSReadsFNs)
        {
            string WGSFileName = tempDir + WGSReadsFNPrefix + "_" + currentPartitionNo + "_" + fip + ".tmp";
            WGSReadsFNs.Add(WGSFileName);
            return new StreamWriter(WGSFileName, false, Encoding.ASCII, 1000000);
        }

        private static void WriteFullUsage()
        {
            Console.WriteLine("usage: Kelpie [-h] [-t #threads] -f forwardPrimer -r reversePrimer readsFNP extendedReadsFN (" + version + ")");
            Console.WriteLine("       -h                - Write out this extended help and exit.");
            Console.WriteLine("       -threads nn       - max parallel threads to allow. max will use all available. Default is #cores/2.");
            Console.WriteLine("       -f forwardPrimer e.g. GTGYCAGCMGCCGCGGTAA");
            Console.WriteLine("       -r reversePrimer e.g. GGACTACNVGGGTWTCTAAT");
            Console.WriteLine("       -filtered         - WGS reads have been pre-filtered down to the genomic region of interest. (default)");
            Console.WriteLine("                           This (small) read set will be kept in memory for faster processing during filter construction.");
            Console.WriteLine("       -unfiltered       - WGS reads were not pre-filtered and will be re-read from files as needed.");
            Console.WriteLine("       -length nn[-mm]   - Expected length of targeted region/amplicon, including primers, e.g. 295 for 16SV4. Either a single length or a range can be specified.");
            Console.WriteLine("       readsFNP          - List of file names/patterns of reads to be processed. e.g. S1_270_reads_anonymous_R?_16S.fasta");
            Console.WriteLine("       extendedReadsFN   - File name for extended between-primer reads. e.g. S1_270_reads_anonymous_16S_V4.fasta");
            Console.WriteLine();
            Console.WriteLine("uncommonly needed options that can usually be ignored");
            Console.WriteLine("       -strict           - kMers used for final reads filtering must be found in both files of a pair. (default)");
            Console.WriteLine("       -loose            - Opposite of -strict, all kMers are use for filtering.");
            Console.WriteLine("       -paired|unpaired  - Use -paired if you have paired-end reads, and -unpaired if you want to force pairs of reads files to be processed individually.");
            Console.WriteLine("                           Default is to assume pairs of reads files are paired. Used when cleaning final kMer filter table with -strict.");
            Console.WriteLine("       -mismatches mm    - Allow up to mm mismatches in the 3' or 5' ends of a primer. Default is 1. Only increase for long imprecise primers. AKA -mm.");
            Console.WriteLine("       -min nn           - Minimum length accepted after assembly/extension, even if terminal primer not found. Default is that all extended 'amplicons' must finish with a terminating primer.");
            Console.WriteLine("       -mindepth nn      - kMers found fewer than this many times are dropped from the filtering table. Only added to help pull out full length 'amplicons' from a very deep dataset.");
            Console.WriteLine("       -qualtrim nn      - FASTQ reads are quality trimmed on the RHS using a sliding window until this qual score is reached. Default is 30. AKA -qt & -tq.");
            Console.WriteLine("       -noLCF            - no low-complexity filter. Low complexity kMers are usually not included in the region filters. This option lets them be included.");
            Console.WriteLine("       -save primerTag   - Save the filtered/trimmed between-primer reads, and add this tag when building the file names. e.g. v4 --> S1_270_reads_anonymous_R?_16S_v4.fasta");
            Console.WriteLine("       -primers          - Save the actual primer sequences found in the reads to XXXX_primers.txt.");
            Console.WriteLine("       -derep            - Dereplicate the extended reads (with size= annotations)");
            Console.WriteLine("       -tmp/-kept        - Used to improve efficiency of processing unfiltered WGS datasets. See documentation for usage.");
            Console.WriteLine("       -log              - Write debugging log.");
        }

        // *filtered reads* version
        private static void FindReadsWithPrimers(List<string> headersToFilter, List<string> readsToFilter, List<string>[] primers, int mismatchesAllowedFP, int mismatchesAllowedRP, Dictionary<char, HashSet<char>> allowedBases,
                                                 int fwdPrimerHeadLength, int rvsPrimerHeadLength, int primerCoreLength, int kMerSize, int shortestReadLength,
                                                 HashSet<ulong> allPrimerCores, HashSet<string>[] forwardPrimerFwd, HashSet<string>[] forwardPrimerRC,
                                                 HashSet<string>[] reversePrimerFwd, HashSet<string>[] reversePrimerRC,
                                                 List<string>[] startingPrimerReadsHeaders, List<string>[] startingPrimerReadsReads, List<string>[] endingPrimerReadsHeaders, List<string>[] endingPrimerReadsReads,
                                                 Dictionary<int, uint>[] primerReadsIdx, List<HashSet<int>> readsAlreadyTiled, Dictionary<string, int>[] primersUsed)
        {
            int firstPossibleCoreIdx = Math.Min(fwdPrimerHeadLength, rvsPrimerHeadLength);
            primerReadsIdx[0] = new Dictionary<int, uint>();

            ulong[] potentialCores = new ulong[1000];
            bool[] potentialCoresViable = new bool[1000];

            for (int r = 0; r < readsToFilter.Count; r++)
            {
                string header = headersToFilter[r];
                string read = readsToFilter[r];

                //if (header.StartsWith("@A00632:395:HLMHJDSX7:4:1623:28881:8578 1:N:0:ATGCGCAG+TACTCCTT"))
                //    Debugger.Break();
                //if (read == "GTTATTCCCTAACAACAGAGTTTTACGACCCGAAAGCCTACATAACTCACGCGGCGTTGCTCCGTCAGACTTTAGTCCATTGCGGAAGATTCCCTACTGCTGCCTCCCGAAGGAGTATGGGCCGAGTCTCAGTACCAGTGTGGCCGATCA")
                //    Debugger.Break();
                //if (r == 1916)
                //    Debugger.Break();

                int possibleCores = kMers.GenerateMersFromRead(read, primerCoreLength, ref potentialCores, ref potentialCoresViable);
                int coreIdx = -1;
                for (int i = firstPossibleCoreIdx; i < possibleCores; i++)
                    if (potentialCoresViable[i] && allPrimerCores.Contains(potentialCores[i]))
                    {
                        coreIdx = i;
                        break;
                    }

                if (coreIdx >= 0)
                {
                    string primerType = CheckReadForPrimers(coreIdx, header, read, primers, mismatchesAllowedFP, mismatchesAllowedRP, allowedBases, fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength,
                                                           forwardPrimerFwd, forwardPrimerRC, reversePrimerFwd, reversePrimerRC, kMerSize, out header, out read, primersUsed);

                    if (primerType != null)
                    {
                        uint primerListIdx = 0;

                        if (read.Length >= shortestReadLength)
                        {
                            if (primerType == ";FP")
                            {
                                primerListIdx = startingFwd | (uint)startingPrimerReadsReads[fwd].Count;
                                startingPrimerReadsHeaders[fwd].Add(header);
                                startingPrimerReadsReads[fwd].Add(read);
                            }
                            if (primerType == ";RP")
                            {
                                primerListIdx = startingRvs | (uint)startingPrimerReadsReads[rvs].Count;
                                startingPrimerReadsHeaders[rvs].Add(header);
                                startingPrimerReadsReads[rvs].Add(read);
                            }

                            if (primerType == ";RP'")
                            {
                                primerListIdx = endingFwd | (uint)endingPrimerReadsReads[fwd].Count;
                                endingPrimerReadsHeaders[fwd].Add(header);
                                endingPrimerReadsReads[fwd].Add(read);
                            }
                            if (primerType == ";FP'")
                            {
                                primerListIdx = endingRvs | (uint)endingPrimerReadsReads[rvs].Count;
                                endingPrimerReadsHeaders[rvs].Add(header);
                                endingPrimerReadsReads[rvs].Add(read);
                            }

                            // remember this is a primer-containing read
                            primerReadsIdx[0].Add(r, primerListIdx);
                            // and that we will have tiled it for kMers early
                            readsAlreadyTiled[0].Add(r);
                            //if (r == 3360)
                            //    Debugger.Break();
                        }
                    }
                }
            }
        }

        // *UNfiltered reads* version
        private static void FindReadsWithPrimers(List<int> readsInPartitions, List<string> WGSReadsFN, List<string>[] primers, int mismatchesAllowedFP, int mismatchesAllowedRP, Dictionary<char, HashSet<char>> allowedBases,
                                                 int fwdPrimerHeadLength, int rvsPrimerHeadLength, int primerCoreLength, int kMerSize, int shortestReadLength,
                                                 HashSet<ulong> allPrimerCores, HashSet<string>[] forwardPrimerFwd, HashSet<string>[] forwardPrimerRC,
                                                 HashSet<string>[] reversePrimerFwd, HashSet<string>[] reversePrimerRC,
                                                 List<string>[] startingPrimerReadsHeaders, List<string>[] startingPrimerReadsReads, List<string>[] endingPrimerReadsHeaders, List<string>[] endingPrimerReadsReads,
                                                 Dictionary<int, uint>[] primerReadsIdx, List<HashSet<int>> readsAlreadyTiled, List<HashSet<int>> readsAlreadyTiledOther, Dictionary<string, int>[] primersUsed, ParallelOptions threads)
        {
            Parallel.For(0, WGSReadsFN.Count, threads, p =>
            //for (int p = 0; p < WGSReadsFN.Count; p++)
            {
                StreamReader WGSReads = new StreamReader(WGSReadsFN[p]);
                int readsInPartition = readsInPartitions[p];

                List<string>[] localStartingPrimerReadsHeaders = new List<string>[2];
                List<string>[] localStartingPrimerReadsReads = new List<string>[2];
                List<string>[] localEndingPrimerReadsHeaders = new List<string>[2];
                List<string>[] localEndingPrimerReadsReads = new List<string>[2];
                foreach (int d in directions)
                {
                    localStartingPrimerReadsHeaders[d] = new List<string>(1000);
                    localStartingPrimerReadsReads[d] = new List<string>(1000);
                    localEndingPrimerReadsHeaders[d] = new List<string>(1000);
                    localEndingPrimerReadsReads[d] = new List<string>(1000);
                }
                Dictionary<string, int>[] localPrimersUsed = new Dictionary<string, int>[4];
                for (int i = 0; i < 4; i++)
                    localPrimersUsed[i] = new Dictionary<string, int>();
                primerReadsIdx[p] = new Dictionary<int, uint>();

                int firstPossibleCoreIdx = Math.Min(fwdPrimerHeadLength, rvsPrimerHeadLength);
                ulong[] potentialCores = new ulong[1000];
                bool[] potentialCoresViable = new bool[1000];

                for (int r = 0; r < readsInPartition; r++)
                {
                    string header = WGSReads.ReadLine();   // read the header for now
                    string read = WGSReads.ReadLine();     // and the sequence (single line), possibly empty (if Ns were present)

                    //if (header == "@SRR13051709.2990024 2990024 length=151")
                    //{
                    //    Console.WriteLine(read);
                    //    Debugger.Break();
                    //}
                    //if (read == "CCTTGTGCATGAGCATCGACAGGTACCAGCCGTTCAGGCCTGCGTTGGAGTTTGCAGTTGCAAGTCCGACTGAGATACCTGATGCTGCGGCAAGGACGGATGCACGCTGAGATCCACCGAAGTGGCTCTCGAGGGCGGTCGGGTATGTCTC")
                    //    Debugger.Break();

                    // look for a possible primer core
                    int possibleCores = kMers.GenerateMersFromRead(read, primerCoreLength, ref potentialCores, ref potentialCoresViable);
                    int coreIdx = -1;
                    for (int i = firstPossibleCoreIdx; i < possibleCores; i++)
                        if (potentialCoresViable[i] && allPrimerCores.Contains(potentialCores[i]))
                        {
                            coreIdx = i;
                            break;
                        }

                    // found  primer core within the read...
                    if (coreIdx >= 0)
                    {
                        string primerType = CheckReadForPrimers(coreIdx, header, read, primers, mismatchesAllowedFP, mismatchesAllowedRP, allowedBases, fwdPrimerHeadLength, rvsPrimerHeadLength, primerCoreLength,
                                                               forwardPrimerFwd, forwardPrimerRC, reversePrimerFwd, reversePrimerRC, kMerSize, out header, out read, localPrimersUsed);

                        if (primerType != null)
                        {
                            uint primerListIdx = 0;     // zero-based for now, adjusted during locked merge later

                            if (read.Length >= shortestReadLength)
                            {
                                if (primerType == ";FP")
                                {
                                    primerListIdx = startingFwd | (uint)localStartingPrimerReadsReads[fwd].Count;
                                    localStartingPrimerReadsHeaders[fwd].Add(header);
                                    localStartingPrimerReadsReads[fwd].Add(read);
                                }
                                if (primerType == ";RP")
                                {
                                    primerListIdx = startingRvs | (uint)localStartingPrimerReadsReads[rvs].Count;
                                    localStartingPrimerReadsHeaders[rvs].Add(header);
                                    localStartingPrimerReadsReads[rvs].Add(read);
                                }

                                if (primerType == ";RP'")
                                {
                                    primerListIdx = endingFwd | (uint)localEndingPrimerReadsReads[fwd].Count;
                                    localEndingPrimerReadsHeaders[fwd].Add(header);
                                    localEndingPrimerReadsReads[fwd].Add(read);
                                    if (readsAlreadyTiledOther != null)
                                        readsAlreadyTiledOther[p].Add(r);   // nature of paired reads mean that the other read for RP' and FP' will be outside the primer-bound region
                                }
                                if (primerType == ";FP'")
                                {
                                    primerListIdx = endingRvs | (uint)localEndingPrimerReadsReads[rvs].Count;
                                    localEndingPrimerReadsHeaders[rvs].Add(header);
                                    localEndingPrimerReadsReads[rvs].Add(read);
                                    if (readsAlreadyTiledOther != null)
                                        readsAlreadyTiledOther[p].Add(r);   // see above
                                }

                                // remember this is a primer-containing read
                                primerReadsIdx[p].Add(r, primerListIdx);
                                // and that it will be pre-loaded into the region filter so don't try to match when extending the filter
                                readsAlreadyTiled[p].Add(r);
                            }
                        }
                    }
                }

                WGSReads.Close();

                lock (startingPrimerReadsReads)
                {
                    // remember starting reads list offsets so these can be added to the local primerReadsIdx values
                    uint startingPrimerReadsFwdOffset = (uint)startingPrimerReadsReads[fwd].Count;
                    uint startingPrimerReadsRvsOffset = (uint)startingPrimerReadsReads[rvs].Count;
                    uint endingPrimerReadsFwdOffset = (uint)endingPrimerReadsReads[fwd].Count;
                    uint endingPrimerReadsRvsOffset = (uint)endingPrimerReadsReads[rvs].Count;

                    // save the primer reads found in this partition to the full set
                    foreach (int d in directions)
                    {
                        foreach (string header in localStartingPrimerReadsHeaders[d])
                            startingPrimerReadsHeaders[d].Add(header);
                        foreach (string read in localStartingPrimerReadsReads[d])
                            startingPrimerReadsReads[d].Add(read);
                        foreach (string header in localEndingPrimerReadsHeaders[d])
                            endingPrimerReadsHeaders[d].Add(header);
                        foreach (string read in localEndingPrimerReadsReads[d])
                            endingPrimerReadsReads[d].Add(read);
                    }

                    // adjust the primerReadsIdx values (from partition-local to global)
                    List<int> primerReadsIdxKeys = new List<int>(primerReadsIdx[p].Keys);
                    primerReadsIdxKeys.Sort();
                    for (int i = 0; i < primerReadsIdx[p].Count; i++)
                    {
                        int partitionReadIdx = primerReadsIdxKeys[i];
                        uint partitionPrimerListIdx = primerReadsIdx[p][partitionReadIdx];
                        uint primerReadSetKey = partitionPrimerListIdx & 0xc0000000;
                        uint idxWithinLocalSet = partitionPrimerListIdx & 0x3fffffff;
                        switch (primerReadSetKey)
                        {
                            case startingFwd:
                                idxWithinLocalSet += startingPrimerReadsFwdOffset;
                                break;
                            case startingRvs:
                                idxWithinLocalSet += startingPrimerReadsRvsOffset;
                                break;
                            case endingFwd:
                                idxWithinLocalSet += endingPrimerReadsFwdOffset;
                                break;
                            case endingRvs:
                                idxWithinLocalSet += endingPrimerReadsRvsOffset;
                                break;
                        }
                        partitionPrimerListIdx = primerReadSetKey | idxWithinLocalSet;
                        primerReadsIdx[p][partitionReadIdx] = partitionPrimerListIdx;
                    }

                    // save the primers found in this partition
                    for (int i = 0; i < 4; i++)
                        foreach (KeyValuePair<string, int> kvp in localPrimersUsed[i])
                        {
                            if (primersUsed[i].ContainsKey(kvp.Key))
                                primersUsed[i][kvp.Key] += kvp.Value;
                            else
                                primersUsed[i].Add(kvp.Key, kvp.Value);
                        }
                } // lock to merge results from this partition into the full set
            //}
            }); // parallel.for over all partitions
        }

        private static string CheckReadForPrimers(int coreIdx, string header, string read,
                                                  List<string>[] primers, int mismatchesAllowedFP, int mismatchesAllowedRP, Dictionary<char, HashSet<char>> allowedBases,
                                                  int fwdPrimerHeadLength, int rvsPrimerHeadLength, int primerCoreLength,
                                                  HashSet<string>[] forwardPrimerFwd, HashSet<string>[] forwardPrimerRC,
                                                  HashSet<string>[] reversePrimerFwd, HashSet<string>[] reversePrimerRC,
                                                  int kMerSize, out string revisedHeader, out string revisedRead, Dictionary<string, int>[] primersUsed)
        {
            string primerType = null;
            revisedHeader = null;
            revisedRead = null;

            //if (header == "@A00632:395:HLMHJDSX7:4:1365:4173:15671 1:N:0:ATGCGCAG+TACTCCTT")
            //    Debugger.Break();
            //if (read == "TTATACAATTTAAAAAAAAATGAGAACTCAAGATTTATTAAGTAATATTAATATATTAAAATTTGATACTATTAGCAAATTAAATATTTATGATTTAAAT")
            //    Debugger.Break();

            int lastPossibleCorePrimerIdx = read.Length - primerCoreLength;

            // (+) ....F.....................(R).....
            // R1  ....F.....................(R).....
            // R2  ....R.....................(F).....
            //
            // (-) ....R.....................(F)..... (after reversal)
            // R1  ....R.....................(F).....
            // R2  ....F.....................(R).....

            // F hhhhcccccccccX, F' Xcccccccccchhhh, R hhhhcccccccccX, R' Xcccccccccchhhh

            for (int i = coreIdx; i <= lastPossibleCorePrimerIdx; i++)
            {
                string possiblePrimerCore = read.Substring(i, primerCoreLength);
                int primerHeadLength = 0;

                bool foundStartingPrimer = false;
                // check for a starting primer (F). Check core first and then head if this matches
                // F hhhhcccccccccX
                if (i >= fwdPrimerHeadLength && forwardPrimerFwd[core].Contains(possiblePrimerCore))
                {
                    string primerHead = read.Substring(i - fwdPrimerHeadLength, fwdPrimerHeadLength);
                    if (forwardPrimerFwd[head].Contains(primerHead))
                    {
                        string primer = primerHead + possiblePrimerCore;
                        //if (primer == "GGTGGTGTTGGATTTACCCAATACGCAACCGC")
                        //    Debugger.Break();
                        if (MismatchesInPrimer(primer, primers[FPpt], allowedBases) <= mismatchesAllowedFP)
                        {
                            // if we found a starting primer (F) in a read, we'll mark the read for pre-trimming
                            primerType = ";FP";
                            revisedHeader = header + primerType;
                            primerHeadLength = fwdPrimerHeadLength;
                            foundStartingPrimer = true;
                            if (primersUsed[FPpt].ContainsKey(primer))
                                primersUsed[FPpt][primer]++;
                            else
                                primersUsed[FPpt].Add(primer, 1);
                        }
                    }
                }

                // check for the other starting primer (R). Check core first and then head if this matches
                // R hhhhcccccccccX
                if (!foundStartingPrimer && i >= rvsPrimerHeadLength && reversePrimerFwd[core].Contains(possiblePrimerCore))
                {
                    string primerHead = read.Substring(i - rvsPrimerHeadLength, rvsPrimerHeadLength);
                    if (reversePrimerFwd[head].Contains(primerHead))
                    {
                        string primer = primerHead + possiblePrimerCore;
                        if (MismatchesInPrimer(primer, primers[RPpt], allowedBases) <= mismatchesAllowedRP)
                        {
                            // if we found a starting primer (R) in a read, we'll mark the read for pre-trimming
                            primerType = ";RP";
                            revisedHeader = header + primerType;
                            primerHeadLength = rvsPrimerHeadLength;
                            foundStartingPrimer = true;
                            if (primersUsed[RPpt].ContainsKey(primer))
                                primersUsed[RPpt][primer]++;
                            else
                                primersUsed[RPpt].Add(primer, 1);
                        }
                    }
                }

                if (foundStartingPrimer)
                {
                    // trim away anything before the start of the primer (and leave the primer)
                    revisedRead = read.Substring(i - primerHeadLength);
                    break;
                }

                bool foundTerminatingPrimer = false;
                int primerLength = 0;
                // check for a terminating primer RC(F). Check core first and then head if this matches
                // F' Xcccccccccchhhh
                if ((i + fwdPrimerHeadLength + primerCoreLength) <= read.Length && forwardPrimerRC[core].Contains(possiblePrimerCore))
                {
                    string primerHead = read.Substring(i + primerCoreLength, fwdPrimerHeadLength);
                    if (forwardPrimerRC[head].Contains(primerHead))
                    {
                        string primer = possiblePrimerCore + primerHead;
                        if (MismatchesInPrimer(primer, primers[rcFPpt], allowedBases) <= mismatchesAllowedFP)
                        {
                            // if we found a starting primer (F') in a read, we'll mark the read 
                            primerType = ";FP'";
                            revisedHeader = header + primerType;
                            foundTerminatingPrimer = true;
                            primerLength = fwdPrimerHeadLength + primerCoreLength;
                            if (primersUsed[rcFPpt].ContainsKey(primer))
                                primersUsed[rcFPpt][primer]++;
                            else
                                primersUsed[rcFPpt].Add(primer, 1);
                        }
                    }
                }

                // check for the other terminating primer RC(R). Check core first and then head if this matches
                // R' Xcccccccccchhhh
                if (!foundTerminatingPrimer && (i + rvsPrimerHeadLength + primerCoreLength) <= read.Length && reversePrimerRC[core].Contains(possiblePrimerCore))
                {
                    string primerHead = read.Substring(i + primerCoreLength, rvsPrimerHeadLength);
                    if (reversePrimerRC[head].Contains(primerHead))
                    {
                        string primer = possiblePrimerCore + primerHead;
                        if (MismatchesInPrimer(primer, primers[rcRPpt], allowedBases) <= mismatchesAllowedRP)
                        {
                            // if we found a terminating primer (R') in a read, we'll mark the read 
                            primerType = ";RP'";
                            revisedHeader = header + primerType;
                            foundTerminatingPrimer = true;
                            primerLength = rvsPrimerHeadLength + primerCoreLength;
                            if (primersUsed[rcRPpt].ContainsKey(primer))
                                primersUsed[rcRPpt][primer]++;
                            else
                                primersUsed[rcRPpt].Add(primer, 1);
                        }
                    }
                }

                if (foundTerminatingPrimer)
                {
                    // trim away anything after the end of the primer (and leave the primer)
                    revisedRead = read.Substring(0, i + primerLength);
                    break;
                }
            }

            return primerType;
        }

        private static int MismatchesInPrimer(string primerFound, List<string> primers, Dictionary<char, HashSet<char>> allowedBases)
        {
            // count mismatches against the specified primers (degenerate bases expanded)
            int minMismatches = int.MaxValue;
            foreach (string primer in primers)
            {
                int mismatches = 0;
                for (int i = 0; i < primerFound.Length; i++)
                {
                    char primerChar = primer[i];
                    char primerFoundChar = primerFound[i];
                    if (!allowedBases[primerChar].Contains(primerFoundChar))
                        mismatches++;
                }
                if (mismatches < minMismatches)
                    minMismatches = mismatches;
            }
            //if (mismatches > 1)
            //    Debugger.Break();
            return minMismatches;
        }

        private static void PartitionReadList(List<string> reads, List<int> partitionStarts, List<int> partitionEnds)
        {
            int partitionsToUse = Environment.ProcessorCount;
            int partitionSize = reads.Count / partitionsToUse;
            if (partitionSize == 0)
                partitionSize = reads.Count;
            int partitionStart = 0;
            int partitionEnd = partitionSize;
            while (partitionEnd < reads.Count)
            {
                partitionStarts.Add(partitionStart);
                partitionEnds.Add(partitionEnd);
                partitionStart += partitionSize;
                partitionEnd += partitionSize;
            }
            partitionStarts.Add(partitionStart);
            partitionEnds.Add(reads.Count);
        }

        private static int PrepStartingReads(List<string> headers, List<string> reads, ScanOpts formWanted, ScanOpts primerLocation, int kMerSize, int primerLength, List<string> preppedReads)
        {
            int readsKept = 0;
            for (int r = 0; r < reads.Count; r++)
            {
                string read = reads[r];

                string primer = primerLocation == ScanOpts.atStart ? read.Substring(0, primerLength) : read.Substring(read.Length - primerLength);

                //if (read == "GTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTACTGGGTTTAAAGGGTGCGTAGGCGGGCGCCTAAGTCGG")
                //    Debugger.Break();

                readsKept++;

                if (formWanted == ScanOpts.RC)
                    read = kMers.ReverseComplement(read);

                preppedReads.Add(read);
            }

            return readsKept;
        }

        // version for normal reads (not dereplicated, so no 'counts')
        private static void CountMersInReads(List<string> reads, int kMerSize, ScanOpts rcWanted, Dictionary<ulong, int> kMerTable)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];

            for (int r = 0; r < reads.Count; r++)
            {
                string read = reads[r];
                // tile for kMers starting anywhere in the read (as long as there's room for at least one kMer)
                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);
                // count all the distinct kMers from the read
                for (int i = 0; i < kMerCount; i++)
                {
                    // add all valid kMers to the filter
                    if (kMersValid[i])
                    {
                        //if (kMersInRead[i] == 0x00E2BBF0BF6B7BC0)
                        //    Debugger.Break();
                        ulong kMer = kMersInRead[i];
                        if (rcWanted == ScanOpts.RC)
                            kMer = kMers.ReverseComplement(kMer, kMerSize);
                        if (kMerTable.ContainsKey(kMer))
                            kMerTable[kMer]++;
                        else
                            kMerTable.Add(kMer, 1);
                    }
                }
            }
        }

        // debugging aid
        private static void ShowDepthsForSeq(string seq, int kMerSize, Dictionary<ulong, int> kMerTable)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];

            // tile for kMers starting anywhere in the read (as long as there's room for at least one kMer)
            int kMerCount = kMers.GenerateMersFromRead(seq, kMerSize, ref kMersInRead, ref kMersValid);
            int[] kMerDepths = new int[kMerCount];
            // count all the distinct kMers from the read
            for (int i = 0; i < kMerCount; i++)
            {
                if (kMersValid[i])
                {
                    kMerDepths[i] = kMerTableLookup(kMersInRead[i], kMerTable, kMerSize);
                }
                else
                    kMerDepths[i] = 0;
            }

            for (int i = 0; i < kMerCount; i++)
            {
                if (i > 0 && i % 10 == 0)
                    Console.WriteLine();
                if (i % 10 == 0)
                    Console.Write(i + "\t");
                Console.Write(kMerDepths[i] + " ");
            }
            Console.WriteLine();
        }

        private static int MeanDepthCount(Dictionary<ulong, int> kMerTable)
        {
            double invSum = 0.0;
            foreach (int depth in kMerTable.Values)
                invSum += 1.0 / (double)depth;

            if (kMerTable.Count == 0)
                return 1;
            else
                return (int)((double)kMerTable.Count / invSum);
        }

        // version that works with dereplicated reads ('counts' gives rep# for associated read)
        private static void CountMersInReads(List<string> reads, List<int> counts, int kMerSize, ScanOpts rcWanted, Dictionary<ulong, int> kMerTable)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];

            for (int r = 0; r < reads.Count; r++)
            {
                string read = reads[r];
                int count = counts[r];
                // tile for kMers starting anywhere in the read (as long as there's room for at least one kMer)
                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);
                // count all the distinct kMers from the read
                for (int i = 0; i < kMerCount; i++)
                {
                    // add all valid kMers to the filter
                    if (kMersValid[i])
                    {
                        //if (kMersInRead[i] == 0x00E2BBF0BF6B7BC0)
                        //    Debugger.Break();
                        ulong kMer = kMersInRead[i];
                        if (rcWanted == ScanOpts.RC)
                            kMer = kMers.ReverseComplement(kMer, kMerSize);
                        if (kMerTable.ContainsKey(kMer))
                            kMerTable[kMer] += count;
                        else
                            kMerTable.Add(kMer, count);
                    }
                }
            }
        }

        private static int InitialiseRegionFilter(List<string> reads, HashSet<ulong> regionFilter, HashSet<ulong>[] filterContextExists, HashSet<ulong>[] filterContexts, int kMerSize, int primerLength, Dictionary<ulong, int> kMerCounts, int minDepth)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int[] kMerDepths = new int[1000];
            int kMersTiled = 0;

            List<string> distinctReads;
            List<int> distinctReadsCounts;
            distinctReads = ReduceToDistinctReads(reads, out distinctReadsCounts);

            for (int r = 0; r < distinctReads.Count; r++)
            {
                string read = distinctReads[r];

                //if (read == "GTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTACTGGGTTTAAAGGGTGCGTAGGCGGGCGCCTAAGTCGGC")
                //    Debugger.Break();

                // tile for kMers starting anywhere in the read (but must have room for a minimum length context)
                int kMerCount = GetDepthsForFilterRead(read, kMerCounts, null, kMerSize, true, ref kMersInRead, ref kMersValid, ref kMerDepths);
                int lastIdx = kMerCount - (shortestContextLength - kMerSize);

                // gather kMers and longest context while the read looks good
                for (int i = 0; i < lastIdx; i++)
                {
                    if (kMerDepths[i] < minDepth / 10)
                        continue;

                    kMersTiled++;
                    ulong kMer = kMersInRead[i];
                    //if (kMer == 0x9778BD941105F6DE)
                    //    Debugger.Break();

                    int lengthIdx;
                    ulong context;
                    bool kMerHasContext = GenerateLongestContext(kMersInRead, kMerCount, kMerSize, i, shortestContextLength, shortestContextLength, out lengthIdx, out context);
                    //if (kMer == 0x4f0fa038503f30ec)
                    //    Console.WriteLine("km@" + i + "(+" + (lengthIdx * contextStride) + ") 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + read.Substring(i, kMerSize + lengthIdx * contextStride));
                    //if (context == 0x96b038991f2ffd21)
                    //    Console.WriteLine("tc@" + i + "(+" + (lengthIdx * contextStride) + ") 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + read.Substring(i, kMerSize + lengthIdx * contextStride));
                    if (kMerHasContext)
                    {
                        regionFilter.Add(kMer);
                        filterContextExists[lengthIdx].Add(kMer);
                        filterContexts[lengthIdx].Add(context);
                    }
                }
            } //  pass through the reads

            return kMersTiled;
        }

        private static int InitialiseEndingFilter(List<string> reads, HashSet<ulong> endingFilter, int kMerSize, int primerLength, Dictionary<ulong, int> kMerCounts, int minDepth)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int[] kMerDepths = new int[1000];
            int kMersTiled = 0;

            List<string> distinctReads;
            List<int> distinctReadsCounts;
            distinctReads = ReduceToDistinctReads(reads, out distinctReadsCounts);

            for (int r = 0; r < distinctReads.Count; r++)
            {
                string read = distinctReads[r];
                //if (read == "TCAGGATGACCAAAAAATCAAAATAAATGTTGATATAAAACTGGATCTCCACCACCAGAAGGATCAAAAAAAGTAGAATTAAAATTTCGATCAAAAAGAAGTATAG")
                //    Debugger.Break();
                read = kMers.ReverseComplement(read);

                // tile for kMers starting anywhere in the read
                int kMerCount = GetDepthsForFilterRead(read, kMerCounts, null, kMerSize, true, ref kMersInRead, ref kMersValid, ref kMerDepths);

                // gather kMers (region) and 48-mer pairs (context) while the read looks good 
                for (int i = 0; i < kMerCount; i++)
                {
                    if (kMerDepths[i] < minDepth / 10)
                        continue;

                    kMersTiled++;
                    ulong kMer = kMersInRead[i];
                    //if (kMer == 0x5d53103440f00353)
                    //    Debugger.Break();

                    int adjacentIdx = i + kMerSize / 2;
                    if (adjacentIdx < kMerCount && kMerDepths[adjacentIdx] >= minDepth)
                    {
                        endingFilter.Add(kMer ^ kMersInRead[adjacentIdx]);
                        //if ((kMer ^ kMersInRead[adjacentIdx]) == 13921541134680753632)
                        //    Debugger.Break();
                    }
                }
            } //  pass through the reads

            return kMersTiled;
        }

        private static void TileReadsForMers(List<string> reads, HashSet<ulong> kMersFromReads, int kMerSize)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            Dictionary<ulong, int> distinctFragments = new Dictionary<ulong, int>();

            int initialMerCount = kMersFromReads.Count;
            kMersFromReads.Clear();

            foreach (string read in reads)
            {
                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);

                //if (read == "GAATAAGTGTTGGTAAAGAATTGGGTCCCCTCCTCCTGCTGGGTCAAAAAATGAGGTGTTTAAGTTTCGGTCTGTTAAAAGTATTGTAATTGCTCCTGCAAGGACAGGAAGAGACAAAAGAAGAA")
                //    Debugger.Break();

                // tile for kMers starting anywhere in the read (as long as there's room for at least one kMer)
                for (int i = 0; i < kMerCount; i++)
                {
                    // skip rest of read if the kMer contains an N (these can block terminating primers and don't want to include them anyway)
                    if (!kMersValid[i])
                        break;

                    // get next kMer tiled from this read
                    ulong kMer = kMersInRead[i];

                    //if (kMer == 0x00E2BBF0BF6B7BC0)
                    //    Debugger.Break();

                    kMersFromReads.Add(kMer);
                }
            }

            Console.WriteLine("retiled " + kMersFromReads.Count + " kMers from " + reads.Count + " reads (was " + initialMerCount + ")");
        }

        private static void ExtendFilterSet(int fip, List<string>[] selectedReads, Dictionary<ulong, int> adapterTrap, HashSet<ulong>[] regionFilter, int kMerSize, bool dropLowComplexity, HashSet<ulong>[] endingFilter, Dictionary<ulong, int>[] kMerCounts,
                                           HashSet<ulong>[] kMersFromSelectedReads, List<int> contextLengths, HashSet<ulong>[][] contextExists, HashSet<ulong>[][] kMerContexts, int[] readsAddingToFilter, int[] endingReadsCount, List<string>[] endingReadsFound)
        {
            List<List<ulong>>[] pairsFromReads = new List<List<ulong>>[2];

            //int droppedReads = 0;
            // sorted distinct selected reads(overlapped reads dropped), with counts of reads subsumed into each distinct read
            List<string>[] distinctReads = new List<string>[2];
            List<int>[] distinctReadsCounts = new List<int>[2];

            // reduce the set of selected reads 
            foreach (int d in directions)
            {
                // reduce the selected reads to distinct, longest ones
                distinctReads[d] = ReduceToDistinctReads(selectedReads[d], out distinctReadsCounts[d]);
                // update kMer counts table with newly selected reads
                CountMersInReads(distinctReads[d], distinctReadsCounts[d], kMerSize, ScanOpts.asIs, kMerCounts[d]);
            }
            int avgDepth = (MeanDepthCount(kMerCounts[fwd]) + MeanDepthCount(kMerCounts[rvs])) / 2;
            int minCount = Math.Max(1, avgDepth / 50);

            // process all the (distinct) selected reads - adding to region filters, contexts, ending filters
            foreach (int d in directions)
            {
                endingReadsFound[d].Clear();
                kMersFromSelectedReads[d].Clear();
                readsAddingToFilter[d] = 0;

                ulong[] kMersInRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int[] kMerDepthsInRead = new int[1000];
                Dictionary<ulong, int> distinctFragments = new Dictionary<ulong, int>(500);
                pairsFromReads[d] = new List<List<ulong>>(distinctReads[d].Count);
                int otherDirection = d == fwd ? rvs : fwd;

                // process every read, adding to region filter + context and (other) ending filter
                for (int r = 0; r < distinctReads[d].Count; r++)
                {
                    string read = distinctReads[d][r];
                    int count = distinctReadsCounts[d][r];
                    bool readAddedToFilter = false;
                    pairsFromReads[d].Add(new List<ulong>());

                    //if (read == "CGGTAATACGTAGGGGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCATTTAAGTCTGGTGTTTAATCCCGGGGCTCAACCCCGGATCGCACTGGAAACTGGGTGACTTGAGTGC")
                    //    Debugger.Break();

                    int kMerCount = GetDepthsForFilterRead(read, kMerCounts[d], kMerCounts[otherDirection], kMerSize, false, ref kMersInRead, ref kMersValid, ref kMerDepthsInRead);
                    int lastIdx = kMerCount - (shortestOtherContextLength - kMerSize);

                    // tile for kMers starting anywhere in the read
                    //for (int i = 0; i < lastIdx; i++)
                    for (int i = 0; i < kMerCount; i++)
                    {
                        // get next kMer tiled from this read
                        ulong kMer = kMersInRead[i];
                        //bool targetkMer = false;

                        //if (kMer == 0x50af5213b42acac2)
                        //{
                        //    targetkMer = true;
                        //    Console.WriteLine("found kMer " + "fip=" + fip + " dir=" + d + " r=" + r + " @" + i);
                        //    //Debugger.Break();
                        //}

                        // skip poor kMers in the middle of a read (GetDepthsForFilterRead will have tail-trimmed)
                        if (kMerDepthsInRead[i] < minCount && i < lastIdx * 3 / 4)
                            continue;

                        // stop tiling if we hit an adapter region
                        if (adapterTrap.ContainsKey(kMer))
                            break;

                        // if we haven't see it before, remember to add it to the filter set later (unless we find we've finished extending already)
                        if (!regionFilter[d].Contains(kMer))
                        {
                            kMersFromSelectedReads[d].Add(kMer);
                            //if (targetkMer)
                            //    Console.WriteLine("added kMer " + "fip=" + fip + " dir=" + d + " r=" + r + " @" + i);

                            // first kMer from this read...
                            if (!readAddedToFilter)
                            {
                                readAddedToFilter = true;
                                readsAddingToFilter[d] += count;
                            }
                        }

                        // hashed 48-mer from (including) current kMer
                        int adjacentIdx = i + kMerSize / 2;
                        if (adjacentIdx < kMerCount && kMerDepthsInRead[adjacentIdx] >= minCount)
                        {
                            ulong adjacentMer = kMersInRead[adjacentIdx];
                            ulong endingPair = kMers.ReverseComplement(kMer, kMerSize) ^ kMers.ReverseComplement(adjacentMer, kMerSize);
                            endingFilter[otherDirection].Add(endingPair);
                            //if (endingPair == 13921541134680753632)
                            //    Debugger.Break();
                            pairsFromReads[d][r].Add(kMer ^ adjacentMer);
                            //if ((kMer ^ adjacentMer) == 13921541134680753632)
                            //    Debugger.Break();
                        }

                        int contextLengthIdx;
                        ulong context;
                        bool kMerHasContext = GenerateLongestContext(kMersInRead, kMerCount, kMerSize, i, shortestOtherContextLength, shortestContextLength, out contextLengthIdx, out context);
                        //if (kMer == 0x4f0fa038503f30ec)
                        //    Console.WriteLine(r + ":km@" + i + "(+" + (lengthIdx * contextStride) + ") 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + read.Substring(i, kMerSize + lengthIdx * contextStride));
                        //if (context == 0x96b038991f2ffd21)
                        //    Console.WriteLine(r + ":tc@" + i + "(+" + (lengthIdx * contextStride) + ") 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + read.Substring(i, kMerSize + lengthIdx * contextStride));
                        if (kMerHasContext)
                        {
                            contextExists[d][contextLengthIdx].Add(kMer);
                            kMerContexts[d][contextLengthIdx].Add(context);
                            RemoveShorterContexts(kMersInRead, kMerSize, i, contextLengthIdx, contextLengths, contextExists[d], kMerContexts[d]);
                        }
                    }
                } // for each selected read
            } // directions

            // finally see how many of these selected reads have run into their respective ending clouds
            foreach (int d in directions)
            {
                endingReadsCount[d] = 0;

                for (int r = 0; r < distinctReads[d].Count; r++)
                {
                    List<ulong> pairsFromRead = pairsFromReads[d][r];
                    int count = distinctReadsCounts[d][r];
                    int matchingPairs = 0;
                    int firstMatch = -1;
                    for (int p = 0; p < pairsFromRead.Count; p++)
                    {
                        ulong pair = pairsFromRead[p];
                        if (endingFilter[d].Contains(pair))
                        {
                            matchingPairs++;
                            if (firstMatch < 0)
                                firstMatch = p;
                        }
                    }
                    // must match most of the remainder of the read (after the first match)
                    if (matchingPairs > 0 && matchingPairs > 3 * (pairsFromRead.Count - firstMatch) / 4)
                    {
                        endingReadsFound[d].Add(distinctReads[d][r]);
                        endingReadsCount[d] += count;
                    }
                }
            }

            //int kMersAddedCount = kMersFromSelectedReads[fwd].Count + kMersFromSelectedReads[rvs].Count;
            //int readsAddingToFilterCount = readsAddingToFilter[fwd] + readsAddingToFilter[rvs];
            //int selectedReadsCount = selectedReads[fwd].Count + selectedReads[rvs].Count;
            //Console.WriteLine("kept " + selectedReadsCount + " (" + selectedReads[fwd].Count + "+" + selectedReads[rvs].Count + ")" + " matching reads. " +
            //                   readsAddingToFilterCount + " (" + readsAddingToFilter[fwd] + "+" + readsAddingToFilter[rvs] + ")" + " reads added " +
            //                   kMersAddedCount + " (" + kMersFromSelectedReads[fwd].Count + "+" + kMersFromSelectedReads[rvs].Count + ")" + " kMers to region filter. " +
            //                   endingReadsCount + " (" + endingReadsFound[fwd].Count + "+" + endingReadsFound[rvs].Count + ")" + " ending reads found. ");
        }

        // generate the longest possible context for a given kMer (idx)
        private static bool GenerateLongestContext(ulong[] kMersInRead, int kMerCount, int kMerSize, int kMerIdx, int smallestLength, int minLength, out int lengthIdx, out ulong context)
        {
            int kMersLeft = kMerCount - kMerIdx;                                          // kMers in the remainder of the read
            int ignoredBases = kMersLeft - (kMersLeft / filterContextStride * filterContextStride);   // how many bases are ignored at the end of the read due to striding
            int finalIdx = kMerCount - 1 - ignoredBases;                                  // index of last kMer in the context
            int contextLength = (finalIdx - kMerIdx) + kMerSize;                          // how long would this context be?

            if (contextLength < minLength)
            {
                lengthIdx = 0;
                context = 0;
                return false;
            }

            lengthIdx = (contextLength - smallestLength) / filterContextStride;

            context = kMersInRead[kMerIdx] ^ kMersInRead[finalIdx];
            for (int m = kMerIdx + kMerSize; m < finalIdx; m += kMerSize)   // add any in-between kMers to the hash
                context ^= kMersInRead[m];

            return true;
        }

        // generate a kMer context of the specified length
        private static bool GenerateWantedContext(ulong[] kMersInRead, bool[] kMersValid, int kMerCount, int kMerSize, int kMerIdx, int wantedLength, out ulong context)
        {
            context = 0;

            int finalIdx = kMerIdx + wantedLength - kMerSize;
            if (finalIdx >= kMerCount)
                return false;
            for (int m = kMerIdx; m > finalIdx; m++)
                if (!kMersValid[m])
                    return false;
            context = kMersInRead[kMerIdx] ^ kMersInRead[finalIdx];
            for (int m = kMerIdx + kMerSize; m < finalIdx; m += kMerSize)   // add any in-between kMers to the hash
                context ^= kMersInRead[m];

            return true;
        }

        // remove any shorter contexts that have been made redundant by the discovery of a longer context
        private static void RemoveShorterContexts(ulong[] kMersInRead, int kMerSize, int kMerIdx, int contextLengthIdx, List<int> contextLengths, HashSet<ulong>[] contextExists, HashSet<ulong>[] kMerContexts/*, int r*/)
        {
            ulong kMer = kMersInRead[kMerIdx];

            for (int li = 0; li < contextLengthIdx; li++)
            {
                if (contextExists[li].Contains(kMer))
                {
                    // index of next (last) kMer to be folded
                    int finalIdx = kMerIdx + contextLengths[li] - kMerSize;
                    ulong context = kMersInRead[kMerIdx] ^ kMersInRead[finalIdx];
                    for (int m = kMerIdx + kMerSize; m < finalIdx; m += kMerSize)   // add any in-between kMers to the hash
                        context ^= kMersInRead[m];
                    //bool removed = false;
                    if (kMerContexts[li].Contains(context))
                    {
                        kMerContexts[li].Remove(context);
                        //removed = true;
                    }
                    //if (removed && kMer == 0x4f0fa038503f30ec)
                    //    Console.WriteLine(r + ":rm-k #(+" + (li*contextStride) +") 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + kMers.ExpandMer(kMersInRead[kMerIdx], kMerSize) + "^" + kMers.ExpandMer(kMersInRead[finalIdx], kMerSize));
                    //if (removed && context == 0x96b038991f2ffd21)
                    //    Console.WriteLine(r + ":rm-c #(+" + (li*contextStride) + " 0x" + kMer.ToString("x16") + "-->0x" + context.ToString("x16") + " " + kMers.ExpandMer(kMersInRead[kMerIdx], kMerSize) + "^" + kMers.ExpandMer(kMersInRead[finalIdx], kMerSize));
                }
            }
        }

        // Mapping low 5 bits of ASCII chars to ACGT or not. Used to detect ambiguous bases. Summing over all translated bases will produce 0 if there are no ambiguous bases
        static readonly int[] baseToNotACGT = new int[] { 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        //                                                @, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, [, \, ], ^, _

        public static bool SeqContainsAmbiguousBases(Sequence seq)
        {
            int nonACGTCount = 0;
            for (int i = 0; i < seq.Length; i++)
                nonACGTCount += baseToNotACGT[seq.Bases[i] & 0x1f];

            return nonACGTCount > 0;
        }

        private static bool CheckReadForLowComplexity(string read, ref ulong[] tripletSet, ref bool[] tripletValid, Dictionary<ulong, int> distinctTriplets)
        {
            distinctTriplets.Clear();

            int tripletsCount = kMers.GenerateMersFromRead(read, 3, ref tripletSet, ref tripletValid);
            int startIdx = (read.Length - 3 + 1) / 2;       // check second half of the read only
            int tripletsConsidered = tripletsCount - startIdx;

            for (int i = startIdx; i < tripletsCount; i++)
            {
                if (tripletValid[i])
                {
                    ulong triplet = tripletSet[i];
                    if (distinctTriplets.ContainsKey(triplet))
                        distinctTriplets[triplet]++;
                    else
                        distinctTriplets.Add(triplet, 1);
                }
            }

            int maxCount = 0;
            foreach (int count in distinctTriplets.Values)
                if (count >= maxCount)
                    maxCount = count;
            int max2ndCount = 0;
            foreach (int count in distinctTriplets.Values)
                if (count < maxCount && count > max2ndCount)
                    max2ndCount = count;
            if (max2ndCount < maxCount / 5)
                max2ndCount = 0;
            int avgCount = tripletsConsidered / distinctTriplets.Count;

            int closeToMaxCount = 0;
            int closeToMaxSum = 0;
            foreach (int count in distinctTriplets.Values)
                if (VeryClose(count, maxCount) || (max2ndCount > 0 && VeryClose(count, max2ndCount)))
                {
                    closeToMaxCount++;
                    closeToMaxSum += count;
                }
            int avgMax = closeToMaxSum / closeToMaxCount;

            // no significant repetition found
            if (closeToMaxCount > tripletsConsidered / 4)
                return false;

            // low average count and mostly the same counts for the triplets
            if (avgCount <= 4 && Close(avgCount, avgMax))
                return false;

            // high rep triplets account for most
            return closeToMaxSum > tripletsConsidered * 2 / 3;
        }

        // filter (filtered) in-memory reads against the agreed, final set of kMers
        private static void FinalFilterReads(Dictionary<int, uint>[][] primerReads, List<int>[] readsToFilterCount, List<string>[] headersToFilter, List<string>[] readsToFilter,
                                             int kMerSize, HashSet<ulong> regionFilter, List<string> selectedHeaders, List<string> selectedReads, int longestReadLength, bool pairedReads, Dictionary<int, int> readPairs)
        {
            ulong[] kMersInRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            //bool chosen = false;
            int fileStride = readsToFilterCount.Length;

            // remember the reads we've chosen to retain from both sets of reads
            Dictionary<int, int>[] retained = null;
            if (pairedReads)
            {
                retained = new Dictionary<int, int>[fileStride];
                for (int fip = 0; fip < fileStride; fip++)
                    retained[fip] = new Dictionary<int, int>(readsToFilter[fip].Count);
            }

            // primerReads is also used for partitioned WGS reads. The first index is fileStride, the second is partition# which always zero here as the in-memory reads are not partitioned

            for (int fip = 0; fip < fileStride; fip++)
                for (int r = 0; r < readsToFilterCount[fip][0]; r++)
                {
                    string header = headersToFilter[fip][r];
                    string read = readsToFilter[fip][r];

                    //if (header.StartsWith("@A00632:395:HLMHJDSX7:4:1365:4173:15671"))
                    //    Debugger.Break();
                    //if (read == "GCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCTAGCGTTGTCCGGAATCATTGGGCGTAAAGCGCGTGTAGGCGGCCCGGTAAGTCCGCTGTGAAAGTCCAGGGCTCAACCCTGGGATGCCGGTGGATACTGTCGGGCTCGAG")
                    //{
                    //    chosen = true;
                    //    Debugger.Break();
                    //}

                    // don't add the primer reads twice - and we have already selected the trimmed versions of these reads
                    if (primerReads[fip][0].ContainsKey(r))
                    {
                        if (retained != null)
                            retained[fip].Add(r, (int)(primerReads[fip][0][r]));
                        continue;
                    }

                    int kMersCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);
                    //if (chosen)
                    //{
                    //    bool[] kMerPresent = new bool[1000];
                    //    for (int i = 0; i < kMersCount; i++)
                    //        kMerPresent[i] = regionFilter.Contains(kMersInRead[i]);
                    //    Debugger.Break();
                    //}

                    // look for matches across all parts of the read
                    int firstThirdEndIdx = kMersCount / 3;
                    int secondThirdEndIdx = firstThirdEndIdx * 2;
                    int thirdThirdEndIdx = firstThirdEndIdx * 3;

                    bool foundInFirstThird = false;
                    for (int i = 0; i < firstThirdEndIdx; i++)
                    {
                        if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                        {
                            foundInFirstThird = true;
                            break;
                        }
                    }

                    bool foundInSecondThird = false;
                    for (int i = firstThirdEndIdx; i < secondThirdEndIdx; i++)
                    {
                        if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                        {
                            //if (header.StartsWith("@RM2|S1|R13243951/1"))
                            //    Debugger.Break();
                            foundInSecondThird = true;
                            break;
                        }
                    }

                    bool foundInThirdThird = false;
                    for (int i = secondThirdEndIdx; i < thirdThirdEndIdx; i++)
                    {
                        if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                        {
                            //if (header.StartsWith("@RM2|S1|R13243951/1"))
                            //    Debugger.Break();
                            foundInThirdThird = true;
                            break;
                        }
                    }

                    if ((foundInFirstThird && foundInSecondThird) || (foundInSecondThird && foundInThirdThird) || (foundInFirstThird && foundInThirdThird))
                    {
                        // save the selected read 
                        selectedHeaders.Add(header);
                        selectedReads.Add(read);

                        // remember we've retained this read and its index in the selectedReads list
                        if (retained != null)
                            retained[fip].Add(r, selectedReads.Count - 1);
                    }
                }

            // build the paired reads data structures
            if (retained != null)
            {
                // find the readsToFilter indexes that were retained from both R1 and R2
                HashSet<int>[] retainedSet = new HashSet<int>[fileStride];
                for (int fip = 0; fip < fileStride; fip++)
                    retainedSet[fip] = new HashSet<int>(retained[fip].Keys);
                retainedSet[0].IntersectWith(retainedSet[1]);

                int trivialPairsFound = 0;

                // build the pair index into selectedReads (R1 to R2 and R2 to R1)
                foreach (int rtfIndex in retainedSet[0])
                {
                    int R1Idx = retained[0][rtfIndex];
                    int R2Idx = retained[1][rtfIndex];

                    // canonical order
                    if (R2Idx < R1Idx)
                    {
                        int temp = R1Idx;
                        R1Idx = R2Idx;
                        R2Idx = temp;
                    }

                    if (TrivialPair(selectedReads[R1Idx], selectedReads[R2Idx], longestReadLength))
                        trivialPairsFound++;
                    else
                        readPairs.Add(R1Idx, R2Idx);
                }

                Console.WriteLine(trivialPairsFound + " trivial paired reads found and ignored");
            }
        }

        // filter (unfiltered) WGS reads against the agreed, final set of kMers
        private static void FinalFilterReads(Dictionary<int, uint>[][] primerReads, List<int>[] readsToFilterCount, List<string>[] WGSReadsFNs,
                                             int kMerSize, HashSet<ulong> regionFilter, List<string> selectedHeaders, List<string> selectedReads, int longestReadLength, bool pairedReads, Dictionary<int, int> readPairs, ParallelOptions threads)
        {
            int fileStride = WGSReadsFNs.Length;
            int trivialPairsFound = 0;

            Parallel.For(0, WGSReadsFNs[0].Count, threads, p =>
            //for (int p = 0; p < WGSReadsFNs[0].Count; p++)
            {
                StreamReader[] WGSReads = new StreamReader[fileStride];
                Dictionary<int, uint>[] primerReadsForPartition = new Dictionary<int, uint>[fileStride];
                int readsInPartition = int.MaxValue;
                for (int fip = 0; fip < fileStride; fip++)
                {
                    WGSReads[fip] = new StreamReader(WGSReadsFNs[fip][(int)p]);
                    primerReadsForPartition[fip] = primerReads[fip][p];
                    // partitions from each file should be the same size but play safe here
                    if (readsToFilterCount[fip][(int)p] < readsInPartition)
                        readsInPartition = readsToFilterCount[fip][(int)p];
                }

                List<string> localMatchingHeaders = new List<string>(1000);
                List<string> localMatchingReads = new List<string>(1000);

                ulong[] kMersInRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];

                // remember the reads we've chosen to retain from both sets of reads (per partition)
                Dictionary<int, int>[] retained = null;                 // value is into localMatchingReads
                Dictionary<int, int>[] retainedForPrimers = null;       // value is into selectedReads
                if (pairedReads)
                {
                    retained = new Dictionary<int, int>[fileStride];
                    retainedForPrimers = new Dictionary<int, int>[fileStride];
                    for (int fip = 0; fip < fileStride; fip++)
                    {
                        retained[fip] = new Dictionary<int, int>(readsInPartition);
                        retainedForPrimers[fip] = new Dictionary<int, int>(readsInPartition);
                    }
                }

                for (int fip = 0; fip < fileStride; fip++)
                {
                    for (int r = 0; r < readsInPartition; r++)
                    {
                        string header = WGSReads[fip].ReadLine();
                        string read = WGSReads[fip].ReadLine();

                        // perhaps it was zapped because it contained Ns
                        if (read.Length == 0)
                            continue;

                        // don't add the primer reads twice - and we have already selected the trimmed versions of these reads
                        if (primerReadsForPartition[fip].ContainsKey(r))
                        {
                            if (retained != null)
                                // these are indexes into the global selectedReads, not partition-relative.
                                retainedForPrimers[fip].Add(r, (int)(primerReads[fip][p][r]));
                            continue;
                        }

                        //if (header[fip].StartsWith("@D00534:185:CBU0KANXX:6:2211:8061:7774"))
                        //    Debugger.Break();
                        //if (read == "TAATATAAGATTTTGACTTCTTCCCCCGTCGTTATCTTTATTATTATCTAGATCAATTATTGAAAATGGAGCAGGAACTGGATGAACAGTTTACCCCCCTCTATCTTCTGGAATTGCTCATGCAG")
                        //    Debugger.Break();

                        // look for matches along the read
                        int kMersCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);
                        // look for matches across all parts of the read
                        int firstThirdEndIdx = kMersCount / 3;
                        int secondThirdEndIdx = firstThirdEndIdx * 2;
                        int thirdThirdEndIdx = firstThirdEndIdx * 3;

                        bool foundInFirstThird = false;
                        for (int i = 0; i < firstThirdEndIdx; i++)
                        {
                            if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                            {
                                foundInFirstThird = true;
                                break;
                            }
                        }

                        bool foundInSecondThird = false;
                        for (int i = firstThirdEndIdx; i < secondThirdEndIdx; i++)
                        {
                            if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                            {
                                //if (header.StartsWith("@RM2|S1|R13243951/1"))
                                //    Debugger.Break();
                                foundInSecondThird = true;
                                break;
                            }
                        }

                        bool foundInThirdThird = false;
                        for (int i = secondThirdEndIdx; i < thirdThirdEndIdx; i++)
                        {
                            if (kMersValid[i] && regionFilter.Contains(kMersInRead[i]))
                            {
                                //if (header.StartsWith("@RM2|S1|R13243951/1"))
                                //    Debugger.Break();
                                foundInThirdThird = true;
                                break;
                            }
                        }

                        if ((foundInFirstThird && foundInSecondThird) || (foundInSecondThird && foundInThirdThird))
                        {
                            localMatchingHeaders.Add(header);
                            localMatchingReads.Add(read);

                            // remember we've retained this read and its index in the local selectedReads list
                            if (retained != null)
                                retained[fip].Add(r, localMatchingReads.Count - 1);
                        }

                    } // for each read in the file
                } // for each file in the pair

                for (int fip = 0; fip < fileStride; fip++)
                    WGSReads[fip].Close();

                // finally copy any selected reads from this partition into the global 'selected' sets
                lock (selectedHeaders)
                {
                    // remember where we started copying from local to global
                    int selectedOffset = selectedReads.Count;

                    // copy the selected reads across
                    for (int i = 0; i < localMatchingReads.Count; i++)
                    {
                        selectedHeaders.Add(localMatchingHeaders[i]);
                        selectedReads.Add(localMatchingReads[i]);
                    }

                    if (retained != null)
                    {
                        /// convert per-partition list indices into global ones (for compatibility with the primer-selected retained entries)
                        for (int fip = 0; fip < fileStride; fip++)
                        {
                            foreach (int r in retained[fip].Keys)
                                retained[fip][r] += selectedOffset;
                            // and add in the (now compatible) primer-selected reads
                            foreach (KeyValuePair<int, int> kvp in retainedForPrimers[fip])
                                retained[fip].Add(kvp.Key, kvp.Value);
                        }

                        // find any reads that were kept from both R1 and R2 
                        HashSet<int>[] retainedSet = new HashSet<int>[fileStride];
                        if (retained != null)
                        {
                            for (int fip = 0; fip < fileStride; fip++)
                                retainedSet[fip] = new HashSet<int>(retained[fip].Keys);
                            retainedSet[0].IntersectWith(retainedSet[1]);
                        }

                        // build the pair index into selectedReads (R1 to R2 and R2 to R1)
                        foreach (int rtfIndex in retainedSet[0])
                        {
                            int R1Idx = retained[0][rtfIndex];
                            int R2Idx = retained[1][rtfIndex];


                            // canonical order
                            if (R2Idx < R1Idx)
                            {
                                int temp = R1Idx;
                                R1Idx = R2Idx;
                                R2Idx = temp;
                            }

                            if (TrivialPair(selectedReads[R1Idx], selectedReads[R2Idx], longestReadLength))
                                trivialPairsFound++;
                            else
                                readPairs.Add(R1Idx, R2Idx);
                        }
                    }
                }
                //}
            });

            Console.WriteLine(trivialPairsFound + " trivial paired reads found and ignored");
        }

        // find (filtered) in-memory reads that match against the (newly-updated) filter kMers, these reads will then be tiled for more kMers to add to the filter
        private static void FindMatchingReads(int iteration, List<string> readsToFilter, HashSet<int> readsAlreadyProcessed, int kMerSize, bool dropLowComplexity,
                                              HashSet<ulong>[] regionFilter, List<int> contextLengths, HashSet<ulong>[][] contextExists, HashSet<ulong>[][] kMerContexts, List<string>[] matchingReads,
                                              out int droppedLowComplexityReads, ParallelOptions threads)
        {
            foreach (int d in directions)
                matchingReads[d].Clear();

            droppedLowComplexityReads = 0;
            // can't just add to out parameters from within the parallel.for loop
            int localDroppedLowComplexityReads = 0;

            List<int> partitionStarts = new List<int>();
            List<int> partitionEnds = new List<int>();
            PartitionReadList(readsToFilter, partitionStarts, partitionEnds);
            HashSet<int>[] partitionReadsAlreadyProcessed = new HashSet<int>[partitionStarts.Count];

            Parallel.For(0, partitionStarts.Count, threads, p =>
            //for (int p = 0; p < partitionStarts.Count; p++)
            {
                List<string>[] localMatchingReads = new List<string>[2];
                foreach (int d in directions)
                    localMatchingReads[d] = new List<string>(matchingReads[d].Capacity);
                partitionReadsAlreadyProcessed[p] = new HashSet<int>();
                int threadDroppedLowComplexityReads = 0;

                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                Dictionary<ulong, int> distinctTriplets = new Dictionary<ulong, int>(500);

                for (int r = partitionStarts[p]; r < partitionEnds[p]; r++)
                {
                    // skip over reads that have already been processed in previous passes
                    if (readsAlreadyProcessed.Contains(r))
                        continue;

                    string read = readsToFilter[r];

                    bool targetRead = false;
                    //if (read == "ACGGAGGATCCGAGCGTTATCCGGATTTACTGGGTTTAAAGGGTGCGTAGGCGGGCGCCTAAGTCGGCGGTGAAAGGCAGCGGCTCAACCGTTGAGCTGCC")
                    //{
                    //    targetRead = true;
                    //    Debugger.Break();
                    //}

                    if (read.Length < kMerSize * 2)
                        continue;

                    // look for a match at the start of a read (F-->)
                    ulong startingMer;
                    bool kMerOK = kMers.CondenseMer(read, 0, kMerSize, out startingMer);
                    if (!kMerOK)
                        continue;
                    // or at the end of an reverse direction read (-->F') if this fails
                    bool matchAtStart = regionFilter[fwd].Contains(startingMer) || regionFilter[rvs].Contains(startingMer);
                    bool matchAtEnd = false;
                    if (!matchAtStart)
                    {
                        ulong endingMerRC;
                        kMerOK = kMers.CondenseMer(read, (read.Length - kMerSize), kMerSize, out endingMerRC);
                        if (!kMerOK)
                            continue;
                        endingMerRC = kMers.ReverseComplement(endingMerRC, kMerSize);

                        matchAtEnd = regionFilter[fwd].Contains(endingMerRC) || regionFilter[rvs].Contains(endingMerRC);
                        if (matchAtEnd)
                        {
                            read = kMers.ReverseComplement(read);
                            startingMer = endingMerRC;
                        }
                    }

                    //if (startingMer == 0xE3020FA34517520A)
                    //    Debugger.Break();

                    if (matchAtStart || matchAtEnd)
                    {
                        //if (targetRead)
                        //{
                        //    Console.WriteLine(iteration + ": matched" + (matchAtStart ? " start" : "") + (matchAtEnd ? " end" : ""));
                        //    //Debugger.Break();
                        //}

                        int d = regionFilter[fwd].Contains(startingMer) ? fwd : rvs;
                        int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);

                        // does the kMer context also match? Try for the longest context first
                        bool readMatches = false;
                        for (int li = kMerContexts[d].Length - 1; li >= 0; li--)
                        {
                            if (contextExists[d][li].Contains(startingMer))
                            {
                                ulong kMerContext;
                                bool contextPossible = GenerateWantedContext(kMersFromRead, kMersValid, kMerCount, kMerSize, 0, contextLengths[li], out kMerContext);
                                if (contextPossible && kMerContexts[d][li].Contains(kMerContext))
                                {
                                    readMatches = true;
                                    break;
                                }
                            }
                        }

                        if (readMatches)
                        {
                            //if (readsToFilter[r] == "ATACGAGCCCCTGGAATATCCTTTGATAAAATACCTCTATTTGTTTGATCAGTTTTAATTACAGCAATTTTATTATTATTGTCTTCACCTGTATTAGCAGGAGCAATTACTATATTATTAACAGA")
                            //    Debugger.Break();

                            // save this matching read for tiling and adding to the filters - unless it appears to be tainted by a low-complexity region
                            if (dropLowComplexity && CheckReadForLowComplexity(read, ref kMersFromRead, ref kMersValid, distinctTriplets))
                                threadDroppedLowComplexityReads++;
                            else
                                // save the selected read
                                localMatchingReads[d].Add(read);

                            // and remember not to process it again 
                            partitionReadsAlreadyProcessed[p].Add(r);
                        }
                    }
                    else
                    {
                        if (targetRead)
                        {
                            Console.WriteLine(iteration + ": no match");
                        }
                    }
                }

                lock (matchingReads)
                {
                    foreach (int d in directions)
                        foreach (string read in localMatchingReads[d])
                            matchingReads[d].Add(read);

                    localDroppedLowComplexityReads += threadDroppedLowComplexityReads;
                }
            });
            //}

            for (int p = 0; p < partitionStarts.Count; p++)
                foreach (int r in partitionReadsAlreadyProcessed[p])
                    readsAlreadyProcessed.Add(r);

            droppedLowComplexityReads = localDroppedLowComplexityReads;
        }

        // find (unfiltered) WGS reads that match against the (newly-updated) filter kMers, these reads will then be tiled for more kMers to add to the filter
        private static void FindMatchingReads(int iteration, List<int> readsInPartitions, List<string> WGSReadsFNs, List<HashSet<int>> readsAlreadyProcessed, int kMerSize, bool dropLowComplexity,
                                              HashSet<ulong>[] regionFilter, List<int> contextLengths, HashSet<ulong>[][] contextExists, HashSet<ulong>[][] kMerContexts, List<string>[] matchingReads,
                                              out int droppedLowComplexityReads, ParallelOptions threads)
        {
            foreach (int d in directions)
                matchingReads[d].Clear();

            droppedLowComplexityReads = 0;
            int localDroppedLowComplexityReads = 0;     // can't add to out param from within parallel.for block

            Parallel.For(0, WGSReadsFNs.Count, threads, p =>
            //for (int p = 0; p < readsInPartitions.Count; p++)
            {
                List<string>[] threadMatchingReads = new List<string>[2];
                foreach (int d in directions)
                    threadMatchingReads[d] = new List<string>(matchingReads[d].Capacity);
                int threadDroppedLowComplexityReads = 0;

                StreamReader WGSReads = new StreamReader(WGSReadsFNs[p]);
                HashSet<int> alreadyProcessedInPartition = readsAlreadyProcessed[p];
                //bool wantedFile = WGSReadsFNs[p].EndsWith("IRE_node4_plus_3_0.tmp");

                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                Dictionary<ulong, int> distinctTriplets = new Dictionary<ulong, int>(500);

                for (int r = 0; r < readsInPartitions[p]; r++)
                {
                    // look for a match at the start of a read
                    string header = WGSReads.ReadLine();
                    string read = WGSReads.ReadLine();
                    //bool wantedRead = false;

                    //if (header == ">SRR13051709.39346134 39346134 length=151")
                    //{
                    //    wantedRead = true;
                    //    Debugger.Break();
                    //}
                    //if (read == "GATCGTATCCGAAGAATCCGAGTCTGCCGTGTGCCTCTTTGTGGAGGTACATGGAGAGATACCATGCAGAGAGACCAGCGTTGCTGTGTCCTGTTGCAAGTGCTGTTGCACCGCCTGCTGCAATTGAGATAACCGTTGCTCTCTGTGAACC")
                    //{
                    //    wantedRead = true;
                    //    Console.WriteLine(p + "@" + r);
                    //}

                    // skip over reads that have already been processed
                    if (alreadyProcessedInPartition.Contains(r))
                        continue;

                    // reads must be at least this long to be matched (allows room for checking more than just the leading kMers)
                    if (read.Length < kMerSize * 2)
                        continue;

                    // look for a match at the start of a read (F-->)
                    ulong startingMer;
                    bool kMerOK = kMers.CondenseMer(read, 0, kMerSize, out startingMer);
                    if (!kMerOK)
                        continue;                    // or at the end of an reverse direction read (-->F') if this fails
                    bool matchAtStart = regionFilter[fwd].Contains(startingMer) || regionFilter[rvs].Contains(startingMer);
                    bool matchAtEnd = false;
                    if (!matchAtStart)
                    {
                        ulong endingMerRC;
                        kMerOK = kMers.CondenseMer(read, (read.Length - kMerSize), kMerSize, out endingMerRC);
                        if (!kMerOK)
                            continue;
                        endingMerRC = kMers.ReverseComplement(endingMerRC, kMerSize);

                        matchAtEnd = regionFilter[fwd].Contains(endingMerRC) || regionFilter[rvs].Contains(endingMerRC);
                        if (matchAtEnd)
                        {
                            read = kMers.ReverseComplement(read);
                            startingMer = endingMerRC;
                        }
                    }

                    //if (read == "AGAAGTATTTAAATTTCGATCTGTTAGTAATATTGTAATTGCTCCTGCTAATACCCCATTATTACTTTTTGTGTAATCTATTTTTTCTCAAAAAGTCAATTCTGACTTGGTGTGTTTTATTACAT")
                    //    Debugger.Break();
                    //if (wantedRead)
                    //    Console.WriteLine("wanted read: " + matchAtStart + ", " + matchAtEnd);

                    if (matchAtStart || matchAtEnd)
                    {
                        //if (wantedRead && p == 7 && r == 1552106)
                        //    Debugger.Break();                        
                        int d = regionFilter[fwd].Contains(startingMer) ? fwd : rvs;
                        int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);

                        //if (wantedRead)
                        //    Debugger.Break();

                        // does the kMer context also match? Try for the longest context first
                        bool readMatches = false;
                        for (int li = kMerContexts[d].Length - 1; li >= 0; li--)
                        {
                            if (contextExists[d][li].Contains(startingMer))
                            {
                                ulong kMerContext;
                                bool contextPossible = GenerateWantedContext(kMersFromRead, kMersValid, kMerCount, kMerSize, 0, contextLengths[li], out kMerContext);

                                if (contextPossible && kMerContexts[d][li].Contains(kMerContext))
                                {
                                    readMatches = true;
                                    break;
                                }
                            }
                        }

                        if (readMatches)
                        {
                            //if (readsToFilter[r] == "ATACGAGCCCCTGGAATATCCTTTGATAAAATACCTCTATTTGTTTGATCAGTTTTAATTACAGCAATTTTATTATTATTGTCTTCACCTGTATTAGCAGGAGCAATTACTATATTATTAACAGA")
                            //    Debugger.Break();

                            // save this matching read for tiling and adding to the filters - unless it appears to be tainted by a low-complexity region
                            if (dropLowComplexity && CheckReadForLowComplexity(read, ref kMersFromRead, ref kMersValid, distinctTriplets))
                                localDroppedLowComplexityReads++;
                            else
                                threadMatchingReads[d].Add(read);

                            // and remember not to process it again 
                            alreadyProcessedInPartition.Add(r);
                        }
                    }
                }

                if (WGSReads != null)
                    WGSReads.Close();

                lock (matchingReads)
                {
                    foreach (int d in directions)
                        foreach (string read in threadMatchingReads[d])
                        {
                            matchingReads[d].Add(read);
                        }

                    localDroppedLowComplexityReads += threadDroppedLowComplexityReads;
                }
                //}
            });

            droppedLowComplexityReads = localDroppedLowComplexityReads;
        }

        private static int WriteFilteredReads(StreamWriter filteredReads, int firstReadToWrite, List<string> selectedHeaders, List<string> selectedReads, int minLength, string tag)
        {
            int readsWritten = 0;

            for (int i = firstReadToWrite; i < selectedHeaders.Count; i++)
            {
                string read = selectedReads[i];
                string header = selectedHeaders[i];

                //if (header.StartsWith("A00632:395:HLMHJDSX7:4:1623:28881:8578 1:N:0:ATGCGCAG+TACTCCTT"))
                //    Debugger.Break();

                if (header[0] == '@')
                    header = ">" + header.Substring(1);

                if (read == null)
                {
                    filteredReads.WriteLine(header);
                    filteredReads.WriteLine("poor primer");
                }
                else
                    SeqFiles.WriteRead(filteredReads, header, read, SeqFiles.formatFNA, 120);
                readsWritten++;
            }

            Console.WriteLine("wrote " + readsWritten + " " + tag + " reads");
            return readsWritten;
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

        private static void AddMerToTable(ulong kMer, Dictionary<ulong, int> kMerTable)
        {
            if (kMerTable.ContainsKey(kMer))
                kMerTable[kMer]++;
            else
                kMerTable.Add(kMer, 1);
        }

        private static ulong Canonical(ulong kMer, int kMerSize)
        {
            ulong rcMer = kMers.ReverseComplement(kMer, kMerSize);
            if (rcMer < kMer)
                kMer = rcMer;
            return kMer;
        }

        private static ulong HashContext(string seq, int m, int contextLength, int kMerSize)
        {
            ReadOnlySpan<char> seqSpan = seq.AsSpan();
            ReadOnlySpan<char> context = seqSpan.Slice(m, contextLength-1);
            int contextHash = string.GetHashCode(context);

            ulong kMer;
            kMers.CondenseMer(seq, m + contextLength - kMerSize, kMerSize, out kMer);

            return ((ulong)contextHash << 32) | (kMer & 0xffffffff);
        }

        private static ulong HashContext(Sequence seq, int m, int contextLength, int kMerSize)
        {
            Span<char> seqSpan = seq.Bases;
            Span<char> context = seqSpan.Slice(m, contextLength-1);
            int contextHash = string.GetHashCode(context);

            ulong kMer;
            Sequence.CondenseMer(seq, m + contextLength - kMerSize, kMerSize, out kMer);

            return ((ulong)contextHash << 32) | (kMer & 0xffffffff);
        }

        private static ulong HashContextVariant(Sequence seq, int m, int contextLength, int kMerSize, char nb)
        {
            seq.Bases[seq.Length] = nb;
            Span<char> seqSpan = seq.Bases;
            Span<char> context = seqSpan.Slice(m, contextLength - 1);
            int contextHash = string.GetHashCode(context);

            ulong kMer;
            Sequence.CondenseMer(seq, m + contextLength - kMerSize, kMerSize, out kMer);

            return ((ulong)contextHash << 32) | (kMer & 0xffffffff);
        }

        private static ulong HashContextVariant(Sequence seq, int m, int contextLength, int kMerSize, ulong kMer)
        {
            Span<char> seqSpan = seq.Bases;
            Span<char> context = seqSpan.Slice(m, contextLength - 1);
            int contextHash = string.GetHashCode(context);

            return ((ulong)contextHash << 32) | (kMer & 0xffffffff);
        }

        private static ulong HashContextVariant(string seq, int m, int contextLength, int kMerSize, ulong kMer)
        {
            ReadOnlySpan<char> seqSpan = seq.AsSpan();
            ReadOnlySpan<char> context = seqSpan.Slice(m, contextLength - 1);
            int contextHash = string.GetHashCode(context);

            return ((ulong)contextHash << 32) | (kMer & 0xffffffff);
        }

        // Tries to remove unreliable (error) kMers from the kMer table. 
        // Reads are scanned, looking for sudden drops in depths and other signs of sequencing errors.
        // Suspected erroneous kMers are marked (kMersToCull) but not actually culled until all reads have
        // been scanned. The random nature of reads means that a kMer can look to be erroneous in one read but OK in another
        // read where there is more context. kMersDeemedOK is used to hold such 'redeemed' kMers, and
        // the kMersToCull and kMersDeemedOK sets are reconciled after all reads have been scanned.
        private static int DeNoiseMerTable(List<string> selectedHeaders, List<string> selectedReads, List<int> partitionStarts, List<int> partitionEnds, HashSet<int> readsStartingWithPrimer, HashSet<int> readsEndingWithPrimer, int longestReadLength, int primerLength,
                                           int kMerSize, Dictionary<ulong, int> kMerTable, int minDepth, HashSet<ulong> kMersWithZero, ReadStats[] statsForReads, StreamWriter log, ParallelOptions threads)
        {
            Console.WriteLine("denoising kMer table");

            // keep kMers to be culled until the end of denoise - avoids rare issue with propagation of cull when first faulty kMer has already been culled
            // also gives them a chance to be redeemed if they're found to be good in another read context
            Dictionary<ulong, int> kMersToCull = new Dictionary<ulong, int>(kMerTable.Count / 2);
            // kMers deemed to be OK in the context of this read. Used to redeem to-be-culled kMers that appeared to be erroneous in the context of another read
            // solves a rare problem where a kMer is culled after it's last legitimate appearance (and so can't be redeemed)
            Dictionary<ulong, int> kMersDeemedOK = new Dictionary<ulong, int>(kMerTable.Count / 2);

            Parallel.For(0, partitionStarts.Count, threads, p =>
            //for (int p = 0; p < partitionStarts.Count; p++)
            {
                List<ulong> kMerVariants = new List<ulong>(kMerSize * 4);
                List<int> kMerVariantDepths = new List<int>(kMerSize * 4);
                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int[] kMerDepths = new int[1000];

                Dictionary<ulong, int> kMersToCullLocal = new Dictionary<ulong, int>(kMerTable.Count / partitionStarts.Count);
                Dictionary<ulong, int> kMersDeemedOKLocal = new Dictionary<ulong, int>(kMerTable.Count / partitionStarts.Count);

                Sequence seq = new Sequence(longestReadLength);                // just needed for CountFollowers but declared out here to avoid allocs on every such read

                // scan through the reads again looking for possible sequencing errors
                // triggers are sudden drops down to noise level. Low depth reads are OK. 
                // it's much more complex than this... 

                for (int r = partitionStarts[p]; r < partitionEnds[p]; r++)
                {
                    string read = selectedReads[r];
                    bool readCopiedToSeq = false;

                    if (read == null)
                        continue;

                    //if (r == 5492)
                    //    Debugger.Break();

                    //if (read == "GCCTTCGCAATCGGTGTTCTTGACGATATCTAAGCATTTCACCGCTACACCGTCAATTCCGCCCACCTCAATTTTACTCAAGCTCTGCAGTATCAATGCCACCCTCCCGGTTAAGCCGGGATATTTCAGCACTGACTTATAGAGCCGCCT")
                    //    Debugger.Break();

                    // never redeem from short reads - too unreliable - must be long enough to recover from errors in first kMer
                    bool readTooShort = read.Length < shortestContextLength;

                    // get the set of kMers from this read
                    int maxDepthFound;
                    int minDepthFound;

                    int kMerCount = GetDepthsForRead(read, kMerTable, kMerSize, ref kMersFromRead, ref kMersValid, ref kMerDepths, out maxDepthFound, out minDepthFound);

                    // find the initial average depth - exclude the start (or end) of a read if we know it contains the heavily conserved forward or reverse primers
                    int initialDepthSum = 0;
                    int initialCount = 0;
                    int initialAvgDepth = minDepth;
                    int skipLength = kMerCount / 20;    // 5%
                    if (skipLength < 1)
                        skipLength = 1;
                    int startForSum = readsStartingWithPrimer.Contains(r) ? skipLength : 0;
                    int endForSum = readsEndingWithPrimer.Contains(r) ? kMerCount - skipLength : kMerCount;

                    maxDepthFound = 0;
                    for (int m = startForSum; m < endForSum; m++)
                    {
                        int kMerDepth = kMerDepths[m];
                        if (kMerDepth >= minDepth && !kMersWithZero.Contains(kMersFromRead[m]))
                        {
                            initialCount++;
                            initialDepthSum += kMerDepth;
                        }
                        if (kMerDepth > maxDepthFound)
                            maxDepthFound = kMerDepth;
                    }
                    if (initialCount > 0)
                        initialAvgDepth = initialDepthSum / initialCount;

                    // find the/a (non-noise) kMer closest in depth to the average (roughly a median) - and use log distance to help compress distances at higher depths
                    double smallestDistance = double.MaxValue;
                    int smallestDistanceDepth = 0;
                    int smallestIdx = 0;
                    int minDepthForNoise = RoundDiv(initialAvgDepth, errorRate);
                    for (int m = 0; m < kMerCount; m++)
                    {
                        int kMerDepth = kMerDepths[m];
                        double distance = Math.Log10(Math.Abs(initialAvgDepth - kMerDepth));
                        if (kMerDepth > minDepthForNoise && distance < smallestDistance)
                        {
                            smallestDistance = distance;
                            smallestDistanceDepth = kMerDepth;
                            smallestIdx = m;
                        }
                    }

                    // what variants exist for this 'median' kMer?
                    GetNextMerVariantsAndCounts(kMersFromRead[smallestIdx], kMerTable, kMerSize, kMerVariants, kMerVariantDepths);

                    // what would 1% of its depth look like?
                    int regulationMinDepth = Math.Max(RoundDiv(smallestDistanceDepth, errorRate), minDepth);
                    // and get an idea of the noise level from this deepest one
                    int lowestVariantOfDeepest = smallestDistanceDepth;
                    int highestCloseToLowest = int.MaxValue;
                    int noiseVariantOfDeepest = 0;

                    // find the lowest variant depth that looks to be noise-like
                    for (int b = 0; b < kMerVariants.Count; b++)
                    {
                        int variantDepth = kMerVariantDepths[b];
                        if (variantDepth > 0 && variantDepth < lowestVariantOfDeepest)
                        { 
                            if (CloseOrLower(variantDepth, regulationMinDepth) && !Close(variantDepth, smallestDistanceDepth))
                                noiseVariantOfDeepest = variantDepth;
                            lowestVariantOfDeepest = variantDepth;
                        }
                    }

                    // found a suitable noise-like lowest, so look for anything very close (and also noise-like)
                    if (noiseVariantOfDeepest > 0 && noiseVariantOfDeepest < smallestDistanceDepth)
                    {
                        // look for the highest depth close to this lowest one 
                        highestCloseToLowest = noiseVariantOfDeepest;
                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            int variantDepth = kMerVariantDepths[b];
                            if (variantDepth > highestCloseToLowest && variantDepth < smallestDistanceDepth && RelativelyClose(variantDepth, regulationMinDepth, smallestDistanceDepth))
                                highestCloseToLowest = variantDepth;
                        }
                    }
                    else
                    // no higher 'noise' level found - choose half the lowest variant if there was a lower variant, otherwise go with the regulation depth
                    {
                        highestCloseToLowest = lowestVariantOfDeepest == smallestDistanceDepth ? regulationMinDepth : lowestVariantOfDeepest / 2;
                    }

                    // scan along the read to find something close but greater than this minimal variant
                    for (int m = 0; m < kMerCount; m++)
                    {
                        int kMerDepth = kMerDepths[m];
                        if (kMerDepth < regulationMinDepth && kMerDepth > highestCloseToLowest && RelativelyClose(kMerDepth, highestCloseToLowest, smallestDistanceDepth))
                            highestCloseToLowest = kMerDepth;
                    }

                    // and set the min depth for the read based on this lowest variant of the deepest
                    // kMer depths must be greater than this min depth
                    int minDepthForRead = highestCloseToLowest;
                    if (highestCloseToLowest < regulationMinDepth && !(readsStartingWithPrimer.Contains(r) || readsEndingWithPrimer.Contains(r)))
                        minDepthForRead = regulationMinDepth;
                    // and remember the highest culled depth as this may become the minReadDepth after the read has been denoised
                    int maxCulledDepth = 0;
                    // and where we started culling in case we need to re-check the start of a read
                    int firstCulledIdx = -1;

                    // calculate the mean depths for the read 
                    // exclude the start (or end) of a read if we know it contains the heavily conserved forward or reverse primers
                    int totalDepthsForRead = 0;
                    int countForRead = 0;
                    double invDepthSumForRead = 0.0;
                    double sumInvFirst = 0.0;
                    int countFirst = 0;
                    double sumInvSecond = 0.0;
                    int countSecond = 0;

                    for (int m = startForSum; m < endForSum; m++)
                        if (kMersValid[m])
                        {
                            int kMerDepth = kMerDepths[m];
                            if (kMerDepth > minDepthForRead)
                            {
                                countForRead++;
                                totalDepthsForRead += kMerDepth;
                                invDepthSumForRead += 1.0f / (double)kMerDepth;

                                if (m < kMerCount / 2)
                                {
                                    countFirst++;
                                    sumInvFirst += 1.0f / (double)kMerDepth;
                                }
                                else
                                {
                                    countSecond++;
                                    sumInvSecond += 1.0f / (double)kMerDepth;
                                }
                            }
                        }

                    int avgDepthForRead = minDepthForRead;
                    int meanDepthForRead = minDepthForRead;
                    if (countForRead > 0)
                    {
                        avgDepthForRead = totalDepthsForRead / countForRead;
                        meanDepthForRead = (int)((double)countForRead / invDepthSumForRead);
                    }
                    int meanDepthFirst = 0;
                    if (countFirst > 0)
                        meanDepthFirst = (int)((double)countFirst / sumInvFirst);
                    int meanDepthSecond = 0;
                    if (countSecond > 0)
                        meanDepthSecond = (int)((double)countSecond / sumInvSecond);

                    // does the start of the read look dubious? (and treat short reads as always suspicious as they may not be long enough to have recovered from errors)
                    bool startOfReadDubious = false;
                    int deepestVariantDepthAtStart = 0;
                    if ((kMerDepths[0] < meanDepthForRead / 2) || readTooShort || kMerDepths[0] < RoundDiv(maxDepthFound, errorRate))
                    {
                        // is there a plausible alternative for the start of the read?
                        deepestVariantDepthAtStart = GetDeepestVariantDepth(VariantsWanted.allSingle, kMersFromRead[0], kMerTable, kMerSize, kMerVariants, kMerVariantDepths, minDepth);

                        // and if single-subs don't look promising, try double-subs
                        if (deepestVariantDepthAtStart < avgDepthForRead && !Close(deepestVariantDepthAtStart, avgDepthForRead))
                            deepestVariantDepthAtStart = GetDeepestVariantDepth(VariantsWanted.allDouble, kMersFromRead[0], kMerTable, kMerSize, kMerVariants, kMerVariantDepths, minDepth);

                        // declare the start dubious if it looks like we could do much better 
                        startOfReadDubious = kMerVariants.Count == 0 ||
                                             (kMerVariants.Count >= 1 && deepestVariantDepthAtStart > kMerDepths[0] && RelativelyClose(minDepthForRead, kMerDepths[0], deepestVariantDepthAtStart) && kMerDepths[0] < RoundDiv(deepestVariantDepthAtStart, errorRate));
                    }

                    // set the initial depth for the read scan
                    int previousGoodDepth = startOfReadDubious ? deepestVariantDepthAtStart : kMerDepths[0];
                    if (previousGoodDepth <= minDepthForRead)
                        previousGoodDepth = meanDepthForRead;

                    // errors taint 'k' consecutive kMers so don't redeem too close to a culled kMer
                    int lastCulledIdx = -1;
                    // was any kMer redeemed while scanning this read (logging)
                    bool kMersRedeemedInRead = false;
                    // track what's been passed as OK (and not culled) - used to determine the endpoint for cull propagation
                    int kMersOK = 0;
                    double invOKSum = 0.0;
                    bool previousCulled = startOfReadDubious;

                    // now look at every kMer to see if it looks OK or may need to be culled
                    for (int m = 0; m < kMerCount; m++)
                    {
                        bool kMerOK = false;

                        ulong kMer = kMersFromRead[m];
                        if (!kMersValid[m])
                            continue;

                        int currentDepth = kMerDepths[m];
                        ulong kMerCanonical = Canonical(kMer, kMerSize);   // just in case we need to check 

                        //if (kMerCanonical == 0xBD98E5E1AC5AB8A0)
                        //    Debugger.Break();

                        // nothing much has changed... or we've gone up rather than down - and never assume a (n,0) kMer is OK
                        if (!kMerOK && currentDepth > minDepthForRead && CloseOrHigher(currentDepth, previousGoodDepth) && !kMersWithZero.Contains(kMerCanonical))
                        {
                            kMerOK = true;
                        }

                        // check against variant kMers and decide if the kMer is dubious enough to deserve culling
                        if (!kMerOK)
                        {
                            // find the deepest alternatives in case they look better
                            if (m == 0)
                            {
                                // may have just calculated generated the alternatives for the initial kMer, so use these if we can
                                if (deepestVariantDepthAtStart == 0)
                                    GetAllMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths, 1);
                            }
                            else
                                GetNextMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths);

                            int deepestVariantDepth = 0;            // depth of deepest variant
                            ulong deepestVariant = 0;               // and which kMer it is
                            int noOfHigherVariants = 0;

                            // and find the highest and total depth of the variants
                            int totalVariantDepths = 0;
                            for (int b = 0; b < kMerVariants.Count; b++)
                            {
                                int variantDepth = kMerVariantDepths[b];
                                totalVariantDepths += variantDepth;
                                if (variantDepth > deepestVariantDepth)
                                {
                                    deepestVariantDepth = variantDepth;
                                    deepestVariant = kMerVariants[b];
                                }
                                if (variantDepth > currentDepth)
                                    noOfHigherVariants++;
                            }

                            // already know that the current kMer isn't close to the previous kMer depth (and greater than the read min depth)
                            // so cull it it if it looks like a low frequency error variant
                            bool redeemingQualities = false;
                            bool looksTooLow = currentDepth <= RoundDiv(previousGoodDepth, errorRate) ||                        // sudden drop 
                                               currentDepth <= RoundDiv(totalVariantDepths, errorRate) ||                       // better alternatives
                                               kMersWithZero.Contains(kMerCanonical) ||                                         // (n,0) counts
                                               (m == 0 && startOfReadDubious) ||                                                // dubious start of read (and we're at the start)
                                               currentDepth < minDepthForRead ||                                                // lower than the estimated min depth 
                                               currentDepth <= maxCulledDepth ||                                                // previous to-be-culled decision
                                               Crater(m, kMerDepths, kMerCount, kMerSize, minDepthForRead, previousGoodDepth);  // start of a kMerSize depth crater

                            // looks like we may want to cull this kMer so check for reasons why we might not want to...
                            if (looksTooLow)
                            {
                                // if there are alternatives, check them for followers in case the status quo looks good by this metric
                                bool statusQuoFollowers = false;
                                // only need to do this if there is an alternative (and we're not too close to the end of the read)
                                if (noOfHigherVariants > 0 && m < kMerCount - 5)
                                {
                                    // ** replace with Span<t> as part of the migration **
                                    if (!readCopiedToSeq)
                                    {
                                        seq.CopyFrom(read);
                                        readCopiedToSeq = true;
                                    }
                                    //List<int> matchesAt = new List<int>();
                                    //List<ulong> matchesFound = new List<ulong>();

                                    int bestFollowers = 0;

                                    bool currentAbandoned;
                                    int kMerFollowers = CountFollowers(seq, m, kMerCount, minDepthForRead, kMer, kMerSize, kMerTable, kMersWithZero, 0, 0, out currentAbandoned/*, matchesAt, matchesFound*/);
                                    if (kMerFollowers > bestFollowers)
                                        bestFollowers = kMerFollowers;
                                    //if (currentAbandoned)
                                    //    kMerFollowers = 0;

                                    bool deepestAbandoned;
                                    int deepestFollowers = CountFollowers(seq, m, kMerCount, minDepthForRead, deepestVariant, kMerSize, kMerTable, kMersWithZero, 0, 0, out deepestAbandoned/*, matchesAt, matchesFound*/);
                                    if (deepestFollowers > bestFollowers)
                                        bestFollowers = deepestFollowers;
                                    //if (deepestAbandoned)
                                    //    deepestFollowers = 0;

                                    //statusQuoFollowers = !currentAbandoned && kMerFollowers > deepestFollowers;
                                    int bestPossibleFollowers = kMerCount - m;
                                    statusQuoFollowers = (bestFollowers > 5) && kMerFollowers >= deepestFollowers;

                                    //DumpFollowingMatches(kMer, currentDepth, kMerSize, kMerFollowers, matchesAt, matchesFound, kMerTable);
                                    //DumpFollowingMatches(deepestVariant, deepestVariantDepth, kMerSize, kMerFollowers, matchesAt, matchesFound, kMerTable);
                                }

                                int meanOK = kMersOK > 0 ? (int)((double)kMersOK / invOKSum) : meanDepthForRead;

                                redeemingQualities = (currentDepth >= minDepthForRead && CloseOrHigher(currentDepth, previousGoodDepth)) ||
                                                     (noOfHigherVariants > 0 && VeryClose(currentDepth, deepestVariantDepth)) ||
                                                     CloserToTarget(currentDepth, minDepthForRead, previousGoodDepth) ||
                                                     (m > 0 && !previousCulled && deepestVariantDepth == currentDepth && currentDepth > 1) ||
                                                     (currentDepth >= minDepthForRead && Close(currentDepth, meanOK)) ||
                                                     statusQuoFollowers;
                            }

                            // should this kMer be culled (as it occurs in the context of this read)
                            if (looksTooLow && !redeemingQualities)
                            {
                                // zap the depth for the broken kMer (for later stats calculation for this read)
                                kMerDepths[m] = 0;
                                previousCulled = true;

                                // starting a culling region (consecutive kMers impacted by the RHS base in this kMer)
                                lastCulledIdx = m;
                                if (currentDepth > maxCulledDepth)
                                    maxCulledDepth = currentDepth;
                                if (currentDepth > minDepthForRead)
                                    minDepthForRead = currentDepth;
                                if (firstCulledIdx < 0)
                                    firstCulledIdx = m;

                                if (log != null)
                                    lock (log)
                                    {
                                        log.Write(r + ":\t" + "kdc@ " + m + "/" + kMerCount + "  " + kMers.ExpandMer(kMer, kMerSize) + " [" + kMers.ExpandMer(Canonical(kMer, kMerSize), kMerSize) + "] (" + currentDepth + ") ");
                                        for (int b = 0; b < kMerVariantDepths.Count; b++)
                                            log.Write(kMerVariantDepths[b] + " ");
                                        log.Write("pgd=" + previousGoodDepth);
                                        log.WriteLine();
                                    }

                                //if (kMerCanonical == 0x0303BE03007A375D)
                                //    Debugger.Break();

                                if (kMersToCullLocal.ContainsKey(kMerCanonical))
                                    kMersToCullLocal[kMerCanonical]++;
                                else
                                    kMersToCullLocal.Add(kMerCanonical, 1);

                                // conservatively propagate cull downstream 
                                // ----------------------------------------
                                // will continue removing kMers while they are very close in depth to the one just culled (could go to the end of the read)
                                int firstPropaguleM = -1;
                                int lastPropaguleM = -1;
                                int firstPropaguleDepth = 0;
                                int lastPropaguleDepth = 0;
                                string firstPropagule = null;
                                string lastPropagule = null;

                                int endOfCrater = 0;

                                // calculate where the crater caused by this faulty kMer ends
                                for (endOfCrater = m + 1; endOfCrater < kMerCount; endOfCrater++)
                                {
                                    int nextDepth = kMerDepths[endOfCrater];
                                    bool stillTooLow = nextDepth <= currentDepth || (RelativelyClose(nextDepth, currentDepth, avgDepthForRead) && nextDepth < deepestVariantDepth);
                                    if (!stillTooLow)
                                        break;
                                }

                                // propagate the cull until we hit the end of the crater
                                for (m++; m < endOfCrater; m++)
                                {
                                    // tile next kmer and get its depth
                                    ulong kMerPropagule = kMersFromRead[m];
                                    kMerCanonical = Canonical(kMerPropagule, kMerSize);
                                    currentDepth = kMerDepths[m];

                                    kMerDepths[m] = 0;

                                    if (currentDepth > maxCulledDepth)
                                        maxCulledDepth = currentDepth;
                                    if (currentDepth > minDepthForRead)
                                        minDepthForRead = currentDepth;

                                    //if (kMerCanonical == 0xCAC2BFF66F9360F0)
                                    //    Debugger.Break();

                                    if (kMersToCullLocal.ContainsKey(kMerCanonical))
                                        kMersToCullLocal[kMerCanonical]++;
                                    else
                                        kMersToCullLocal.Add(kMerCanonical, 1);

                                    if (log != null)
                                    {
                                        if (firstPropaguleM < 0)
                                        {
                                            firstPropaguleM = m;
                                            firstPropagule = kMers.ExpandMer(kMerPropagule, kMerSize);
                                            firstPropaguleDepth = currentDepth;
                                        }
                                        lastPropaguleM = m;
                                        lastPropagule = kMers.ExpandMer(kMerPropagule, kMerSize);
                                        lastPropaguleDepth = currentDepth;
                                    }

                                } // propagation loop

                                // 'm' finishes at first base after the end of the crater and the main loop increments 'm' immediately, so decrement 'm' now
                                // so that the main loop looks next at the next unculled kMer
                                m--;

                                if (maxCulledDepth > minDepthForRead)
                                    minDepthForRead = maxCulledDepth;

                                if (log != null && firstPropaguleM >= 0)
                                    lock (log)
                                    {
                                        log.WriteLine(r + ":\t" + "kdp: " + firstPropagule + "(" + firstPropaguleDepth + ") [@" + firstPropaguleM + "] - " + lastPropagule + "(" + lastPropaguleDepth + ") [@" + lastPropaguleM + "] " + "pgd=" + previousGoodDepth);
                                    }
                            } // kMer looked dubious and we decided to cull it (and propagate the cull)
                            else
                            {
                                // decided the kMer was actually OK after checking, so move on
                                kMerOK = true;
                            }
                        } // kMer needed to be examined a bit more closely

                        if (kMerOK)
                        {
                            //if (kMerCanonical == 0x0303BE03007A375D)
                            //    Debugger.Break();

                            //if (r == 4211 || r == 10367)
                            //    Debugger.Break();

                            // decided the kMer is OK in the context of this read - could have been seen as poor in another context, so possibly redeem it in that case.
                            // never redeem from short reads as it's too easy for an error to make the whole read look like a low abundance strain/variant
                            // and don't redeem if we're too close to the last culled kMer
                            int distanceFromLastCull = lastCulledIdx >= 0 ? m - lastCulledIdx : kMerSize;
                            if (!readTooShort && distanceFromLastCull >= kMerSize / 4)
                            {
                                if (kMersDeemedOKLocal.ContainsKey(kMerCanonical))
                                    kMersDeemedOKLocal[kMerCanonical]++;
                                else
                                    kMersDeemedOKLocal.Add(kMerCanonical, 1);

                                //if (kMerCanonical == 0x874A960957F46F61)
                                //    Debugger.Break();

                                if (log != null && kMersToCullLocal.ContainsKey(kMerCanonical))
                                {
                                    lock (log)
                                    {
                                        log.WriteLine(r + ":\t" + "kdm: " + kMers.ExpandMer(kMerCanonical, kMerSize) + " [" + kMers.ExpandMer(kMers.ReverseComplement(kMerCanonical, kMerSize), kMerSize) + "] @" + m + " d=" + currentDepth);
                                    }
                                    kMersRedeemedInRead = true;
                                }
                            }

                            if (startOfReadDubious && m == 0)
                                startOfReadDubious = false;

                            if (m > startForSum)
                            {
                                kMersOK++;
                                invOKSum += 1.0f / (double)currentDepth;
                            }
                            previousGoodDepth = currentDepth;
                            previousCulled = false;

                            if (currentDepth < minDepthForRead && currentDepth > maxCulledDepth)
                                minDepthForRead = currentDepth;
                        }
                    } // for every kMer in the read

                    if (kMersRedeemedInRead && log != null)
                        lock (log)
                        {
                            log.WriteLine(r + ":\t" + "rdm: " + read);
                        }

                    // it's possible we let some kMers pass through at the start of the read that would have been culled if they'd been encountered later
                    if (startOfReadDubious)
                        for (int m = 0; m < firstCulledIdx; m++)
                        {
                            int kMerDepth = kMerDepths[m];
                            if (kMerDepth <= maxCulledDepth && kMerDepth > 0)
                            {
                                kMerDepths[m] = 0;
                                ulong kMerCanonical = kMers.Canonical(kMersFromRead[m], kMerSize);
                                if (kMersDeemedOKLocal.ContainsKey(kMerCanonical))
                                {
                                    if (kMersDeemedOKLocal[kMerCanonical] > 0)
                                        kMersDeemedOKLocal[kMerCanonical]--;
                                    if (log != null)
                                        lock (log)
                                        {
                                            log.WriteLine(r + ":\t" + "kmx: " + kMers.ExpandMer(kMerCanonical, kMerSize) + " [" + kMers.ExpandMer(kMers.ReverseComplement(kMerCanonical, kMerSize), kMerSize) + "] @" + m + " d=" + kMerDepth);
                                        }
                                }
                                if (kMersToCullLocal.ContainsKey(kMerCanonical))
                                    kMersToCullLocal[kMerCanonical]++;
                                else
                                    kMersToCullLocal.Add(kMerCanonical, 1);
                            }
                        }

                    // recalculate the stats for this read and save them for later steps
                    totalDepthsForRead = 0;
                    countForRead = 0;
                    invDepthSumForRead = 0.0;
                    minDepthForRead = Math.Max(maxCulledDepth, minDepthForRead);
                    for (int m = 0; m < kMerCount; m++)
                    {
                        int kMerDepth = kMerDepths[m];
                        if (kMerDepth > minDepthForRead)
                        {
                            countForRead++;
                            totalDepthsForRead += kMerDepth;
                            invDepthSumForRead += 1.0f / (double)kMerDepth;
                        }
                    }

                    if (countForRead > 0)
                    {
                        avgDepthForRead = totalDepthsForRead / countForRead;
                        meanDepthForRead = (int)((double)countForRead / invDepthSumForRead);
                    }
                    int initialGoodDepth = Math.Max(meanDepthForRead, startOfReadDubious ? deepestVariantDepthAtStart : kMerDepths[0]);

                    // save the per-read stats for now. These will be recalculated after the cull/redeem step has driven some kMer depths to zero
                    ReadStats statsForRead = new ReadStats();
                    statsForRead.avgDepth = avgDepthForRead;
                    statsForRead.meanDepth = meanDepthForRead;
                    statsForRead.minDepthAllowed = minDepthForRead;
                    statsForRead.initialGoodDepth = initialGoodDepth;
                    statsForReads[r] = statsForRead;

                } // for each read in the partition

                // merge per-thread kMers to cull or keep into the global sets
                lock (kMersToCull)
                {
                    foreach (KeyValuePair<ulong, int> kvp in kMersToCullLocal)
                    {
                        if (kMersToCull.ContainsKey(kvp.Key))
                            kMersToCull[kvp.Key] += kvp.Value;
                        else
                            kMersToCull.Add(kvp.Key, kvp.Value);
                    }
                    foreach (KeyValuePair<ulong, int> kvp in kMersDeemedOKLocal)
                    {
                        if (kMersDeemedOK.ContainsKey(kvp.Key))
                            kMersDeemedOK[kvp.Key] += kvp.Value;
                        else
                            kMersDeemedOK.Add(kvp.Key, kvp.Value);
                    }
                }
            });
            //}
            //for each partition

            // remove the detected unsound kMers from the table (unless they're redeemed)
            int kMersCulled = 0;
            foreach (KeyValuePair<ulong, int> kvp in kMersToCull)
            {
                ulong kMerToCull = kvp.Key;
                int cullCount = kvp.Value;
                int redeemedCount;
                kMersDeemedOK.TryGetValue(kMerToCull, out redeemedCount);

                //if (kMerToCull == 0x046B521D71A8A492)
                //    Debugger.Break();

                if (cullCount > redeemedCount * 5 && CloseToNoise(redeemedCount, Math.Max(minDepth, RoundDiv(redeemedCount, errorRate))))
                {
                    kMerTable[kMerToCull] = 0;
                    kMersCulled++;
                }
            }

            return kMersCulled;
        }

        // finds the (tagged) primer reads from selectedReads - called at the start of the extending phase to get the sets of starting reads 
        private static void ExtractPrimerReads(List<string> selectedHeaders, List<string> selectedReads, List<int> partitionStarts, List<int> partitionEnds,
                                               int startingReadTrim, int forwardPrimerLength, List<int> readsWithStartingPrimers, HashSet<int> startingReadRCed, Dictionary<string, int> startsOfFPReads, Dictionary<int, string> FPReadPrimers,
                                               List<int> nonStartingReads, HashSet<int> readsStartingWithPrimer, HashSet<int> readsEndingWithPrimer, out int longestSelectedRead, ParallelOptions threads)
        {
            Console.WriteLine("find reads with primers");

            int longestSelectedReadFound = 0;

            Parallel.For(0, partitionStarts.Count, threads, p =>
            //for (int p = 0; p < partitionStarts.Count; p++)
            {
                int longestSelectedReadLocal = 0;
                List<int> readsWithStartingPrimersLocal = new List<int>();
                HashSet<int> readsStartingWithPrimerLocal = new HashSet<int>();
                HashSet<int> readsEndingWithPrimerLocal = new HashSet<int>();
                HashSet<int> startingReadRCedLocal = new HashSet<int>();
                List<int> nonStartingReadsLocal = new List<int>();

                for (int r = partitionStarts[p]; r < partitionEnds[p]; r++)
                {
                    string header = selectedHeaders[r];
                    string read = selectedReads[r];

                    if (read == null)
                        continue;

                    if (read.Length > longestSelectedReadLocal)
                        longestSelectedReadLocal = read.Length;

                    //if (read == "GGGGCAAGCGTTATCCGGAATTACTGGGCGTAAAGGGTCCGTAGGCGGCTAGTTAAGTCGAGGTTAAAAGGCAGTAGCTCAACTACTGTTGGGCCTTGAAACTAATTAGCTTGAGTATAGGAGAGGAAAGTGGAATTCCCGGAGTAGCGG")
                    //    Debugger.Break();

                    if (header.EndsWith(";FP"))
                    {
                        readsWithStartingPrimersLocal.Add(r);
                        readsStartingWithPrimerLocal.Add(r);
                    }
                    else if (header.EndsWith(";FP'"))
                    {
                        selectedReads[r] = kMers.ReverseComplement(read);
                        readsWithStartingPrimersLocal.Add(r);
                        readsStartingWithPrimerLocal.Add(r);
                        startingReadRCedLocal.Add(r);
                    }
                    else
                    {
                        if (header.EndsWith(";RP"))
                            readsStartingWithPrimerLocal.Add(r);
                        if (header.EndsWith(";RP'"))
                            readsEndingWithPrimerLocal.Add(r);
                        nonStartingReadsLocal.Add(r);
                    }
                } // for each read in the partition

                // merge the per-partition tables/lists into the global ones
                lock (readsWithStartingPrimers)
                {
                    foreach (int r in readsWithStartingPrimersLocal)
                        readsWithStartingPrimers.Add(r);
                    foreach (int r in readsStartingWithPrimerLocal)
                        readsStartingWithPrimer.Add(r);
                    foreach (int r in readsEndingWithPrimerLocal)
                        readsEndingWithPrimer.Add(r);
                    foreach (int r in nonStartingReadsLocal)
                        nonStartingReads.Add(r);
                    foreach (int r in startingReadRCedLocal)
                        startingReadRCed.Add(r);
                    if (longestSelectedReadLocal > longestSelectedReadFound)
                        longestSelectedReadFound = longestSelectedReadLocal;
                }
            });
            //} // for each partition

            longestSelectedRead = longestSelectedReadFound;

            // remember what's at the start of the starting reads so we can look for other reads with only partial forward primers
            foreach (int r in readsWithStartingPrimers)
            {
                string startingRead = selectedReads[r];

                //if (r == 12813)
                //    Debugger.Break();

                // tile for all shortest-context fragments, starting from the second base of the forward primer and going until half-way through the primer
                // remember the fragment and where in the primer it started (these will always contain at least the second half of the primer)
                for (int i = 1; i < startingReadTrim; i++)
                {
                    if (i + shortestContextSize > startingRead.Length)
                        break;

                    string startOfRegion = startingRead.Substring(i, shortestContextSize);
                    if (!startsOfFPReads.ContainsKey(startOfRegion))
                        startsOfFPReads.Add(startOfRegion, i);
                    //if (startOfRegion == "CAGGGTATCTAATCTTGTTCGCTCCCCATGCTTTCGCTCC")
                    //    Debugger.Break();
                }
                string startingReadPrimer = startingRead.Substring(0, forwardPrimerLength);
                FPReadPrimers.Add(r, startingReadPrimer);

                // partially trim the forward primer and save the rest of the read (at least short kMer context in length)
                startingRead = startingRead.Substring(startingReadTrim);
                selectedReads[r] = startingRead;
            }
        }

        // build initial kMer table, containing all kMers tiled from all selected reads. This initial table is then denoised.
        private static void GenerateInitialkMerTable(List<string> selectedReads, HashSet<int> startingReadRCed, List<int> partitionStarts, List<int> partitionEnds, Dictionary<ulong, int> kMerTable, int kMerSize, HashSet<int> nonACGTReads, HashSet<ulong> kMersWithZero, ParallelOptions threads)
        {
            Console.WriteLine("building initial kMer table for extension phase");

            Dictionary<ulong, int> kMerTableAsRead = new Dictionary<ulong, int>(1000000);

            Parallel.For(0, partitionStarts.Count, threads, p =>
            //for (int p = 0; p < partitionStarts.Count; p++)
            {
                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int[] kMerDepths = new int[1000];

                // non-canonical kMers+counts
                Dictionary<ulong, int> kMerTableLocal = new Dictionary<ulong, int>(100000);
                HashSet<int> nonACGTReadsLocal = new HashSet<int>();

                for (int r = partitionStarts[p]; r < partitionEnds[p]; r++)
                {
                    string read = selectedReads[r];

                    if (read == null)
                        continue;

                    // some starting reads have been reversed already - so reverse them again to get accurate per-form counts
                    if (startingReadRCed.Contains(r))
                        read = kMers.ReverseComplement(read);

                    //if (read == "GGGGCAAGCGTTATCCGGAATTACTGGGCGTAAAGGGTCCGTAGGCGGCTAGTTAAGTCGAGGTTAAAAGGCAGTAGCTCAACTACTGTTGGGCCTTGAAACTAATTAGCTTGAGTATAGGAGAGGAAAGTGGAATTCCCGGAGTAGCGG")
                    //    Debugger.Break();

                    // initial set of kMers for this partition
                    int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);
                    bool allACGT = true;

                    for (int k = 0; k < kMerCount; k++)
                    {
                        if (kMersValid[k])
                        {
                            //if (kMersFromRead[k] == 0x667F1952C3D60419 || kMersFromRead[k] == 0x9BEF683C7A9B0266 || kMersFromRead[k] == 0x5C6667F1952C3D60 || kMersFromRead[k] == 0xF683C7A9B02666CA)
                            //{
                            //    Console.WriteLine(kMersFromRead[k].ToString("X16") + " " + kMers.ExpandMer(kMersFromRead[k], 32) + " from " + header);
                            //    Debugger.Break();
                            //}
                            AddMerToTable(kMersFromRead[k], kMerTableLocal);
                            //if (kMersFromRead[k] == 0xEC079FFB1A0600A7)
                            //    Debugger.Break();
                        }
                        else
                            allACGT = false;
                    }

                    if (!allACGT)
                    {
                        nonACGTReadsLocal.Add(r);
                        continue;
                    }

                } // for each read in the partition

                // merge the per-partition tables/lists into the global ones
                lock (kMerTableAsRead)
                {
                    foreach (KeyValuePair<ulong, int> kvp in kMerTableLocal)
                    {
                        if (kMerTableAsRead.ContainsKey(kvp.Key))
                            kMerTableAsRead[kvp.Key] += kvp.Value;
                        else
                            kMerTableAsRead.Add(kvp.Key, kvp.Value);
                    }

                    foreach (int r in nonACGTReadsLocal)
                        nonACGTReads.Add(r);
                }
            });
            //} // for each partition

            // build the canonical table from the as-read table, and remember any kMers only found in one form (N,0 counts)
            foreach (KeyValuePair<ulong, int> kvp in kMerTableAsRead)
            {
                //if (kvp.Key == 0x25FF6F5B10092FC4)
                //    Debugger.Break();
                ulong canonicalMer = Canonical(kvp.Key, kMerSize);
                if (!kMerTable.ContainsKey(canonicalMer))
                    kMerTable.Add(canonicalMer, 0);
                kMerTable[canonicalMer] += kvp.Value;
                ulong rcMer = kMers.ReverseComplement(kvp.Key, kMerSize);
                if (!kMerTableAsRead.ContainsKey(rcMer))
                    kMersWithZero.Add(canonicalMer);
            }

            // only exclude (<=2,0) kMers
            List<ulong> keepers = new List<ulong>();
            foreach (ulong kMer in kMersWithZero)
                if (kMerTable[kMer] > 2)
                    keepers.Add(kMer);
            foreach (ulong kMer in keepers)
                kMersWithZero.Remove(kMer);
        }

        // Generate the set of kMer context from the selected reads - for a given context-length - returns the average coverage depth
        private static int GenerateContextTable(int contextLength, List<string> selectedHeaders, List<string> selectedReads, List<int> selectedPartitionStarts, List<int> selectedPartitionEnds,
                                                 HashSet<int> nonACGTReads, Dictionary<ulong, int> kMerTable, int kMerSize, ReadStats[] statsForReads, Dictionary<int, Dictionary<ulong, int>> kMerContexts, ParallelOptions threads)
        {
            Dictionary<ulong, int> contextTable = kMerContexts[contextLength];

            //for (int srp = 0; srp < selectedPartitionStarts.Count; srp++)
            Parallel.For(0, selectedPartitionStarts.Count, threads, srp =>
            {
                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int[] kMerDepths = new int[1000];

                Dictionary<ulong, int> contextTableLocal = new Dictionary<ulong, int>();

                for (int r = selectedPartitionStarts[srp]; r < selectedPartitionEnds[srp]; r++)
                {
                    string header = selectedHeaders[r];
                    string read = selectedReads[r];

                    if (read == null)
                        continue;

                    if (!nonACGTReads.Contains(r))
                    {
                        //if (r == 13630)
                        //    Debugger.Break();

                        int minDepthFound;
                        int maxDepthFound;
                        int kMerCount = GetDepthsForRead(read, kMerTable, kMerSize, ref kMersFromRead, ref kMersValid, ref kMerDepths, out maxDepthFound, out minDepthFound);
                        // zap any kMer depths lower than the (revised) min for the read
                        for (int k = 0; k < kMerCount; k++)
                            if (kMerDepths[k] < statsForReads[r].minDepthAllowed)
                                kMerDepths[k] = 0;

                        int contextCount = read.Length - contextLength + 1;

                        // forward read contexts
                        for (int p = 0; p < contextCount; p++)
                        {
                            // don't generate contexts that appear to have dubious kMers
                            for (int m = p; m < contextLength - kMerSize + 1; m++)
                                if (kMerDepths[m] == 0)
                                    break;

                            ulong hashedContext = HashContext(read, p, contextLength, kMerSize);
                            if (contextTableLocal.ContainsKey(hashedContext))
                                contextTableLocal[hashedContext]++;
                            else
                                contextTableLocal.Add(hashedContext, 1);
                        }

                        // and for the RCed read
                        string rcRead = kMers.ReverseComplement(read);
                        for (int p = 0; p < contextCount; p++)
                        {
                            // don't generate contexts that appear to have dubious kMers
                            for (int m = p; m < contextLength - kMerSize + 1; m++)
                                if (kMerDepths[m] == 0)
                                    break;

                            ulong hashedContext = HashContext(rcRead, p, contextLength, kMerSize);
                            if (contextTableLocal.ContainsKey(hashedContext))
                                contextTableLocal[hashedContext]++;
                            else
                                contextTableLocal.Add(hashedContext, 1);
                        }
                    }
                }

                lock (contextTable)
                {
                    foreach (KeyValuePair<ulong, int> kvp in contextTableLocal)
                    {
                        if (contextTable.ContainsKey(kvp.Key))
                            contextTable[kvp.Key] += kvp.Value;
                        else
                            contextTable.Add(kvp.Key, kvp.Value);
                    }
                }
            });

            long totalCoverage = 0;
            long totalCounted = 0;
            int avgCoverage = 0;

            foreach (int count in contextTable.Values)
            {
                if (count > 1)
                {
                    totalCoverage += count;
                    totalCounted++;
                }
            }

            if (totalCounted > 0)
                avgCoverage = (int)(totalCoverage / totalCounted);

            return avgCoverage;
        }

        // core Refresh Stats method. Called only from wrapper methods, never directly
        private static int RefreshStatsForRead(string read, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth,
                                                ReadStats readStats, ulong[] kMersFromRead, bool[] kMersValid, int[] kMerDepths)
        {
            int deepestDepth;
            int lowestAllowedDepth;
            int lowestDepthFound = 0;
            int initialMinDepth = minDepth;
            int kMerCount = GetDepthsForRead(read, kMerTable, kMerSize, ref kMersFromRead, ref kMersValid, ref kMerDepths, out deepestDepth, out lowestDepthFound);
            readStats.minDepthFound = lowestDepthFound;
            lowestAllowedDepth = Math.Max(lowestDepthFound/2, readStats.avgDepth / 100);
            if (lowestAllowedDepth < minDepth)
                lowestAllowedDepth = minDepth;

            // check for 0-depth kMers - indicators of an error crater - and we should be less forgiving about including low-depth kMers in the stats
            // and calculate an average that doesn't include these really low depths as well
            int notLowestSum = 0;
            int notLowestCount = 0;
            double invDepthSumForRead = 0.0;

            for (int m = 0; m < kMerCount; m++)
            {
                int kMerDepth = kMerDepths[m];

                if (!CloseOrLower(kMerDepth, lowestAllowedDepth))
                {
                    notLowestSum += kMerDepth;
                    notLowestCount++;
                    invDepthSumForRead += 1.0f / (double)kMerDepth;
                }
            }
            int avgDepth = notLowestCount > 0 ? notLowestSum / notLowestCount : Math.Max(readStats.avgDepth, minDepth);
            int meanNonZero = (int)((double)(notLowestCount) / invDepthSumForRead);
            initialMinDepth = Math.Max(RoundDiv(avgDepth, 20), minDepth);

            // does the start of the read look dubious?
            int deepestVariantDepthAtStart = kMerDepths[0];
            bool startDubious = false;
            bool endDubious = false;
            List<ulong> kMerVariants = new List<ulong>(kMerSize * 4);
            List<int> kMerVariantDepths = new List<int>(kMerSize * 4);

            if (kMerDepths[0] <= RoundDiv(deepestDepth, errorRate))
            {
                // is there a plausible alternative for the start of the read?
                GetAllMerVariantsAndCounts(kMersFromRead[0], kMerTable, kMerSize, kMerVariants, kMerVariantDepths, initialMinDepth);
                for (int b = 0; b < kMerVariants.Count; b++)
                {
                    int variantDepth = kMerVariantDepths[b];
                    if (variantDepth > deepestVariantDepthAtStart)
                        deepestVariantDepthAtStart = variantDepth;
                }
            }

            if (deepestVariantDepthAtStart > kMerDepths[0] && !Close(deepestVariantDepthAtStart, kMerDepths[0]))
            {
                startDubious = true;
            }

            // calculate mean and average for the first half of the read
            int totalDepthsForRead = 0;
            int countForRead = 0;
            invDepthSumForRead = 0.0;

            for (int m = 0; m < kMerCount / 2; m++)
            {
                int kMerDepth = kMerDepths[m];
                if (kMerDepth >= initialMinDepth)
                {
                    countForRead++;
                    totalDepthsForRead += kMerDepth;
                    invDepthSumForRead += 1.0f / (double)kMerDepth;
                }
            }
            int avgDepthFirst = 0;
            int meanDepthFirst = 0;
            if (countForRead > 0)
            {
                avgDepthFirst = totalDepthsForRead / countForRead;
                meanDepthFirst = (int)((double)countForRead / invDepthSumForRead);
            }
            else
                startDubious = true;

            // and the same calcs for the second half
            totalDepthsForRead = 0;
            countForRead = 0;
            invDepthSumForRead = 0.0;
            for (int m = kMerCount / 2; m < kMerCount; m++)
            {
                int kMerDepth = kMerDepths[m];
                if (kMerDepth >= initialMinDepth)
                {
                    countForRead++;
                    totalDepthsForRead += kMerDepth;
                    invDepthSumForRead += 1.0f / (double)kMerDepth;
                }
            }
            int avgDepthSecond = 0;
            int meanDepthSecond = 0;
            if (countForRead > 0)
            {
                avgDepthSecond = totalDepthsForRead / countForRead;
                meanDepthSecond = (int)((double)countForRead / invDepthSumForRead);
            }
            else
                endDubious = true;

            // choose an appropriate mean and average
            if ((startDubious && !Close(meanDepthFirst, avgDepthFirst)) || (endDubious && !Close(meanDepthSecond, avgDepthSecond)))
            {
                readStats.avgDepth = Math.Max(avgDepthFirst, avgDepthSecond);
                readStats.meanDepth = Math.Max(meanDepthFirst, meanDepthSecond);
            }
            else
            {
                readStats.avgDepth = Math.Min(avgDepthFirst, avgDepthSecond);
                readStats.meanDepth = Math.Min(meanDepthFirst, meanDepthSecond);
            }

            readStats.initialGoodDepth = Math.Max(meanDepthFirst, meanDepthSecond); 

            // initial guess at a min depth
            readStats.minDepthAllowed = RoundDiv(readStats.meanDepth,20);       // 5% of harmonic mean
            // reduce to be clear of the lowest (non-zero) depth found in the read 
            // (only if we don't think there's a bad patch and the lowest is dubious)
            int lowCount = kMerCount - notLowestCount;
            if (lowCount < kMerCount / 3 && readStats.minDepthAllowed > initialMinDepth && CloseOrHigher(readStats.minDepthAllowed, lowestAllowedDepth))
                readStats.minDepthAllowed = Math.Min(readStats.minDepthAllowed, lowestAllowedDepth - 1);

            if (readStats.minDepthAllowed < minDepth)
                readStats.minDepthAllowed = minDepth;

            return kMerCount;
        }

        // called on a read-by-read basis after a read has been cleaned or trimmed
        private static int RefreshStatsForRead(Sequence readSeq, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth, ReadStats readStats)
        {
            string read = readSeq.ToString();
            int kMersInRead = read.Length - kMerSize + 1;
            ulong[] kMersFromRead = new ulong[kMersInRead];
            bool[] kMersValid = new bool[kMersInRead];
            int[] kMerDepths = new int[kMersInRead];

            return RefreshStatsForRead(read, kMerTable, kMerSize, minDepth, readStats, kMersFromRead, kMersValid, kMerDepths);
        }

        // called after read extension  
        private static int RefreshStatsForRead(string read, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth, ReadStats readStats)
        {
            int kMersInRead = read.Length - kMerSize + 1;
            ulong[] kMersFromRead = new ulong[kMersInRead];
            bool[] kMersValid = new bool[kMersInRead];
            int[] kMerDepths = new int[kMersInRead];

            return RefreshStatsForRead(read, kMerTable, kMerSize, minDepth, readStats, kMersFromRead, kMersValid, kMerDepths);
        }

        // called post-denoise to better reflect counts change there
        private static void RefreshStatsForSelectedReads(List<string> selectedReads, List<int> selectedPartitionStarts, List<int> selectedPartitionEnds, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth, Dictionary<int, int> readPairs, ReadStats[] statsForReads, ParallelOptions threads)
        {
            // overlapping (trivial) read pairs will have non-zero coverage of 2 (1,1) for single errors. Min depth allowed for such reads has to be > 2.
            // Trivial pairs won't appear in readPairs - but this table is canonical so needs to be expanded
            HashSet<int> readsWithPairs = new HashSet<int>(readPairs.Count * 2);
            foreach(KeyValuePair<int, int> kvp in readPairs)
            {
                readsWithPairs.Add(kvp.Key);
                readsWithPairs.Add(kvp.Value);
            }

            //ulong[] kMersFromRead = new ulong[1000];
            //bool[] kMersValid = new bool[1000];
            //int[] kMerDepths = new int[1000];

            //for (int srp = 0; srp < selectedPartitionStarts.Count; srp++)
            Parallel.For(0, selectedPartitionStarts.Count, threads, srp =>
            {
                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int[] kMerDepths = new int[1000];

                for (int r = selectedPartitionStarts[srp]; r < selectedPartitionEnds[srp]; r++)
                {
                    //if (r == 57203 || r == 88495)
                    //    Debugger.Break();

                    // refresh the per-read stats now that the denoising has been completed
                    string read = selectedReads[r];
                    bool trivialPair = false;
                    if (readsWithPairs.Count > 0)
                        trivialPair = !readsWithPairs.Contains(r);
                    int minDepthForRead = minDepth;
                    if (trivialPair)
                        minDepthForRead++;
                    if (read != null)
                        RefreshStatsForRead(read, kMerTable, kMerSize, minDepthForRead, statsForReads[r], kMersFromRead, kMersValid, kMerDepths);
                }
            });
            //}
        }

        // called to generate final stats for trimmed starting reads
        private static void RefreshStatsForStartingReads(List<int> trimmedReads, List<string> selectedReads, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth, ReadStats[] statsForReads)
        {
            ulong[] kMersFromRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int[] kMerDepths = new int[1000];

            for (int r = 0; r < trimmedReads.Count; r++)
            {
                //if (r == 25592)
                //    Debugger.Break();

                // refresh the per-read stats now that the denoising has been completed
                string read = selectedReads[r];
                if (read != null)
                    RefreshStatsForRead(read, kMerTable, kMerSize, minDepth, statsForReads[r], kMersFromRead, kMersValid, kMerDepths);
            }
        }

        // look up kMer depths for a read in CANONICAL kMerTable
        private static int GetDepthsForRead(string read, Dictionary<ulong, int> kMerTable, int kMerSize, ref ulong[] kMersFromRead, ref bool[] kMersValid, ref int[] kMerDepths, out int deepestDepth, out int lowestDepth)
        {
            deepestDepth = 0;
            lowestDepth = int.MaxValue;
            int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);
            if (kMerCount > kMerDepths.Length)
                Array.Resize<int>(ref kMerDepths, kMerCount);

            // find all the depths
            for (int m = 0; m < kMerCount; m++)
                if (kMersValid[m])
                {
                    int kMerDepth = kMerTableLookup(kMersFromRead[m], kMerTable, kMerSize);
                    kMerDepths[m] = kMerDepth;
                    if (kMerDepth > deepestDepth)
                        deepestDepth = kMerDepth;
                    if (kMerDepth < lowestDepth)
                        lowestDepth = kMerDepth;
                }
                else
                    kMerDepths[m] = 0;

            return kMerCount;
        }

        // lookup depths for NON-canonical kMerTable (and trim kMers based on depth) 
        private static int GetDepthsForFilterRead(string read, Dictionary<ulong, int> kMerCounts, Dictionary<ulong, int> kMerCountsOtherDirection, int kMerSize, bool initialRead, ref ulong[] kMersFromRead, ref bool[] kMersValid, ref int[] kMerDepths)
        {
            int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);
            if (kMerCount > kMerDepths.Length)
                Array.Resize<int>(ref kMerDepths, kMerCount);

            // find all the depths 
            double invDepthSum = 0.0;
            int maxDepth = 0;
            for (int m = 0; m < kMerCount; m++)
                if (kMersValid[m])
                {
                    int kMerDepth = 0;
                    kMerCounts.TryGetValue(kMersFromRead[m], out kMerDepth);
                    kMerDepths[m] = kMerDepth;
                    invDepthSum += 1.0 / (double)kMerDepth;
                    if (kMerDepth > maxDepth)
                        maxDepth = kMerDepth;
                    //if (kMersFromRead[m] == 0xE83E74E4A27DEF8D || kMersFromRead[m] == 0x8D048275E4E243D4)
                    //    Debugger.Break();
                }
                else
                    kMerDepths[m] = 0;

            // scan backwards until we hit a depth > min
            int hmDepth = (int)((double)kMerCount / invDepthSum);
            int minDepth = Math.Max(hmDepth / 4, 1);
            int noiseDepth = hmDepth / 100;
            int trimmedBases = 0;

            // initial reads - treat them less harshly as they're shorter (after trimming) - and so more likely to fade away towards the end
            if (initialRead)
            {
                for (int i = kMerCount - 1; i >= 0; i--)
                {
                    if (kMerDepths[i] <= noiseDepth)
                        trimmedBases++;
                    else
                        break;
                }
            }
            else
            // all other reads
            {
                for (int i = kMerCount - 1; i >= 0; i--)
                {
                    // really looks like a sequencing or PCR error 
                    if (kMerDepths[i] <= noiseDepth)
                    {
                        trimmedBases++;
                        continue;
                    }

                    // just suspiciously low - so trim it if we haven't also found it in the collection being built in the opposite direction
                    if (kMerDepths[i] < minDepth)
                    {
                        int otherDepth;
                        bool foundOther = kMerCountsOtherDirection.TryGetValue(kMers.ReverseComplement(kMersFromRead[i], kMerSize), out otherDepth);
                        if (foundOther && otherDepth >= minDepth)
                        {
                            //Debugger.Break();
                            break;
                        }
                        else
                            trimmedBases++;
                    }
                    else
                        break;
                }
            }

            //if (maxDepth > 1000 && trimmedBases > 0 && !initialRead)
            //    Debugger.Break();

            if (trimmedBases == kMerCount)
                trimmedBases = 0;
            if (trimmedBases > kMerSize)
                trimmedBases -= kMerSize / 2;
            kMerCount = kMerCount - trimmedBases;

            return kMerCount;
        }

        // generate all 4 possible 'next' kMers from a starting kMer and retrieve their depths
        private static int GetNextMerVariantsAndCounts(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize, List<ulong> kMerVariants, List<int> kMerCounts)
        {
            ulong lbMask = 0xffffffffffffffff << (32 - kMerSize + 1) * 2;
            ulong kMerPrefix = kMer & lbMask;
            kMerVariants.Clear();
            kMerCounts.Clear();
            int deepestDepth = 0;

            // vary last base in kMer
            for (int b = 0; b < 4; b++)
            {
                ulong kMerVariant = kMerPrefix | ((ulong)b << (32 - kMerSize) * 2);
                int depth = kMerTableLookup(kMerVariant, kMerTable, kMerSize);
                kMerCounts.Add(depth);
                // always add variant to list - even if doesn't exist in the kMer table
                // the extension code uses position in the list to determine the base involved
                kMerVariants.Add(kMerVariant);

                if (depth > deepestDepth)
                    deepestDepth = depth;
            }

            return deepestDepth;
        }

        // generate all single-base-substitution variants of a kMer (and their counts). Variants with depth less than minDepth are dropped
        private static void GetAllMerVariantsAndCounts(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize, List<ulong> kMerVariants, List<int> kMerCounts, int minDepth)
        {
            ulong baseMask = 0xc000000000000000;
            kMerVariants.Clear();
            kMerCounts.Clear();
            HashSet<ulong> distinctVariants = new HashSet<ulong>();

            for (int m = 0; m < kMerSize; m++)
            {
                ulong kMerWithHole = kMer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong kMerVariant = kMerWithHole | newBase;
                    int depth = kMerTableLookup(kMerVariant, kMerTable, kMerSize);
                    if (depth >= minDepth)
                    {
                        if (distinctVariants.Add(kMerVariant))
                        {
                            kMerCounts.Add(depth);
                            kMerVariants.Add(kMerVariant);
                        }
                    }
                }
            }
        }

        // generate all double-substitution variants of a kMer (and their counts). Variants with depth less than minDepth are dropped

        private static void GetAllDoubleMerVariantsAndCounts(List<ulong> kMerVariants, List<int> kMerCounts, Dictionary<ulong, int> kMerTable, int kMerSize, int minDepth)
        {
            ulong baseMask = 0xc000000000000000;

            List<ulong> startingVariants = new List<ulong>(kMerVariants);
            HashSet<ulong> distinctVariants = new HashSet<ulong>();

            kMerVariants.Clear();
            kMerCounts.Clear();

            foreach (ulong kMer in startingVariants)
                for (int m = 0; m < kMerSize; m++)
                {
                    ulong kMerWithHole = kMer & ~(baseMask >> (m * 2));
                    for (ulong b = 0; b <= 3; b++)
                    {
                        ulong newBase = b << (64 - (m + 1) * 2);
                        ulong kMerVariant = kMerWithHole | newBase;
                        int depth = kMerTableLookup(kMerVariant, kMerTable, kMerSize);
                        if (depth >= minDepth)
                        {
                            if (distinctVariants.Add(kMerVariant))
                            {
                                kMerCounts.Add(depth);
                                kMerVariants.Add(kMerVariant);
                            }
                        }
                    }
                }
        }

        // clean (error correct) a starting read. Read depths are checked against the de-noised kMers, and only culled kMers (depth 0) are corrected. 
        private static bool CleanStartingRead(int kMerSize, Dictionary<ulong, int> kMerTable, List<int> contextLengths, Dictionary<int, Dictionary<ulong, int>> kMerContexts, HashSet<ulong> kMersWithZero, 
                                              Dictionary<string, List<int>> startsOfStartingReads, List<string> selectedReads, Sequence seq, int minDepth, int maxReadLength, ReadStats readStats, int r, List<string> traces, out int lastGoodBase, out bool readChanged)
        {
            int cumulativeChanges = 0;
            bool readClean = true;
            string failureReason = null;
            readChanged = false;
            lastGoodBase = -1;
            bool targetRead = false;

            List<ulong> kMerVariants = new List<ulong>(kMerSize * 4);
            List<int> kMerVariantDepths = new List<int>(kMerSize * 4);
            List<int> kMerContextDepths = new List<int>(kMerSize * 4);
            List<bool> variantAlive = new List<bool>(kMerSize * 4);

            //if (seq.ToString() == "CCTGGCTCAGAACGAACGCTGGCGGCGTGCTGAACACATGCAAGTCGAACGTGA")
            //    Debugger.Break();
            //if (r == 4629)
            //{
            //    targetRead = true;
            //    Debugger.Break();
            //}

            string originalRead = (traces == null) ? null : seq.ToString();

            int previousDepth = Math.Max(readStats.initialGoodDepth, readStats.meanDepth);
            int kMerCount = seq.Length - kMerSize + 1;

            // scan the (possibly updated) read looking for aberrant kMers. As the kMerTable has been denoised, only kMers with depth=0 will be considered for replacement.
            for (int m = 0; m < kMerCount; m++)
            {
                // look for kMers that appear to be suspect - poor kMer depth or no context (when possible)
                ulong kMer;
                Sequence.CondenseMer(seq, m, kMerSize, out kMer);
                int currentDepth = kMerTableLookup(kMer, kMerTable, kMerSize);
                bool kMerHasZero = kMersWithZero.Contains(Canonical(kMer, kMerSize));
                int startsFound = -1;

                //if (kMer == 0xA980978E490659BB)
                //    Debugger.Break();
                //if (currentDepth == 1)
                //    Debugger.Break();

                // look for sudden drops - indicating a reason for checking - or a kMer that's been disappeared - or one that's not found in RC form as well.
                if (currentDepth < RoundDiv(previousDepth, 3) || CloseOrLower(currentDepth, readStats.minDepthAllowed) || kMerHasZero || (startsFound = CountStartsWithinStartingReads(seq, m, startsOfStartingReads, selectedReads, kMerSize)) < 1)
                {
                    // generate alternatives to this dubious kMer
                    if (m == 0)
                    {
                        // check all single-sub variants for the first kMer (and make sure current kMer is included if its non-zero)
                        int deepestInitialVariant = GetDeepestVariantDepth(VariantsWanted.allSingle, kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths, readStats.minDepthAllowed);
                        // and if single-subs don't look promising, try double-subs
                        if (deepestInitialVariant < previousDepth || deepestInitialVariant < readStats.minDepthAllowed)
                            deepestInitialVariant = GetDeepestVariantDepth(VariantsWanted.allDouble, kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths, readStats.minDepthAllowed);
                        // ensure the current kMer is always in the list of variants
                        kMerVariants.Add(kMer);
                        kMerVariantDepths.Add(currentDepth);
                        // reset previousDepth if the first kMer looks to be an error
                        if (kMerVariants.Count >= 1 && RelativelyClose(previousDepth, readStats.minDepthAllowed, deepestInitialVariant))
                            previousDepth = deepestInitialVariant;
                    }
                    else
                    {
                        // just vary the last base for all but first kMer in the read
                        GetNextMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths);
                    }

                    kMerContextDepths.Clear();
                    for (int i = 0; i < kMerVariantDepths.Count; i++)
                        kMerContextDepths.Add(0);

                    // find the 'variant' that is not actually a 'variant'
                    int currentVariantIdx = -1;
                    for (int v = 0; v < kMerVariants.Count; v++)
                    {
                        if (kMerVariants[v] == kMer)
                        {
                            currentVariantIdx = v;
                            break;
                        }
                    }

                    // check for a longest context that will fit into the (possibly truncated) read
                    int lastContextLengthIdx = 4;
                    while (m + kMerSize < contextLengths[lastContextLengthIdx] && lastContextLengthIdx > 0)
                        lastContextLengthIdx--;

                    int viableVariants = CheckVariantViability(seq, m, kMerSize, kMerVariants, kMerVariantDepths, kMerContextDepths, readStats.minDepthAllowed, variantAlive, kMerContexts, contextLengths, lastContextLengthIdx, out failureReason);

                    // never declare the current kMer as unviable at this stage (unless it has a zero count etc.)
                    bool looksBroken = currentDepth < readStats.minDepthAllowed || kMerHasZero || startsFound == 0;
                    if (looksBroken)
                        variantAlive[currentVariantIdx] = false;
                    else
                    {
                        variantAlive[currentVariantIdx] = true;
                        kMerVariantDepths[currentVariantIdx] = currentDepth;
                    }

                    // make sure all variants are actually found in the starting reads (once we've moved away from the start of the read)
                    if (m > 0)
                        for (int v = 0; v < kMerVariants.Count; v++)
                        {
                            if (variantAlive[v])
                            {
                                int startsFoundForVariant = 0;
                                if (v == currentVariantIdx)
                                {
                                    if (startsFound == -1)
                                        startsFoundForVariant = CountStartsWithinStartingReads(seq, m, startsOfStartingReads, selectedReads, kMerSize);
                                    else
                                        startsFoundForVariant = startsFound;
                                }
                                else
                                {
                                    char currentBase = seq.Bases[m + kMerSize - 1];
                                    seq.Bases[m + kMerSize - 1] = kMers.baseToChar[v];
                                    startsFoundForVariant = CountStartsWithinStartingReads(seq, m, startsOfStartingReads, selectedReads, kMerSize);
                                    seq.Bases[m + kMerSize - 1] = currentBase;
                                }

                                if (startsFoundForVariant == 0)
                                    variantAlive[v] = false;
                            }
                        }

                    // count how many downstream kMers match for all the viable variants
                    List<int> followers = new List<int>(4);
                    List<bool> abandonedVariants = new List<bool>(4);
                    int maxFollowers = 0;
                    int maxFollowersIdx = 0;
                    int currentVariantFollowers = 0;
                    //List<List<int>> followingMatchesAt = new List<List<int>>();
                    //List<List<ulong>> followingMatchesFound = new List<List<ulong>>();
                    //bool dumpWanted = false;

                    for (int v = 0; v < kMerVariants.Count; v++)
                    {
                        int matches = 0;
                        bool abandoned = true;
                        //List<int> matchesAt = new List<int>();
                        //List<ulong> matchesFound = new List<ulong>();

                        if (variantAlive[v])
                            matches = CountFollowers(seq, m, kMerCount, readStats.minDepthAllowed, kMerVariants[v], kMerSize, kMerTable, kMersWithZero, 0, 0, out abandoned/*, matchesAt, matchesFound*/);

                        if (matches > maxFollowers)
                        {
                            maxFollowers = matches;
                            maxFollowersIdx = v;
                        }

                        followers.Add(matches);
                        if (matches == 0 && m < kMerCount - 1)
                            abandoned = true;
                        abandonedVariants.Add(abandoned);
                        //followingMatchesAt.Add(matchesAt);
                        //followingMatchesFound.Add(matchesFound);

                        if (v == currentVariantIdx)
                            currentVariantFollowers = matches;
                    }
                    //if (dumpWanted)
                    //    DumpFollowingMatches(kMerVariants, kMerVariantDepths, kMerSize, followingMatches, followingMatchesAt, followingMatchesFound, kMerTable);

                    // scan the remaining variants to find the deepest 
                    int deepestVariantDepth = currentDepth;
                    ulong deepestVariant = kMer;
                    int deepestVariantFollowers = currentVariantFollowers;
                    bool allAbandoned = true;
                    for (int v = 0; v < kMerVariantDepths.Count; v++)
                    {
                        int variantDepth = kMerVariantDepths[v];
                        if (variantAlive[v] && variantDepth > deepestVariantDepth)
                        {
                            deepestVariantDepth = variantDepth;
                            deepestVariant = kMerVariants[v];
                            deepestVariantFollowers = followers[v];
                        }
                        allAbandoned &= abandonedVariants[v]/* && followers[v] < (seq.Length - m + kMerSize)/5*/;
                    }

                    // choose the best alternative - deepest or most following matches
                    ulong bestVariant = deepestVariant;
                    int bestVariantDepth = deepestVariantDepth;
                    int bestVariantFollowers = deepestVariantFollowers;

                    if (maxFollowers > deepestVariantFollowers)
                    {
                        bestVariant = kMerVariants[maxFollowersIdx];
                        bestVariantDepth = kMerVariantDepths[maxFollowersIdx];
                        bestVariantFollowers = maxFollowers;
                    }

                    // correct the kMer if it looks like it's likely to be an error coming from a dominant organism
                    bool replaceCurrent = kMer != bestVariant && bestVariantDepth > 0 && bestVariantFollowers > currentVariantFollowers; // && currentDepth < RoundDiv(deepestVariantDepth, errorRate); 

                    bool keepCurrent = currentVariantFollowers >= bestVariantFollowers || VeryClose(currentVariantFollowers, bestVariantFollowers);

                    if (replaceCurrent && !keepCurrent)
                    {
                        //if (kMer == 0x626FB5A0F3EA6C08)
                        //    Debugger.Break();
                        if (targetRead)
                            Console.WriteLine(r + ": @" + m + " " + kMers.ExpandMer(kMer, kMerSize) + "->" + kMers.ExpandMer(bestVariant, kMerSize));

                        // does it look like the repairs are getting out of control?
                        cumulativeChanges++;
                        if (cumulativeChanges > 2)
                        {
                            readClean = false;
                            break;
                        }

                        // decided to replace the current kMer
                        seq.Replace(m, bestVariant, kMerSize);
                        readChanged = true;

                        IncrementMerTable(m, seq, kMerTable, kMerSize);
                        IncrementContextTables(m, seq, contextLengths, kMerContexts, kMerSize);

                        RefreshStatsForRead(seq, kMerTable, kMerSize, minDepth, readStats);
                        currentDepth = bestVariantDepth;
                    }
                    else
                    {
                        // decided that the kMers was broken and unfixable
                        if (currentDepth < readStats.minDepthAllowed || allAbandoned)
                        {
                            readClean = false;
                            break;
                        }
                        else
                        {
                            // decided the current kMer didn't need repair after all
                            lastGoodBase = m;
                            if (m == 0)
                                readStats.initialGoodDepth = currentDepth;
                            if (cumulativeChanges > 0)
                                cumulativeChanges--;
                        }
                    }

                    previousDepth = currentDepth;
                }
                else
                {
                    lastGoodBase = m;
                    previousDepth = currentDepth;
                    if (cumulativeChanges > 0)
                        cumulativeChanges--;
                    continue;
                }
            }

            if (traces != null && readChanged)
                traces.Add(r + ": cleaned" + " min=" + readStats.minDepthAllowed + " pd=" + previousDepth + " " + originalRead + " --> " + seq.ToString());

            lastGoodBase += kMerSize;
            return readClean;
        }

        private static int CountStartsWithinStartingReads(Sequence seq, int m, Dictionary<string, List<int>> startsOfStartingReads, List<string> selectedReads, int kMerSize)
        {
            string seqKey = seq.ToString(20, 12);
            string seqSoFar = seq.ToString(0, m+kMerSize);

            if (!startsOfStartingReads.ContainsKey(seqKey))
                return 0;

            List<int> seqsForkey = startsOfStartingReads[seqKey];
            int startsFound = 0;
            foreach (int r in seqsForkey)
            {
                string read = selectedReads[r];
                if (read != null && read.StartsWith(seqSoFar, StringComparison.Ordinal))
                    startsFound++;
            }
            return startsFound;
        }

        private static void DumpFollowingMatches(List<ulong> kMerVariants, List<int> kMerVariantDepths, int kMerSize, List<int> followingMatches, List<List<int>> followingMatchesAt, List<List<ulong>> followingMatchesFound, Dictionary<ulong, int> kMerTable)
        {
            for (int i = 0; i < kMerVariants.Count; i++)
            {
                Console.Write(kMers.ExpandMer(kMerVariants[i], kMerSize) + " (" + kMerVariantDepths[i] + ") [" + followingMatches[i] + "] ");
                for (int f = 0; f < followingMatchesAt[i].Count; f++)
                    Console.Write(" " + followingMatchesAt[i][f] + "=" + kMers.ExpandMer(followingMatchesFound[i][f], kMerSize) + "(" + kMerTableLookup(followingMatchesFound[i][f], kMerTable, kMerSize) + ")");
                Console.WriteLine();
            }
        }

        private static void DumpFollowingMatches(ulong kMerVariant, int kMerVariantDepth, int kMerSize, int followingMatches, List<int> followingMatchesAt, List<ulong> followingMatchesFound, Dictionary<ulong, int> kMerTable)
        {
            Console.Write(kMers.ExpandMer(kMerVariant, kMerSize) + " (" + kMerVariantDepth + ") [" + followingMatches + "] ");
            for (int f = 0; f < followingMatchesAt.Count; f++)
                Console.Write(" " + followingMatchesAt[f] + "=" + kMers.ExpandMer(followingMatchesFound[f], kMerSize) + "(" + kMerTableLookup(followingMatchesFound[f], kMerTable, kMerSize) + ")");
            Console.WriteLine();
        }

        private static string DumpHashSet(HashSet<ulong> hs)
        {
            string hss = "";
            foreach (ulong km in hs)
                hss += kMers.ExpandMer(km, 32) + ";";
            return hss;
        }

        private static bool TrimStartingRead(int kMerSize, Sequence seq, int lastGoodBaseFwd, HashSet<ulong> kMersWithZero, int readNo, out bool readTrimmed, List<string> traces)
        {
            int startingLength = seq.Length;

            string priorSeq = traces != null ? seq.ToString() : null;

            seq.Length = lastGoodBaseFwd;

            if (seq.Length >= kMerSize)
            {
                // trim back further if the seq ends in (n,0) kMers
                ulong kMerAtEnd;
                Sequence.CondenseMer(seq, seq.Length - kMerSize, kMerSize, out kMerAtEnd);
                while (kMersWithZero.Contains(Canonical(kMerAtEnd, kMerSize)))
                {
                    seq.Length--;
                    if (seq.Length >= kMerSize)
                        Sequence.CondenseMer(seq, seq.Length - kMerSize, kMerSize, out kMerAtEnd);
                    else
                        break;
                }
            }

            if (traces != null)
                traces.Add(readNo + ": trimmed " + priorSeq + "(" + startingLength + ") --> " + seq.ToString() + "(" + seq.Length + ")");

            readTrimmed = seq.Length < startingLength;

            return seq.Length >= kMerSize;
        }

        // generate final set of starting reads for extension. Reads are cleaned/trimmed. Shorter reads that are included in longer reads are subsumed. 
        private static void GenerateFinalStartingReads(List<int> readsWithStartingPrimers, List<string> selectedHeaders, List<string> selectedReads, ReadStats[] statsForReads, List<int> finalStartingReads, int minReadLength, int maxReadLength,
                                                       int kMerSize, int minDepth, Dictionary<ulong, int> kMerTable, List<int> contextLengths, Dictionary<int, Dictionary<ulong, int>> kMerContexts, HashSet<ulong> kMersWithZero, Dictionary<int, int> readPairs,
                                                       out int extendedShortStartingReads, StreamWriter log)
        {
            extendedShortStartingReads = 0;

            if (log != null)
            {
                log.WriteLine("Initial starting reads");
                for (int i = 0; i < readsWithStartingPrimers.Count; i++)
                {
                    int r = readsWithStartingPrimers[i];
                    string read = selectedReads[r];
                    log.WriteLine(r + "\t" + read.Length + "\t" + read);
                }
            }

            // partition the starting reads on their starts - each partition can be processed in parallel
            Dictionary<string, List<int>> partitionedStartingReads = new Dictionary<string, List<int>>(1000);
            for (int i = 0; i < readsWithStartingPrimers.Count; i++)
            {
                string startOfRead = selectedReads[readsWithStartingPrimers[i]].Substring(0, Math.Min(selectedReads[readsWithStartingPrimers[i]].Length, kMerSize));
                if (!partitionedStartingReads.ContainsKey(startOfRead))
                    partitionedStartingReads.Add(startOfRead, new List<int>());
                partitionedStartingReads[startOfRead].Add(readsWithStartingPrimers[i]);
            }

            // build an index over the starting reads to reduce the size of the linear searches that are about to come
            Dictionary<string, List<int>> startingReadsByStarts = new Dictionary<string, List<int>>(100);
            for (int i = 0; i < readsWithStartingPrimers.Count; i++)
            {
                int r = readsWithStartingPrimers[i];
                string readKey = selectedReads[r].Substring(20, 12);
                if (!startingReadsByStarts.ContainsKey(readKey))
                    startingReadsByStarts.Add(readKey, new List<int>());
                startingReadsByStarts[readKey].Add(r);
            }

            // clean, trim, extend each partition in parallel
            int readsCleanedTotal = 0;
            int readsTrimmedTotal = 0;
            int readsDroppedTotal = 0;
            int readsExtendedTotal = 0;

            // will check later if long-enough read starts are coming from trivial read pairs
            HashSet<int> readsWithPairs = new HashSet<int>(readPairs.Count * 2);
            foreach (KeyValuePair<int, int> kvp in readPairs)
            {
                readsWithPairs.Add(kvp.Key);
                readsWithPairs.Add(kvp.Value);
            }

            //foreach (KeyValuePair<string, List<int>> p in partitionedStartingReads)
            Parallel.ForEach(partitionedStartingReads, p =>
            {
                List<int> partition = p.Value;

                int readsCleaned = 0;
                int readsTrimmed = 0;
                int readsDropped = 0;
                int readsExtended = 0;

                Dictionary<string, string> cleanedTrimmedCache = new Dictionary<string, string>(readsWithStartingPrimers.Count);
                HashSet<int> cleaned = new HashSet<int>();
                HashSet<int> trimmed = new HashSet<int>();
                HashSet<int> dropped =  new HashSet<int>();

                List<string> traces = null;

                if (log != null)
                {
                    traces = new List<string>();
                    traces.Add("Cleaning/trimming starting reads for partition [" + p.Key + "] (" + p.Value.Count + ")");
                }

                //if (p.Key == "AGGCAGCAGTGGGGAATTTTGCGCAATGGGCT")
                //    Debugger.Break();
                for (int pi = 0; pi < partition.Count; pi++)
                {
                    int r = partition[pi];

                    bool noNext = ReadNotExtendable(selectedReads[r], kMerTable, kMerSize, kMersWithZero, statsForReads[r].minDepthAllowed);

                    //if (r == 57203)
                    //    Debugger.Break();

                    // only try to correct starting reads where there are known poor quals or kMers
                    if (selectedHeaders[r].Contains(";pq=", StringComparison.Ordinal) || CloseOrLower(statsForReads[r].minDepthFound, statsForReads[r].minDepthAllowed) || noNext)
                    {
                        //if (r == 4629)
                        //    Debugger.Break();
                        //if (selectedReads[r] == "GCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATAGCTAAGTCTGATGTGAAAACTGGGGGCTC")
                        //    Debugger.Break();

                        if (cleanedTrimmedCache.ContainsKey(selectedReads[r]))
                        {
                            bool readWasChanged = true;
                            bool readWasTrimmed = true;

                            string ctRead = cleanedTrimmedCache[selectedReads[r]];

                            if (ctRead != null)
                            {
                                readWasChanged = !selectedReads[r].StartsWith(ctRead, StringComparison.Ordinal);
                                readWasTrimmed = ctRead.Length != selectedReads[r].Length;
                            }

                            if (readWasChanged)
                            {
                                readsCleaned++;
                                cleaned.Add(r);
                            }
                            if (readWasTrimmed)
                            {
                                readsTrimmed++;
                                trimmed.Add(r);
                            }

                            if (readWasChanged || readWasTrimmed)
                            {
                                selectedReads[r] = ctRead;
                                selectedHeaders[r] = selectedHeaders[r].Replace(";pq=", ";cpq=");
                                if (ctRead == null)
                                    readsDropped++;
                                else
                                    RefreshStatsForRead(ctRead, kMerTable, kMerSize, minDepth, statsForReads[r]);
                            }
                        }
                        else
                        {
                            int lastGoodBaseFwd;
                            bool readWasChanged;
                            bool readWasTrimmed = false;
                            Sequence seq = new Sequence(selectedReads[r]);
                            string startingRead = selectedReads[r];

                            bool readIsClean = CleanStartingRead(kMerSize, kMerTable, contextLengths, kMerContexts, kMersWithZero, startingReadsByStarts, selectedReads, seq, minDepth, maxReadLength, statsForReads[r], r, traces, out lastGoodBaseFwd, out readWasChanged);
                            if (readWasChanged)
                            {
                                readsCleaned++;
                                cleaned.Add(r);
                            }

                            if (!readIsClean)
                                readIsClean = TrimStartingRead(kMerSize, seq, lastGoodBaseFwd, kMersWithZero, r, out readWasTrimmed, traces);
                            if (readWasTrimmed)
                            {
                                readsTrimmed++;
                                trimmed.Add(r);
                            }

                            string cleanedSeq = seq.ToString();

                            if (cleanedSeq.Length >= kMerSize && !CheckExtensionViability(cleanedSeq, kMerTable, kMerSize, cleanedSeq.Length+50, statsForReads[r].minDepthAllowed))
                            {
                                int priorLength = cleanedSeq.Length;
                                cleanedSeq = TrimEndOfRead(r, cleanedSeq, kMerTable, kMerSize, minReadLength, traces);
                                if (cleanedSeq.Length < priorLength)
                                {
                                    if (!readWasTrimmed)
                                    {
                                        readsTrimmed++;
                                        trimmed.Add(r);
                                    }
                                    readWasTrimmed = true;
                                }
                            }

                            if ((readWasChanged || readWasTrimmed) && seq.Length >= minReadLength)
                                RefreshStatsForRead(seq, kMerTable, kMerSize, minDepth, statsForReads[r]);

                            if (readIsClean && cleanedSeq.Length >= minReadLength)
                            {
                                selectedReads[r] = cleanedSeq;
                                selectedHeaders[r] = selectedHeaders[r].Replace(";pq=", ";cpq=");
                            }
                            else
                            {
                                selectedReads[r] = null;
                                dropped.Add(r);
                                readsDropped++;
                            }

                            //if (startingRead == "GCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGAATCATTGGGCGTAAAGGGCGCGTAGGCGGATAGCTAAGTCTGATGTGAAAACTGGGGGCTC")
                            //    Debugger.Break();
                            cleanedTrimmedCache.Add(startingRead, selectedReads[r]);
                        }
                    }
                }

                if (log != null)
                {
                    traces.Add("Cleaned/trimmed starting reads:" + " C=" + cleaned.Count + ", T=" + trimmed.Count + ", D=" + dropped.Count);
                    for (int pi = 0; pi < partition.Count; pi++)
                    {
                        int r = partition[pi];
                        string read = selectedReads[r];
                        traces.Add(r + "\t" + (read == null ? "T\tnull" : (read.Length + "\t" + (cleaned.Contains(r) ? "C" : "") + 
                                                                                                (trimmed.Contains(r) ? "T" : "") + 
                                                                                                (dropped.Contains(r) ? "D" : "") + "\t" + read)));
                    }
                }

                // process each partition - find distinct reads, reduce by finding covering reads and extend shorter reads to distinct covering reads

                // get a distinct set of starting reads in this partition
                Dictionary<string, int> distinctStartingReads = new Dictionary<string, int>(partition.Count);
                Dictionary<string, int> distinctStartingReadsSource = new Dictionary<string, int>(partition.Count);
                for (int pi = 0; pi < partition.Count; pi++)
                {
                    int r = partition[pi];
                    string distinctRead = selectedReads[r];

                    if (distinctRead == null)
                        continue;

                    // don't include reads that have been corrected or trimmed - may be unreliable
                    if (cleaned.Contains(r) || trimmed.Contains(r))
                        continue;
                    // and don't always include trivial reads as these will look better than they are - being duplicated
                    bool trivialRead = false;
                    if (readsWithPairs.Count > 0)
                        trivialRead = !readsWithPairs.Contains(r);
                    if (trivialRead && statsForReads[r].meanDepth > statsForReads[r].minDepthAllowed)
                        trivialRead = false;
                    if (trivialRead)
                        continue;

                    if (!distinctStartingReads.ContainsKey(distinctRead))
                    {
                        distinctStartingReads.Add(selectedReads[r], 0);
                        distinctStartingReadsSource.Add(selectedReads[r], r);
                    }

                    distinctStartingReads[distinctRead]++;
                }

                //if (distinctStartingReads.ContainsKey("AGGCAGCAGTAGGGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAG"))
                //{
                //    Debugger.Break();
                //}

                if (log != null)
                {
                        traces.Add("Distinct starting reads");
                        foreach (KeyValuePair<string, int> kvp in distinctStartingReads)
                        {
                            traces.Add("(" + kvp.Value + ")\t" + "[" + distinctStartingReadsSource[kvp.Key] + "]\t" + kvp.Key);
                        }
                }

                // find the distinct set of long-enough starting read starts in this partition
                Dictionary<string, int> longEnoughReads = new Dictionary<string, int>(partition.Count);
                int longEnoughLength = Math.Max(maxReadLength / 2, 70);
                foreach (KeyValuePair<string, int> kvp in distinctStartingReads)
                {
                    string distinctRead = kvp.Key;
                    if (distinctRead.Length < longEnoughLength)
                        continue;
                    string startOfRead = kvp.Key.Substring(0, longEnoughLength);
                    if (!longEnoughReads.ContainsKey(startOfRead))
                        longEnoughReads.Add(startOfRead, 0);
                    longEnoughReads[startOfRead] += kvp.Value;
                }

                // cull any distinct long-enough that are rare or will fail on simple extension
                int totalLengthEnough = 0;
                foreach (int count in longEnoughReads.Values)
                    totalLengthEnough += count;

                List<string> seqsToCull = new List<string>();
                foreach (KeyValuePair<string, int> kvp in longEnoughReads)
                    if (kvp.Value * 100 / totalLengthEnough < 5 || !CheckExtensionViability(kvp.Key, kMerTable, kMerSize, maxReadLength, minDepth))
                        seqsToCull.Add(kvp.Key);

                foreach (string seq in seqsToCull)
                    longEnoughReads.Remove(seq);

                // go through all of the shorter-than-enough reads and try to extend by finding a suitable long-enough seq
                List<int> extendedStartingReads = new List<int>(partition.Count);
                Dictionary<int, string> extensions = null;
                if (log != null)
                    extensions = new Dictionary<int, string>(partition.Count);
                for (int pi = 0; pi < partition.Count; pi++)
                {
                    int r = partition[pi];
                    string read = selectedReads[r];
                    string extendedRead = null;

                    if (read == null)
                        continue;

                    //if (r == 1257)
                    //    Debugger.Break();

                    if (read.Length < longEnoughLength)
                        extendedRead = ChooseLongEnough(read, longEnoughReads);

                    if (extendedRead != null)
                    {
                        readsExtended++;
                        selectedReads[r] = extendedRead;
                        RefreshStatsForRead(extendedRead, kMerTable, kMerSize, minDepth, statsForReads[r]);
                        if (extensions != null)
                            extensions.Add(r, read + " --> " + extendedRead);
                    }

                    if (selectedReads[r].Length >= minReadLength)
                        extendedStartingReads.Add(r);
                    else
                        readsDropped++;
                }

                lock (finalStartingReads)
                {
                    foreach (int r in extendedStartingReads)
                        finalStartingReads.Add(r);
                    readsCleanedTotal += readsCleaned;
                    readsTrimmedTotal += readsTrimmed;
                    readsDroppedTotal += readsDropped;
                    readsExtendedTotal += readsExtended;
                }

                if (log != null)
                { 
                    traces.Add("Extended short reads");
                    foreach (KeyValuePair<int, string> kvp in extensions)
                        traces.Add(kvp.Key + "(E):\t " + kvp.Value);
                }

                if (log != null)
                    lock (log)
                    {
                        foreach (string line in traces)
                            log.WriteLine(line);
                    }
            //}
            });

            Console.WriteLine("cleaned " + readsCleanedTotal + ", trimmed " + readsTrimmedTotal + ", dropped " + readsDroppedTotal + " poor starting reads");
            extendedShortStartingReads = readsExtendedTotal;
        }

        private static string ChooseLongEnough(string shortRead, Dictionary<string, int> longEnoughReads)
        {
            List<string> possibleExtensions = new List<string>(1000);
            foreach (string seq in longEnoughReads.Keys)
                if (seq.StartsWith(shortRead, StringComparison.Ordinal))
                    for (int i = 0; i < longEnoughReads[seq]; i++)
                        possibleExtensions.Add(seq);

            if (possibleExtensions.Count > 0)
            {
                Random rnd = new Random();
                int roll = rnd.Next(possibleExtensions.Count);
                return possibleExtensions[roll];
            }
            else
                return null;
        }

        // trim a read with no viable following path. Find the depth at the (dubious) end of the read, and cut back until a better depth is found
        private static string TrimEndOfRead(int r, string cleanedSeq, Dictionary<ulong, int> kMerTable, int kMerSize, int minReadLength, List<string> traces)
        {
            ulong[] kMersFromRead = new ulong[1000];
            int[] kMerDepths = new int[1000];
            bool[] kMersValid = new bool[1000];
            int deepestDepth;
            int lowestDepth;
            int kMerCount = GetDepthsForRead(cleanedSeq, kMerTable, kMerSize, ref kMersFromRead, ref kMersValid, ref kMerDepths, out deepestDepth, out lowestDepth);

            int finalDepth = kMerDepths[kMerCount - 1];
            int basesToCull = 0;
            int trimmedLength = cleanedSeq.Length;
            for (int i = kMerCount - 1; i >= 0; i--) 
            {
                if (kMerDepths[i] <= finalDepth && trimmedLength > minReadLength)
                {
                    basesToCull++;
                    trimmedLength--;
                }
                else
                    break;
            }

            if (traces != null)
                traces.Add(r + ": trimmed " + basesToCull + " bases from non-viable end. " + cleanedSeq + " --> " + cleanedSeq.Substring(0, cleanedSeq.Length - basesToCull));

            return cleanedSeq.Substring(0, cleanedSeq.Length-basesToCull);
        }

        // check that there is at least one 'next' kMer possible at the end of this read (no next or not all 0,n)
        private static bool ReadNotExtendable(string read, Dictionary<ulong, int> kMerTable, int kMerSize, HashSet<ulong> kMersWithZero, int minDepthAllowed)
        {
            int nonZero = 0;
            int unbalanced = 0;
            ulong kMerAtEnd;
            kMers.CondenseMer(read, read.Length - kMerSize, kMerSize, out kMerAtEnd);
            ulong lbMask = 0xffffffffffffffff << (32 - kMerSize + 1) * 2;
            ulong nextKMerPrefix = kMerAtEnd << 2;

            for (int b = 0; b < 4; b++)
            {
                ulong nextMerVariant = (nextKMerPrefix & lbMask) | ((ulong)b << (32 - kMerSize) * 2);
                int nextVariantDepth = kMerTableLookup(nextMerVariant, kMerTable, kMerSize);
                if (nextVariantDepth > minDepthAllowed)
                {
                    nonZero++;
                    if (kMersWithZero.Contains(nextMerVariant))
                        unbalanced++;
                }
            }

            return nonZero == 0 || nonZero == unbalanced;
        }

        private static void DecrementMerTable(int m, Sequence seq, Dictionary<ulong, int> kMerTable, int kMerSize)
        {
            int craterLength = Math.Min(kMerSize, (seq.Length - (m + kMerSize)));
            int endOfCraterM = m + craterLength;
            for (int i = m; i < endOfCraterM; i++)
            {
                ulong kMer;
                Sequence.CondenseMer(seq, i, kMerSize, out kMer);
                ulong kMerCanonical = Canonical(kMer, kMerSize);

                //if (kMerCanonical == 0xCCFC64FF167113A0)
                //    Debugger.Break();

                if (kMerTable.ContainsKey(kMerCanonical) && kMerTable[kMerCanonical] > 0)
                    kMerTable[kMerCanonical]--;
            }
        }

        private static void IncrementMerTable(int m, Sequence seq, Dictionary<ulong, int> kMerTable, int kMerSize)
        {
            int craterLength = Math.Min(kMerSize, (seq.Length - (m + kMerSize)));
            int endOfCraterM = m + craterLength;
            for (int i = m; i < endOfCraterM; i++)
            {
                ulong kMer;
                Sequence.CondenseMer(seq, i, kMerSize, out kMer);
                ulong kMerCanonical = Canonical(kMer, kMerSize);

                //if (kMerCanonical == 0x14B549577F410F71)
                //    Debugger.Break();

                if (kMerTable.ContainsKey(kMerCanonical) && kMerTable[kMerCanonical] > 0)
                    kMerTable[kMerCanonical]++;
            }
        }

        private static void DecrementContextTables(int m, Sequence seq, List<int> contextLengths, Dictionary<int, Dictionary<ulong, int>> kMerContexts, int kMerSize)
        {
            int craterLength = Math.Min(kMerSize, (seq.Length - (m + kMerSize)));
            int endOfCraterM = m + craterLength;

            for (int i = m; i < endOfCraterM; i++)
            {
                ulong kMer;
                Sequence.CondenseMer(seq, m, kMerSize, out kMer);

                for (int p = 0; p < contextLengths.Count; p++)
                {
                    int contextLength = contextLengths[p];
                    Dictionary<ulong, int> kMerContextTable = kMerContexts[contextLength];
                    int endOfContextM = m - kMerSize + contextLength;

                    if (endOfContextM + kMerSize > seq.Length)
                        break;

                    // private static ulong HashContext(Sequence seq, int m, int contextLength, int kMerSize)
                    ulong hashedContext = HashContext(seq, m, contextLength, kMerSize);
                    //if (hashedContext == 9380718699553647176)
                    //    Debugger.Break();
                    if (kMerContextTable.ContainsKey(hashedContext) && kMerContextTable[hashedContext] > 0)
                        kMerContextTable[hashedContext]--;
                }
            }
        }

        private static void IncrementContextTables(int m, Sequence seq, List<int> contextLengths, Dictionary<int, Dictionary<ulong, int>> kMerContexts, int kMerSize)
        {
            int craterLength = Math.Min(kMerSize, (seq.Length - (m + kMerSize)));
            int endOfCraterM = m + craterLength;

            for (int i = m; i < endOfCraterM; i++)
            {
                ulong kMer;
                Sequence.CondenseMer(seq, m, kMerSize, out kMer);

                for (int p = 0; p < contextLengths.Count; p++)
                {
                    int contextLength = contextLengths[p];
                    Dictionary<ulong, int> kMerContextTable = kMerContexts[contextLength];
                    int endOfContextM = m - kMerSize + contextLength;

                    if (endOfContextM + kMerSize > seq.Length)
                        break;

                    ulong hashedContext = HashContext(seq, m, contextLength, kMerSize);
                    //if (hashedContext == 13687008140760238416)
                    //    Debugger.Break();
                    if (kMerContextTable.ContainsKey(hashedContext) && kMerContextTable[hashedContext] > 0)
                        kMerContextTable[hashedContext]++;
                }
            }
        }

        private static int RoundDiv(int x, int y)
        {
            return (int)Math.Round((double)x / (double)y);
        }

        private static double HarmonicMean(List<int> depths)
        {
            double invSum = 0.0;
            int depthsCounted = 0;
            for (int b = 0; b < depths.Count; b++)
                if (depths[b] != 0)
                {
                    invSum += 1.0 / (double)depths[b];
                    depthsCounted++;
                }
            double meanDepth = (double)depthsCounted / invSum;
            return meanDepth;
        }

        private static int AverageActiveDepths(List<int> depths)
        {
            int count = 0;
            int sum = 0;
            for (int i = 0; i < depths.Count; i++)
                if (depths[i] != 0)
                { 
                    count++;
                    sum += depths[i];
                }
            return sum/count;
        }

        private static bool DepthsAreClose(List<int> depths)
        {
            int highestNonZeroIdx = 0;
            int highestNonZero = 0;
            int nonZeroCount = 0;
            for (int i = 0; i < depths.Count; i++)
                if (depths[i] != 0)
                {
                    if (depths[i] > highestNonZero)
                    {
                        highestNonZeroIdx = i;
                        highestNonZero = depths[i];
                    }
                    nonZeroCount++;
                }
            if (nonZeroCount == 1)
                return false;
            bool allClose = true;
            for (int i = 0; i < depths.Count; i++)
                if (depths[i] > 0 && i != highestNonZeroIdx)
                    allClose &= (VeryClose(depths[i], highestNonZero) || (highestNonZero - depths[i] == 1));
            return allClose;
        }
        private static bool Close(int currentDepth, int previousDepth)
        {
            if (currentDepth == 0 || previousDepth == 0)
                return false;

            int lowestDepth = Math.Min(currentDepth, previousDepth);
            int highestDepth = Math.Max(currentDepth, previousDepth);

            return (Math.Log10(highestDepth + 10) - Math.Log10(lowestDepth + 10)) <= 0.15;
        }

        private static bool CloserToTarget(int currentDepth, int noise, int target)
        {
            int distToNoise = Math.Abs(currentDepth - noise);
            int distToTarget = Math.Abs(currentDepth - target);

            // move distances to log10 space to get relative closeness
            double distToNoiseLog = Math.Log10(distToNoise);
            double distToTargetLog = Math.Log10(distToTarget);

            return (distToTargetLog < distToNoiseLog);
        }

        private static bool CloseOrHigher(int firstDepth, int secondDepth)
        {
            if (firstDepth == 0 || secondDepth == 0)
                return false;
            if (firstDepth > secondDepth)
                return true;
            return Close(firstDepth, secondDepth);
        }

        private static bool HigherAndNotClose(int firstDepth, int secondDepth)
        {
            if (firstDepth == 0 || secondDepth == 0)
                return false;
            if (firstDepth > secondDepth)
                return !Close(firstDepth, secondDepth);
            return false;
        }

        private static bool CloseOrLower(int firstDepth, int secondDepth)
        {
            if (secondDepth == 0)
                return false;
            if (firstDepth < secondDepth)
                return true;
            return Close(firstDepth, secondDepth);
        }

        private static bool VeryClose(int firstDepth, int secondDepth)
        {
            if (firstDepth == 0 || secondDepth == 0)
                return false;

            int lowestDepth = Math.Min(firstDepth, secondDepth);
            int highestDepth = Math.Max(firstDepth, secondDepth);

            return (Math.Log10(highestDepth + 10) - Math.Log10(lowestDepth + 10)) <= 0.05;
        }

        private static bool VeryCloseOrHigher(int firstDepth, int secondDepth)
        {
            if (firstDepth == 0 || secondDepth == 0)
                return false;

            if (firstDepth > secondDepth)
                return true;

            int lowestDepth = Math.Min(firstDepth, secondDepth);
            int highestDepth = Math.Max(firstDepth, secondDepth);

            return (Math.Log10(highestDepth + 10) - Math.Log10(lowestDepth + 10)) <= 0.05;
        }

        private static bool CloseToNoise(int currentDepth, int noiseDepth)
        {
            if (noiseDepth == 0)
                return currentDepth == 0;
            if (currentDepth <= noiseDepth)
                return true;

            return (Math.Log10(currentDepth + 10) - Math.Log10(noiseDepth + 10)) <= 0.05 || (currentDepth - noiseDepth) <= 1;
        }

        private static bool RelativelyClose(int currentDepth, int previousDepth, int relativeToDepth)
        {
            if (currentDepth == 0 || previousDepth == 0)
                return false;

            int difference = Math.Abs(previousDepth - currentDepth);
            return (double)difference / (double)relativeToDepth < 0.05;
        }

        private static bool Crater(int m, int[] kMerDepths, int kMerCount, int kMerSize, int noiseThreshold, int previousGoodDepth)
        {
            bool foundCrater = false;
            // index of first kMer following the region that would be impacted by this (erroneous) kMer - the potential crater
            int afterIdx = m + kMerSize;
            // expected depth of floor of crater
            int floorDepth = kMerDepths[m];
            // extend end-of-crater as needed (allows for additional errors within the crater)
            while (afterIdx < kMerCount && VeryClose(kMerDepths[afterIdx], floorDepth))
               afterIdx++;
            // how many kMers to consider after the crater
            int postCraterSampleLength = 0;
            if (afterIdx < kMerCount)
                postCraterSampleLength = Math.Min(5, kMerCount - afterIdx);

            // depth after crater (if possible)
            if (postCraterSampleLength > 0)
            {
                double invSum = 0.0;
                int count = 0;
                for (int i = afterIdx; i < afterIdx + postCraterSampleLength; i++)
                {
                    if (i == afterIdx)
                    invSum += 1.0 / (double)kMerDepths[i];
                    count++;
                }
                int postCraterMean = (int)((double)count / invSum);

                // calculate mean of the potential crater region
                double invSumCrater = 0.0;
                int countCrater = 0;
                for (int i = m; i < afterIdx; i++)
                {
                    invSumCrater += 1.0 / (double)kMerDepths[i];
                    countCrater++;
                }
                int craterMean = (int)((double)countCrater / invSumCrater);

                int currentDepth = kMerDepths[m];

                foundCrater = RelativelyClose(craterMean, noiseThreshold, postCraterMean) &&
                                  (CloseOrHigher(postCraterMean, previousGoodDepth)) &&
                                  !CloseOrHigher(currentDepth, postCraterMean);
            }

            return foundCrater;
        }

        // Count how many following kMers in a read have counts above the minDepthForRead threshold. This function is recursive to allow it to progress past erroneous bases (that will be cleaned up soon)
        // Consecutive mismatches and too many mismatches are penalised. Provides RHS read context for decisions about alternative kMers (Denoise & CleanRead)
        private static int CountFollowers(Sequence seq, int m, int lastToCheck, int minDepthForRead, ulong kMer, int kMerSize, Dictionary<ulong, int> kMerTable, HashSet<ulong> kMersWithZero, int consecutiveNoMatches, int totalMisMatches, out bool abandoned/*, List<int> matchesAt, List<ulong> matchesFound*/)
        {
            int matches = 0;
            int bestFollowers = 0;
            abandoned = false;
            //List<int> bestMatchesAt = new List<int>();
            //List<ulong> bestMatchesFound = new List<ulong>();

            ulong lbMask = 0xffffffffffffffff << (32 - kMerSize + 1) * 2;
            ulong nextMer = kMer;

            // check every kMer in the remainder of the read
            for (int nm = m + 1; nm < lastToCheck; nm++)
            {
                int nbi = nm + kMerSize - 1;
                char nextBaseChar = seq.Bases[nbi];
                long nextBaseInt = kMers.BaseCharToInt(nextBaseChar);
                nextMer = nextMer << 2 | (ulong)nextBaseInt;
                int nextDepth = kMerTableLookup(nextMer, kMerTable, kMerSize);

                // acceptable kMer so increment the counter and move on
                if (nextDepth >= minDepthForRead)
                {
                    matches++;
                    consecutiveNoMatches = 0;
                    //matchesAt.Add(nm);
                    //matchesFound.Add(nextMer);
                }
                else
                {
                    // hit a poor kMer
                    consecutiveNoMatches++;
                    totalMisMatches++;
                    // if we're still allowed to recover from this error, generate the 4 possible last-base variants if the error kMer, and see how far each of them can get downstream
                    if (consecutiveNoMatches <= 2 && totalMisMatches <= 3)
                    {
                        ulong nextMerPrefix = nextMer & lbMask;
                        for (int b = 0; b < 4; b++)
                        {
                            ulong nextMerVariant = nextMerPrefix | ((ulong)b << (32 - kMerSize) * 2);
                            int nextVariantDepth = kMerTableLookup(nextMerVariant, kMerTable, kMerSize);

                            // found a viable alternative kMer, so call recursively to get the downstream count - and find the best of the viable alternatives
                            if (nextVariantDepth >= minDepthForRead)
                            {
                                List<int> variantMatchesAt = new List<int>();
                                List<ulong> variantMatchesFound = new List<ulong>();
                                int followers = CountFollowers(seq, nm, lastToCheck, minDepthForRead, nextMerVariant, kMerSize, kMerTable, kMersWithZero, consecutiveNoMatches, totalMisMatches, out abandoned/*, variantMatchesAt, variantMatchesFound*/);
                                if (followers > bestFollowers)
                                {
                                    bestFollowers = followers;
                                    //bestMatchesAt = variantMatchesAt;
                                    //bestMatchesFound = variantMatchesFound;
                                }
                            }
                        }
                    }
                    else
                        // only 'abandon' if we haven't run into an unreliable (n,0) end of read
                        abandoned = nextDepth < minDepthForRead || kMersWithZero.Contains(Canonical(nextMer, kMerSize));

                    // recursive calls will have finished the loop
                    break;
                }
            }

            //foreach (int fm in bestMatchesAt)
            //    matchesAt.Add(fm);
            //foreach (ulong fmk in bestMatchesFound)
            //    matchesFound.Add(fmk);
            return matches + bestFollowers;
        }

        // final step in the extension process - trimming away the FP stub, and the TP. This method also decides whether an extended read is 'good' or should be discarded
        private static int TrimExtendedReads(ConcurrentDictionary<int, string> extendedReads, int fwdPrimerHeadLength, int rvsPrimerHeadLength, int primerCoreLength, HashSet<string>[] terminatingPrimers,
                                             int minLength, int maxLength, Dictionary<int, string> trimmedGoodReads, Dictionary<string, int> distinctGoodReads, Dictionary<int, string> TPsFound, Dictionary<string, int> discards, Dictionary<string, int> discardsR, StreamWriter log)
        {
            // trimmedGoodReads just the fully-extended reads - discarded reads get dereplicated and saved in discards. The int is the record number and the string is the trimmed extended read
            int droppedReadCount = 0;
            int partialForwardpPrimerLength = (fwdPrimerHeadLength + primerCoreLength) - (fwdPrimerHeadLength + primerCoreLength) / 2;

            foreach (KeyValuePair<int, string> kvp in extendedReads)
            {
                int r = kvp.Key;
                string extendedRead = kvp.Value;
                bool foundTerminatingPrimer = false;

                //if (extendedRead == "TAAAGATATTGGAACTTTATATTTTATTTTTGGTGTATGAGCAGGAATAATTGGAACTTCTTTAAGTTTATTAATTCGTACAGAATTAGGTAATCCTGGATCACTAATTGGAGATGATCAAATTTATAATACTATTGTTACAGCTCATGCCTTTATTATAATTTTTTTTATAGTTATGCCAATTATAATTGGAGGATTTGGTAATTGATTAATTCCTTTGATATTAGGAGCCCCTGATATAGCTTTCCCACGTATAAATAA")
                //    Debugger.Break();
                //if (extendedRead.Contains("TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACCGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTCAAAACTATCGGTCTGGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAAATCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGACGACTGTTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"))
                //    Debugger.Break();

                // forward primer at the start has already been partially trimmed, leaving a forwardPrimer/2 stub behind
                // so trim away this primer remnant
                extendedRead = extendedRead.Substring(partialForwardpPrimerLength);

                // check for a terminal primer at the end
                string TPFound = "noTPFound";
                int lastPossibleRPIdx = extendedRead.Length - (rvsPrimerHeadLength + primerCoreLength);
                if (lastPossibleRPIdx >= 0)
                {
                    string possibleRPHead = extendedRead.Substring(lastPossibleRPIdx + primerCoreLength);
                    string possibleRPCore = extendedRead.Substring(lastPossibleRPIdx, primerCoreLength);

                    // and trim away the primer if we found one
                    if (terminatingPrimers[head].Contains(possibleRPHead) && terminatingPrimers[core].Contains(possibleRPCore))
                    {
                        TPFound = possibleRPCore + possibleRPHead;
                        extendedRead = extendedRead.Substring(0, lastPossibleRPIdx);
                        foundTerminatingPrimer = true;
                    }
                }

                // no terminal primer found but the extended read is longer than the minimum, but longer than the max then trim the read to this max length
                if (!foundTerminatingPrimer && minLength > 0 && extendedRead.Length > minLength && extendedRead.Length >= maxLength)
                {
                    extendedRead = extendedRead.Substring(0, maxLength);
                }

                if (foundTerminatingPrimer || (minLength > 0 && extendedRead.Length >= minLength))
                {
                    trimmedGoodReads.Add(r, extendedRead);
                    TPsFound.Add(r, TPFound);

                    if (distinctGoodReads.ContainsKey(extendedRead))
                        distinctGoodReads[extendedRead]++;
                    else
                        distinctGoodReads.Add(extendedRead, 1);
                }
                else
                {
                    if (discards.ContainsKey(extendedRead))
                        discards[extendedRead]++;
                    else
                    {
                        discards.Add(extendedRead, 1);
                        discardsR.Add(extendedRead, r);
                    }
                    droppedReadCount++;

                    if (log != null)
                        log.WriteLine("dropped: " + r + " (" + extendedRead.Length + "b," + (foundTerminatingPrimer ? "TP" : "NTP") + ") " + extendedRead);
                }
            }

            return trimmedGoodReads.Count;
        }

        // Loop over all the starting reads, trying to extend each one in turn. Successfully extended reads are cached and the cache value is re-used if possible (no coin toss was used)
        private static ConcurrentDictionary<int, string> ExtendFromStartingReads(List<string> selectedHeaders, List<string> selectedReads, List<int> startingReads, HashSet<int> startingReadRCed, int minReadLength, int maxLength, int longestReadLength, ReadStats[] readStats, int kMerSize, Dictionary<ulong, int> kMerTable,
                                                                       List<int> contextLengths, Dictionary<int, Dictionary<ulong, int>> kMerContexts, int[] contextLengthStats,
                                                                       HashSet<ulong>[] extensionTerminatingPrimers, int fwdPrimerHeadLength, int rvsPrimerHeadLength, int primerCoreLength, Dictionary<int,int> readPairs, ParallelOptions threads, StreamWriter log)
        {
            // Caches - all shared across all extending threads
            // extended starting reads - int is read idx (from startingReads), string is extended read
            ConcurrentDictionary<int, string> extendedStartingReads = new ConcurrentDictionary<int, string>();
            // read prefixes that have been extended previously - no need to repeat
            ConcurrentDictionary<string, Extension> cachedExtensions = new ConcurrentDictionary<string, Extension>();
            // remember paired reads we've already resolved (all 3 updated together under lock protection)
            Dictionary<string, List<int>> cachedReadsContainingTarget = new Dictionary<string, List<int>>();
            Dictionary<string, List<int>> cachedPairedReads = new Dictionary<string, List<int>>(1000);
            Dictionary<string, List<int>> cachedReadsTargetIdxs = new Dictionary<string, List<int>>(1000);
            // remember full pair coverage 
            ConcurrentDictionary<char[], int> cachedFullPairCoverage = new ConcurrentDictionary<char[], int>();

            int readsExtendedProgress = 0;

            int reportingInterval = startingReads.Count / 10;
            if (reportingInterval == 0)
                reportingInterval = 10;

            Stopwatch extendingSW = null;
            if (log != null)
                extendingSW = new Stopwatch();

            // convert the minReadLength to the max context length that will fit inside such a short read (index)
            int minViableContextIdx = 0;
            for (int p = 0; p < contextLengths.Count; p++)
            {
                if (contextLengths[p] <= longestReadLength / 3)
                    minViableContextIdx = p;
                else
                    break;
            }

            int partialForwardPrimerLength = (fwdPrimerHeadLength + primerCoreLength) - (fwdPrimerHeadLength + primerCoreLength) / 2;
            int maxExtendedLength = Int32.MaxValue;
            if (maxLength != 0)
                maxExtendedLength = maxLength + partialForwardPrimerLength + rvsPrimerHeadLength + primerCoreLength;

            if (log != null)
                log.WriteLine("extending reads");

            // extending each of the 'starting' reads
            Parallel.ForEach(startingReads, threads, r =>
            //foreach (int r in startingReads)
            {
                List<string> logees = null;
                if (log != null)
                    logees = new List<string>();

                string startingRead = selectedReads[r];
                if (startingRead != null)
                {
                    if (logees != null)
                        logees.Add(r + ":\t" + startingRead + " hm=" + readStats[r].meanDepth + " min=" + readStats[r].minDepthAllowed);

                    Interlocked.Increment(ref readsExtendedProgress);
                    if (readsExtendedProgress % reportingInterval == 0)
                        Console.WriteLine("extended " + readsExtendedProgress + "/" + startingReads.Count + " reads");

                    bool breakOnRead = false;

                    //if (startingRead == "AGGCAGCAGTAGGGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGCGATGAAGGTCTTCG" || startingRead == "AGGCAGCAGTAGGGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGATGAAGGTTTTCG" || startingRead == "AGGCAGCAGTAGGGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCG")
                    //{
                    //    //Debugger.Break();
                    //    breakOnRead = true;
                    //}
                    if (r == traceOnReadNo)
                    {
                        breakOnRead = true;
                        //startingRead = "AGGCAGCAGTAGGGAATCTTCGGCAATGGACGAAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCG";
                        //Debugger.Break();
                    }

                    Extension cachedExtension;
                    if (cachedExtensions.TryGetValue(startingRead, out cachedExtension) && !breakOnRead)
                    {
                        //Extension cachedExtension = cachedExtensions[startingRead];
                        string cachedExtendedRead = cachedExtension.extendedSeq.ToString();
                        extendedStartingReads.TryAdd(r, cachedExtendedRead);
                        contextLengthStats[statsCachedResultIdx]++;
                        if (log != null)
                            logees.Add("cached:" + " tp=" + cachedExtension.terminatingPrimerReached + " cost=" + cachedExtension.cost + " md=" + cachedExtension.meanDepth + " avg=" + cachedExtension.avgDepth + " " + cachedExtendedRead);
                    }
                    else
                    {
                        HashSet<ulong> loopTrap = new HashSet<ulong>();
                        int loopTrapLength = contextLengths[0];
                        for (int p = 1; p < contextLengths.Count; p++)
                        {
                            if (contextLengths[p] < startingRead.Length)
                                loopTrapLength = contextLengths[p];
                            else
                                break;
                        }
                        bool terminatingPrimerReached;
                        int levelAtTP;
                        bool coinTossed = false;
                        bool extensionAbandoned = false;
                        int costOfExtension = 0;
                        int meanDepth = 0;
                        int avgDepth = 0;
                        int fewestRecursiveCalls = int.MaxValue;
                        if (log != null)
                        {
                            extendingSW.Restart();
                        }

                        Sequence extendedSeq = ExtendRead(1, new List<int>(), new Sequence(startingRead, Math.Min(maxExtendedLength, 5000)), r, selectedHeaders, selectedReads, startingReadRCed, maxExtendedLength, longestReadLength, kMerSize, kMerTable, contextLengths, minViableContextIdx, kMerContexts, contextLengthStats, -1.0f, -1,
                                                          fewestRecursiveCalls, 0, cachedExtensions, loopTrap, loopTrapLength, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength, readPairs, cachedReadsContainingTarget, cachedReadsTargetIdxs, cachedPairedReads, cachedFullPairCoverage,
                                                          readStats[r], out terminatingPrimerReached, out coinTossed, out extensionAbandoned, out costOfExtension, out meanDepth, out avgDepth, out levelAtTP, breakOnRead, logees);

                        string extendedRead = extendedSeq.ToString();

                        //if (extendedRead.ToString() == "GGGACTGGTTGAACTGTTTATCCTCCTTTATCTTCTAGTATTGCTCATAATGGAGCTTCTGTAGATTTAGCAATTTTTTCTCTTCATCTAGCTGGAATTTCAAGAATCTTAGGAGCAGTAAATTTTATTACTACAATTATTAATATACGATCAACAGGTATTACTTTTG")
                        //    Debugger.Break();

                        extendedStartingReads.TryAdd(r, extendedRead);

                        if (cacheThings && !coinTossed)
                        {
                            //if (startingRead == "GCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGCT")
                            //    Debugger.Break();
                            cachedExtension = new Extension();
                            cachedExtension.extendedSeq = extendedSeq;
                            cachedExtension.terminatingPrimerReached = terminatingPrimerReached;
                            cachedExtension.levelAtTP = levelAtTP;
                            cachedExtension.abandoned = extensionAbandoned;
                            cachedExtension.cost = costOfExtension;
                            cachedExtension.meanDepth = meanDepth;
                            cachedExtension.avgDepth = avgDepth;
                            cachedExtensions.TryAdd(startingRead, cachedExtension);
                        }

                        if (log != null)
                        {
                            extendingSW.Stop();
                            int pqIdx = selectedHeaders[r].IndexOf(";cpq=");
                            string pqPart = "";
                            if (pqIdx > 0)
                            {
                                pqPart = selectedHeaders[r].Substring(pqIdx + ";cpq=".Length);
                                int semiIdx = pqPart.IndexOf(';');
                                if (semiIdx != -1)
                                    pqPart = pqPart.Substring(0, semiIdx);
                                pqPart = " pq=" + pqPart;
                            }

                            logees.Add("extended:" + " tp=" + terminatingPrimerReached + " ct=" + coinTossed + " ab=" + extensionAbandoned + " md=" + meanDepth + " avg=" + avgDepth + pqPart + " (" + extendingSW.ElapsedMilliseconds + "ms) " + extendedRead.ToString());
                        }
                    }
                    
                    if (logees != null)
                    {
                        lock (log)
                        {
                            foreach (string logee in logees)
                                log.WriteLine(logee);
                        }
                    }
                }
            //}
            });

            return extendedStartingReads;
        }

        // The core Kelpie recursive extension algorithm as described in the paper. Takes a starting read and extends it base-at-a-time. Multiple viable paths are handled by recursively exploring each one.
        private static Sequence ExtendRead(int level, List<int> forkPath, Sequence seq, int readNo, List<string> selectedHeaders, List<string> selectedReads, HashSet<int> startingReadRCed, int maxExtendedLength, int longestReadLength, int kMerSize, Dictionary<ulong, int> kMerTable, 
                                           List<int> contextLengths, int minViableContextIdx, Dictionary<int, Dictionary<ulong, int>> kMerContexts, int[] contextLengthStats, double invDepthSum, int depthSum, int fewestRecursiveCalls, int lastAcceptedDepth,
                                           ConcurrentDictionary<string, Extension> cachedExtensions, HashSet<ulong> loopTrap, int loopTrapLength, HashSet<ulong>[] extensionTerminatingPrimers, int rvsPrimerHeadLength, int primerCoreLength,
                                           Dictionary<int, int> readPairs, Dictionary<string, List<int>> cachedReadsContainingTarget, Dictionary<string, List<int>> cachedReadsTargetIdxs, Dictionary<string, List<int>> cachedReadPairs, ConcurrentDictionary<char[], int> cachedFullPairCoverage,
                                           ReadStats readStats, out bool terminatingPrimerReached, out bool tossedCoin, out bool abandoned, out int cost, out int meanDepthExtendedRead, out int avgDepthExtendedRead, out int levelAtTP, bool breakOnRead, List<string> logees)
        {
            //bool traceRead = false;
            bool extending = true;
            terminatingPrimerReached = false;
            tossedCoin = false;
            abandoned = false;
            cost = 0;
            meanDepthExtendedRead = readStats.meanDepth;
            avgDepthExtendedRead = readStats.avgDepth;
            levelAtTP = -1;
            int kMersInRead = seq.Length - kMerSize + 1;
            int basesAdded = 0;
            string failureReason = null;
            List<int> LADs = new List<int>();
            List<int> pathToExtension = new List<int>(forkPath);
            int basesSinceLastFork = 0;

            // undo mean/avg calcs so we can do running means/avgs as bases are added
            if (invDepthSum < 0)
                invDepthSum = (1.0f / (double)readStats.meanDepth) * kMersInRead;
            if (depthSum < 0)
                depthSum = readStats.avgDepth * kMersInRead;

            ulong kMer;
            Sequence.CondenseMer(seq, seq.Length - kMerSize, kMerSize, out kMer);

            // may need to stop immediately if the current kMer has hit a TP. This could happen on a recursive call - even it is unlikely
            if (kMerIsTerminating(kMer, kMerSize, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength))
            {
                terminatingPrimerReached = true;
                meanDepthExtendedRead = (int)((double)kMersInRead / invDepthSum);
                avgDepthExtendedRead = depthSum / kMersInRead;
                levelAtTP = level;
                return seq;
            }

            List<ulong> kMerVariants = new List<ulong>(4);
            List<int> kMerVariantDepths = new List<int>(4);
            List<int> kMerContextDepths = new List<int>(4);
            List<bool> variantAlive = new List<bool>(4);
            List<bool> variantPreviouslyAlive = new List<bool>(4);
            List<int> previousContextDepths = new List<int>(4);
            List<int> currentContextDepths = new List<int>(4);
            List<int> followingContextDepths = new List<int>(4);
            int midContextLength = contextLengths[kMerContexts.Count / 2];

            for (int v = 0; v < 4; v++)
            {
                variantAlive.Add(false);
                variantPreviouslyAlive.Add(false);
                previousContextDepths.Add(0);
                currentContextDepths.Add(0);
                followingContextDepths.Add(0);
                kMerContextDepths.Add(0);
            }

            if (logees != null && breakOnRead)
                logees.Add(new string(' ', level * 4) + level + ": " + "extending " + seq.ToString());

            // keep on adding bases, one at a time, while there is a clear choice
            while (extending)
            {
                // if we've reached the defined maxLength for an extended read, stop extending as well
                if (seq.Length >= maxExtendedLength && maxExtendedLength > 0)
                {
                    if (logees != null)
                        failureReason = "maxLength";
                    break;
                }

                int viableAlternatives = 0;
                int chosenIdx = -1;
                contextLengthStats[statsNextIdx]++;

                //if (kMer == 0x20055E1BB4ABF96B && breakOnRead)
                //    Debugger.Break();
                //if (kMer == 0xAFB0247F49AA2829)
                //    Debugger.Break();
                //if (kMer == 0xFFFFF6F003CFFFF4 && seq.Length == 860)
                //    Debugger.Break();
                //if (breakOnRead && seq.Length == 241)
                //    Debugger.Break();
                //if (breakOnRead && level == 1 && seq.Length == 77)
                //    Debugger.Break();
                //if (kMer == 0xC5CF30FA28FE80F8)
                //{
                //    Sequence.CondenseMer(seq, seq.Length - 2 * kMerSize - 8, kMerSize, out ulong preceding);
                //    if (preceding == 0xBC49D39FF3CC3FFF && seq.Bases[seq.Length - 32 - 2] == 'T')
                //        Debugger.Break();
                //}
                //if (breakOnRead && seq.ToString() == "CCCGGCTCAGAATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAGCGAGAAAGTTCCCTTCGGGGAGCGAGTACAGCGGCGCACGGGTGAGTAACGCGTGGGTAACCTACCCTGGAGTGGGGTATAACTCGCCGAAAGGCGGGCTAATCCCGCATACGTCCTGCAGGCACATGACTGCAGGAGAAAGATGGCCTCTGCTTGCAAGCTATCGCTCTTGGATGGGCCCGCGTTAGATTAGCTTGTTGGTGAGGTAACGGCTCACCAAGGCTACGATCTATAGCTGGTTTGAGAGGACGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGAGGAATATTGCGCAATGGGGGAAACCCTGACGCAGCGACGCCGCGTGAGTGATGAAGGCCTTCGGGTTGTAAAGCTCTGTCAGAGGGAAAGAAAGGGGATGCGGAATAATACCCGCAGCCATTGACGGTACCCTCGGAGGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGTGTTATTCGGAATCACTGGGCGTAAAGAGCACGTAGGCGGATGGGTAAGTCAGATGTGAAATCCCGGGGCTCAACTTCGGAACTGCATTTGATACTGCCCGTCTTGAGTGTGGTAGAGGGGGATGGAATTCCCGGTGTAGAGGTGAAATTCGTAGATATCGGGAGGAACACCAGAGGCGAAGGCGATCCCCTGGGCCATTACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATG")
                //    Debugger.Break();

                basesSinceLastFork++;
                int runningMeanDepth = (int)((double)kMersInRead / invDepthSum);
                int runningAvgDepth = depthSum / kMersInRead;
                int minDepthForRead = lastAcceptedDepth > 0 ? Math.Min(runningMeanDepth / 10, lastAcceptedDepth / 10) : runningMeanDepth / 10;
                if (minDepthForRead < 2)
                    minDepthForRead = 2;

                // generate and sanity check all 4 possible following kMers 
                GetNextMerVariantsAndCounts(kMer << 2, kMerTable, kMerSize, kMerVariants, kMerVariantDepths);
                viableAlternatives = CheckVariantViability(seq, seq.Length - kMerSize + 1, kMerSize, kMerVariants, kMerVariantDepths, kMerContextDepths, minDepthForRead, variantAlive, kMerContexts, contextLengths, minViableContextIdx, out failureReason);

                //if (kMer == 0x0000FFD380000CFF && viableAlternatives == 0)
                //    Debugger.Break();

                // none of the kMers were viable so stop extending
                if (viableAlternatives == 0)
                    break;

                // no decision to make if there's only one choice 
                if (viableAlternatives == 1)
                {
                    chosenIdx = FindChosenVariant(variantAlive);
                    contextLengthStats[statsSingleIdx]++;
                    if (breakOnRead && logees != null)
                    {
                        GetNextMerVariantsAndCounts(kMer << 2, kMerTable, kMerSize, kMerVariants, kMerVariantDepths);
                        logees.Add(new string(' ', level * 4) + level + ": " + "1 viable choice@" + seq.Length + " '" + kMers.ExpandMer(kMerVariants[chosenIdx], kMerSize) + 
                                "' from [" + kMerVariantDepths[0] + "," + kMerVariantDepths[1] + "," + kMerVariantDepths[2] + "," + kMerVariantDepths[3] + "]" + 
                                " rmd=" + runningMeanDepth + " rad=" + runningAvgDepth +
                                " (" + kMerContextDepths[0] + "," + kMerContextDepths[1] + "," + kMerContextDepths[2] + "," + kMerContextDepths[3] + ")");
                    }
                }

                int lastContextLengthChecked = 0;   // debug use 

                //if (breakOnRead && level == 1 && viableAlternatives > 1)
                //    Debugger.Break();

                // no obvious next base so check the alternatives for context support (longer kMers)
                if (viableAlternatives > 1)
                {
                    // demand moderately long contexts - longer contexts are less common so don't be too demanding
                    int minLengthForContext = Math.Min(longestReadLength / 3, seq.Length);
                    if (minLengthForContext < shortestContextLength)
                        minLengthForContext = shortestContextLength;
                    // round this length to a full context boundary
                    int contextIdx = contextLengths.Count - 1;
                    int minContextLength = contextLengths[contextIdx];
                    while (minContextLength > minLengthForContext)
                        minContextLength = contextLengths[--contextIdx];

                    // iterate down through the context lengths, until we run into multiple choices or past the shortest context length
                    int pli = contextLengths.Count - 1;
                    int previousViableAlternatives = viableAlternatives;
                    viableAlternatives = 0;

                    while (pli >= 0)
                    {
                        int contextLength = contextLengths[pli];
                        lastContextLengthChecked = contextLength;

                        int maxContextDepth = GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, pli, kMerContexts, currentContextDepths);

                        // continue if there are no matching contexts and we haven't run out of contexts to check
                        if (maxContextDepth == 0 && contextLength >= minContextLength)
                        {
                            pli--;
                            continue;
                        }

                        // peek ahead because we want non-trivial and increasing counts
                        if (pli - 1 > 0)
                            GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, pli - 1, kMerContexts, followingContextDepths);
                        else
                            for (int i = 0; i < currentContextDepths.Count; i++)
                                followingContextDepths[i] = currentContextDepths[i] == 0 ? 0 : currentContextDepths[i] + 1;

                        // how many active contexts are there (ignore depths of 1 here)
                        int activeVariantsCount = 0;
                        int activeVariantIdx = -1;

                        // look for active contexts - >1 depth and increasing or same but one non-zero and non-trivial
                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            if (variantAlive[b] && currentContextDepths[b] > 1 && followingContextDepths[b] >= currentContextDepths[b])
                            {
                                activeVariantsCount++;
                                activeVariantIdx = b;
                            }
                        }

                        // prepare to check whether the next highest count is also looking viable
                        int secondHighestDepth = 0;
                        int viableFollowers = 0;
                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            if (variantAlive[b] && currentContextDepths[b] > secondHighestDepth && b != activeVariantIdx)
                                secondHighestDepth = currentContextDepths[b];
                            if (variantAlive[b] && currentContextDepths[b] > 0 && followingContextDepths[b] > 1)
                                viableFollowers++;
                        }

                        // other variants are close, so trigger the >1 path
                        if (activeVariantsCount == 1 && ((secondHighestDepth > 0 && viableFollowers > 1) || DepthsAreClose(followingContextDepths)))
                            activeVariantsCount++;

                        // stopped at this context length with just one viable next kMer
                        if (activeVariantsCount == 1)
                        {
                            contextLengthStats[pli + statsContextsIdx]++;
                            chosenIdx = activeVariantIdx;
                            viableAlternatives = 1;

                            if (logees != null && breakOnRead)
                            {
                                string logee = "";
                                logee += new string(' ', level * 4) + level + ": " + contextLength + "-mer context chose '" + kMers.ExpandMer(kMerVariants[chosenIdx], kMerSize) + "' at " + seq.Length + " with [ ";
                                for (int b = 0; (b < kMerVariants.Count); b++)
                                    logee += "(" + currentContextDepths[b] + "," + followingContextDepths[b] + ") ";
                                logee += "] from [";
                                for (int b = 0; (b < kMerVariants.Count); b++)
                                    logee += kMerVariantDepths[b] + " ";
                                logee += "]";
                                logee += " rdcvg=" + CheckReadCoverage(seq, bases[chosenIdx], kMerSize, startsOfReads, level, breakOnRead, logees) + "/" + (seq.Length + 1);
                                logees.Add(logee);
                                //List<int> contextDepths = new List<int>(currentContextDepths);
                                //for (int p = 0; p < kMerContexts.Count; p++)
                                //{
                                //    GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, p, kMerContexts, contextDepths);
                                //    log.Write(contextLengths[p] + "\t");
                                //    foreach (int cd in contextDepths)
                                //        log.Write(cd + "\t");
                                //    log.WriteLine();
                                //}
                            }

                            break;
                        }

                        // multiple viable contexts - looks like we've gone too far - so clean up the 'alive' flags and move to the next stage of variant viability checking
                        if (activeVariantsCount > 1)
                        {
                            if (logees != null && breakOnRead)
                            {
                                string logee = "";
                                GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, pli - 1, kMerContexts, followingContextDepths);
                                logee += (new string(' ', level * 4) + level + ": " + contextLength + "-mer context indecisive for '" + kMers.ExpandMer(kMer, kMerSize) + "' at " + seq.Length + " with [ ");
                                for (int b = 0; (b < kMerVariants.Count); b++)
                                    logee += "(" + currentContextDepths[b] + "," + followingContextDepths[b] + ") ";
                                logee += "]";
                                logees.Add(logee);

                                List<int> contextDepths = new List<int>(currentContextDepths);
                                for (int p = 0; p < kMerContexts.Count; p++)
                                {
                                    if (contextLengths[p] <= seq.Length)
                                    {
                                        string logContext = "";
                                        logContext += new string(' ', level * 4);
                                        GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, p, kMerContexts, contextDepths);
                                        logContext += contextLengths[p] + "\t";
                                        foreach (int cd in contextDepths)
                                            logContext += cd + "\t";
                                        logees.Add(logContext);
                                    }
                                }
                            }
                            viableAlternatives = 0;
                            for (int b = 0; (b < kMerVariants.Count); b++)
                                if (variantAlive[b])
                                {
                                    variantAlive[b] = false;
                                    if (followingContextDepths[b] > 0)
                                    {
                                        viableAlternatives++;
                                        variantAlive[b] = true;
                                    }
                                }
                            break;
                        }

                        // move on to the next context length unless we've hit the lower limit
                        if (contextLength > minContextLength)
                            pli--;
                        else
                            break;
                    }

                    //if (kMer == 0xAFB0247F49AA2829)
                    //  Debugger.Break();
                }

                // contexts check indicated that there were no good variants so abandon this extension attempt
                if (viableAlternatives == 0)
                {
                    meanDepthExtendedRead = runningMeanDepth;
                    avgDepthExtendedRead = runningAvgDepth;
                    failureReason = "no context";
                    if (logees != null && breakOnRead)
                    {
                        List<int> contextDepths = new List<int>(currentContextDepths);
                        for (int p = 0; p < kMerContexts.Count; p++)
                        {
                            if (contextLengths[p] <= seq.Length)
                            {
                                string logContext = "";
                                logContext += new string(' ', level * 4);
                                GetContextDepths(seq, kMerVariants, kMerSize, variantAlive, 1, contextLengths, p, kMerContexts, contextDepths);
                                logContext += contextLengths[p] + (contextLengths[p] == lastContextLengthChecked ? "*" : "") + "\t";
                                foreach (int cd in contextDepths)
                                logContext += cd + "\t";
                                logees.Add(logContext);
                            }
                        }
                    }
                    break;
                }

                if (viableAlternatives > 1)
                {
                    int[] coveragesFromReads = new int[4];
                    int wellCoveredCount = 0;
                    int wellCoveredIdx = 0;
                    string coverageType = "single";

                    //if (breakOnRead)
                    //    Debugger.Break();

                    // first check whether the current seq variant is completely covered by a single read
                    for (int b = 0; b < 4; b++)
                        if (variantAlive[b] && seq.Length < longestReadLength)
                        {
                            coveragesFromReads[b] = CheckFullReadCoverage(seq, bases[b], kMerSize, startsOfReads, level, breakOnRead, logees);
                            if (coveragesFromReads[b] == seq.Length + 1)
                            {
                                wellCoveredCount++;
                                wellCoveredIdx = b;
                            }
                        }

                    // now check if the seq can be covered by multiple reads
                    if (wellCoveredCount != 1)
                    {
                        coverageType = "multiple";
                        for (int b = 0; b < 4; b++)
                            if (variantAlive[b])
                            {
                                coveragesFromReads[b] = CheckReadCoverage(seq, bases[b], kMerSize, startsOfReads, level, breakOnRead, logees);
                                if (coveragesFromReads[b] == seq.Length + 1)
                                {
                                    wellCoveredCount++;
                                    wellCoveredIdx = b;
                                }
                            }
                    }

                    if (wellCoveredCount == 1)
                    {
                        viableAlternatives = 1;
                        chosenIdx = wellCoveredIdx;
                        contextLengthStats[statsReadsCoveredIdx]++;
                        if (logees != null && breakOnRead)
                            logees.Add(new string(' ', level * 4) + level + ": " + coverageType + " read coverage chose '" + kMers.ExpandMer(kMerVariants[chosenIdx], kMerSize) 
                                           + "' at " + seq.Length + " [" + coveragesFromReads[0] + "," + coveragesFromReads[1] + "," + coveragesFromReads[2] + "," + coveragesFromReads[3] + "]");
                    }
                }

                // resolve with backwards paired reads if needed (and available)
                if (readPairs.Count > 0 && viableAlternatives > 1)
                {
                    //if (breakOnRead)
                    //    Debugger.Break();

                    if (seq.Length > longestReadLength)
                        CheckViabilityWithPairedReads(level, seq, selectedHeaders, selectedReads, startingReadRCed, longestReadLength, kMerSize, contextLengthStats, readPairs, cachedReadsContainingTarget, cachedReadsTargetIdxs, cachedReadPairs,
                                                      pathToExtension, kMerVariants, kMerVariantDepths, variantAlive, variantPreviouslyAlive, ref viableAlternatives, ref chosenIdx, breakOnRead, logees);

                    if (breakOnRead && logees != null)
                    {
                        int[] fpc = new int[4];
                        for (int b = 0; b < 4; b++)
                        {
                            if (variantAlive[b])
                                fpc[b] = CheckSeqAgainstPairedReads(new Sequence(seq.ToString() + bases[b]), selectedHeaders, selectedReads, longestReadLength, readPairs, startingReadRCed, cachedFullPairCoverage, breakOnRead, level, logees);
                        }
                        logees.Add(new string(' ', level * 4) + level + ": " + " fpc[" + fpc[0] + "," + fpc[1] + "," + fpc[2] + "," + fpc[3] + "]");
                    }
                }

                // still have multiple alternatives, so see how far downstream each of them can get

                // don't do a further recursive call if we're already at the maximum call depth
                // be generous if we coming from a long run of easy next-base decisions
                int recursiveBoost = basesSinceLastFork / longestReadLength;
                recursiveBoost = Math.Max(recursiveBoost, 2);
                int effectiveLevel = Math.Max(level - recursiveBoost, 1);
                // or if we're almost finished (maxExtendedLength was boosted by 10% for filtering)
                if (maxExtendedLength != int.MaxValue && (seq.Length * 100 / (maxExtendedLength - maxExtendedLength/10)) > 90)
                    effectiveLevel--;
                if (viableAlternatives > 1 && (effectiveLevel >= maxRecursion))
                {
                    viableAlternatives = 0;
                    abandoned = true;
                    break;
                }

                if (viableAlternatives > 1)
                {
                    //if (kMer == 0xAFB0247F49AA2829)
                    //    Debugger.Break();
                    //if (kMerVariants[1] == 0xF127F3C7FCF37F15)
                    //    Debugger.Break();

                    Extension[] variantExtensions = new Extension[4];
                    int tpReachedCount = 0;
                    string reason = null;

                    // remember the path taken to get to this fork
                    pathToExtension.Add(seq.Length - 1);

                    cost = 1;
                    if (basesAdded <= 1)
                        cost = 2;

                    contextLengthStats[statsDownstreamIdx]++;     // count how many times we needed to look downstream to resolve ambiguity

                    if (logees != null && breakOnRead)
                        logees.Add(new string(' ', level * 4) + level + "(" + effectiveLevel + "): trialling @" + seq.Length + " (hm=" + runningMeanDepth + ", avg=" + runningAvgDepth + ", min=" + minDepthForRead + ", blf=" + basesSinceLastFork + ") " + seq.ToString());                   

                    // recursively explore each of the potential paths and choose one of them
                    // if multiple paths can reach a terminal primer at the same lowest cost, choose amongst these alternatives in proportion to their abundance
                    // as all of these paths are equally 'correct'... 

                    // try the deepest path first - and if that gets to a terminal primer, use the length of this extended sequence to set an upper limit on laser extensions.
                    // done to try to tame some aberrant tree explorations that are eventually abandoned anyway.
                    int deepestVariantDepth = 0;
                    int deepestVariantIdx = 0;
                    for (int b = 0; b < 4; b++)
                    {
                        if (variantAlive[b] && kMerVariantDepths[b] > deepestVariantDepth)
                        {
                            deepestVariantDepth = kMerVariantDepths[b];
                            deepestVariantIdx = b;
                        }
                    }

                    if (logees != null && breakOnRead)
                        logees.Add(new string(' ', level * 4) + level + ": " + "append " + bases[deepestVariantIdx] + "@" + seq.Length + " frc=" + fewestRecursiveCalls + " (" + kMerVariantDepths[deepestVariantIdx] + ")");

                    // explore the deepest variant first
                    variantExtensions[deepestVariantIdx] = ExplorePossibleExtension(seq, bases[deepestVariantIdx], selectedHeaders, selectedReads, startingReadRCed, cachedExtensions, maxExtendedLength, longestReadLength, invDepthSum, depthSum, fewestRecursiveCalls, lastAcceptedDepth,
                                                                            kMerSize, kMerTable, contextLengths, minViableContextIdx, kMerContexts,
                                                                            new HashSet<ulong>(loopTrap), loopTrapLength, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength, readPairs, cachedReadsContainingTarget, cachedReadsTargetIdxs, cachedReadPairs, cachedFullPairCoverage,
                                                                            readStats, level, pathToExtension, readNo, kMerVariantDepths[deepestVariantIdx], contextLengthStats, breakOnRead, logees);
                    if (variantExtensions[deepestVariantIdx].terminatingPrimerReached && maxExtendedLength == int.MaxValue)
                        maxExtendedLength = variantExtensions[deepestVariantIdx].extendedSeq.Length + variantExtensions[deepestVariantIdx].extendedSeq.Length / 10;
                    if (variantExtensions[deepestVariantIdx].terminatingPrimerReached)
                        if (variantExtensions[deepestVariantIdx].levelAtTP < fewestRecursiveCalls)
                            fewestRecursiveCalls = variantExtensions[deepestVariantIdx].levelAtTP;

                    // and then the others now we've set maxExtendedLength
                    for (int b = 0; b < 4; b++)
                    {
                        if (variantAlive[b])
                        {
                            if (b != deepestVariantIdx)
                            {
                                if (logees != null && breakOnRead)
                                    logees.Add(new string(' ', level * 4) + level + ": " + "append " + bases[b] + "@" + seq.Length + "(" + " frc=" + fewestRecursiveCalls + " (" + kMerVariantDepths[b] + ")");

                                variantExtensions[b] = ExplorePossibleExtension(seq, bases[b], selectedHeaders, selectedReads, startingReadRCed, cachedExtensions, maxExtendedLength, longestReadLength, invDepthSum, depthSum, fewestRecursiveCalls, lastAcceptedDepth, 
                                                                                kMerSize, kMerTable, contextLengths, minViableContextIdx, kMerContexts,
                                                                                new HashSet<ulong>(loopTrap), loopTrapLength, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength, readPairs, cachedReadsContainingTarget, cachedReadsTargetIdxs, cachedReadPairs, cachedFullPairCoverage,
                                                                                readStats, level, pathToExtension, readNo, kMerVariantDepths[b], contextLengthStats, breakOnRead, logees);
                            }

                            if (variantExtensions[b].terminatingPrimerReached)
                                if (variantExtensions[b].levelAtTP < fewestRecursiveCalls)
                                    fewestRecursiveCalls = variantExtensions[b].levelAtTP;

                            if (variantExtensions[b].terminatingPrimerReached)
                            {
                                tpReachedCount++;
                                chosenIdx = b;
                            }
                        }
                    }

                    // cache any extensions that are safe to cache (tp reached without coin toss)
                    for (int b = 0; b < 4; b++)
                    {
                        if (variantAlive[b])
                        {
                            if (cacheThings && !variantExtensions[b].alreadyCached && !variantExtensions[b].coinTossed)
                            {
                                cachedExtensions.TryAdd(variantExtensions[b].startingRead, variantExtensions[b]);

                                //if (variantExtensions[b].startingRead == "GCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGGGCGTGTAGGTGGTCAAACAAGTCAGGTGTGAAACCCCGGGGCTTAACTTCGGGCATGCATTTGAAACTGTATGACTTGAGGGCAAGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCA")
                                //    Debugger.Break();

                                if (logees != null && breakOnRead)
                                {
                                    logees.Add(new string(' ', level * 4) + level + ": +cache+ " + "[" + bases[b] + "]" + " tp=" + variantExtensions[b].terminatingPrimerReached + " d=" + kMerVariantDepths[b] + " ct=" + variantExtensions[b].coinTossed +
                                                  " ab=" + variantExtensions[b].abandoned + " cost=" + variantExtensions[b].cost + " md=" + variantExtensions[b].meanDepth + " avg=" + variantExtensions[b].avgDepth + " " +
                                                  variantExtensions[b].startingRead + " --> " + variantExtensions[b].extendedSeq.ToString());
                                }
                            }
                        }
                    }

                    // choose amongst the alternatives
                    //      ones that lead to a terminal primer are always chosen if any exist
                    //      if there was only one we would have chosen its index on the way
                    //      if more than one such TP path exists..
                    //          choose the lowest cost extension (if one choice lower)
                    //          otherwise choose among the paths in proportion to their kMer counts
                    //      if no TP paths were found, just choose the longest path - which will probably be rejected as too short later anyway

                    if (tpReachedCount == 1)
                    {
                        if (logees != null)
                            reason = "1tp";
                        contextLengthStats[statsSingleDSTIdx]++;
                    }

                    // mark any no-tp-reached variants as unviable
                    if (tpReachedCount > 1)
                        for (int b = 0; b < kMerVariants.Count; b++)
                            variantAlive[b] = variantExtensions[b] != null && variantExtensions[b].terminatingPrimerReached;

                    // multiple good choices after trialling extensions...
                    //if (log != null && tpReachedCount > 1)
                    //{
                    //    log.Write(new string(' ', level * 4) + level + ": choices @" + seq.Length + " min=" + readStats.minDepth + " avg=" + readStats.avgDepth + " mean=" + readStats.meanDepth + " last=" + lastAcceptedDepth + "\t");
                    //    for (int b = 0; b < kMerVariants.Count; b++)
                    //    {
                    //        if (variantAlive[b] && variantExtensions[b].terminatingPrimerReached)
                    //            log.Write(bases[b] + "(kd=" + kMerVariantDepths[b] + ",m=" + variantExtensions[b].meanDepth + ",cost=" + variantExtensions[b].cost + ") ");
                    //    }
                    //    log.WriteLine();
                    //}

                    // perhaps we can separate the alternatives using coverage of full extended seq using all paired reads
                    if (tpReachedCount > 1 && readPairs.Count > 0)
                    {
                        //if (breakOnRead)
                        //    Debugger.Break();

                        Stopwatch sw = Stopwatch.StartNew();
                        int[] fullCoverages = new int[kMerVariants.Count];
                        for (int b = 0; b < kMerVariants.Count; b++)
                            if (variantAlive[b])
                                fullCoverages[b] = CheckSeqAgainstPairedReads(variantExtensions[b].extendedSeq, selectedHeaders, selectedReads, longestReadLength, readPairs, startingReadRCed, cachedFullPairCoverage, breakOnRead, level, logees);

                        sw.Stop();

                        if (logees != null && breakOnRead)
                        {
                            string logLine = "";
                            logLine += new string(' ', level * 4) + level + ": full pairs " + "(" + sw.ElapsedMilliseconds + "ms) ";
                            for (int b = 0; b < kMerVariants.Count; b++)
                                if (variantAlive[b])
                                    logLine += " " + bases[b] + " " + fullCoverages[b] + "/" + variantExtensions[b].extendedSeq.Length;
                            logees.Add(logLine);
                        }

                        // did only one of the variants get well covered by paired-read matches
                        int previousTPReachedCount = tpReachedCount;
                        tpReachedCount = 0;
                        int bestCoveragePct = 0;
                        int secondBestCoveragePct = 0;
                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            if (variantAlive[b])
                            {
                                int coveragePct = fullCoverages[b] * 100 / variantExtensions[b].extendedSeq.Length;
                                if (coveragePct > bestCoveragePct)
                                {
                                    if (bestCoveragePct > secondBestCoveragePct)
                                        secondBestCoveragePct = bestCoveragePct;
                                    bestCoveragePct = coveragePct;
                                }
                            }
                        }
                        int targetCoveragePct = 100;
                        if (bestCoveragePct < 100)
                            targetCoveragePct = Math.Max(bestCoveragePct, 90);

                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            variantPreviouslyAlive[b] = variantAlive[b];
                            if (variantAlive[b])
                            { 
                                if (fullCoverages[b] * 100 / variantExtensions[b].extendedSeq.Length >= targetCoveragePct)
                                {
                                    tpReachedCount++;
                                    variantAlive[b] = true;
                                    chosenIdx = b;
                                }
                                else
                                    variantAlive[b] = false;
                            }
                        }

                        if (logees != null && tpReachedCount == 1)
                            reason = "fullPR";

                        if (tpReachedCount == 1)
                            contextLengthStats[statsPairedReadsFullIdx]++;

                        // didn't work, restore previous values and choose best available PR cover if there is a single 'better' choice
                        if (tpReachedCount == 0)
                        {
                            // remove any paths with no paired-read support at all
                            for (int b = 0; b < kMerVariants.Count; b++)
                            {
                                variantAlive[b] = variantPreviouslyAlive[b];
                                // and select any alternative which is better than any of the others
                                if (variantPreviouslyAlive[b] && bestCoveragePct > 0 && fullCoverages[b] > 0 && !Close(secondBestCoveragePct, bestCoveragePct))
                                {
                                    variantAlive[b] = true;
                                    chosenIdx = b;
                                    reason = "bestPaired";
                                    tpReachedCount++;
                                }
                                else
                                    variantAlive[b] = false;
                            }

                            // no paired-read help - just restore everything
                            if (tpReachedCount == 0)
                            {
                                for (int b = 0; b < kMerVariants.Count; b++)
                                    variantAlive[b] = variantPreviouslyAlive[b];
                                tpReachedCount = previousTPReachedCount;
                            }
                        }
                    }

                    // ChooseClosest option set... force choice of closest alternative if we're some way down a tree
                    if (tpReachedCount > 1 && level > 1 && chooseClosest)
                    {
                        int closestVariantDistance = int.MaxValue;
                        int closeVariantIdx = 0;
                        for (int b = 0; b < kMerVariants.Count; b++)
                        {
                            if (variantAlive[b] && variantExtensions[b].terminatingPrimerReached)
                            {
                                int distance = Math.Abs(kMerVariantDepths[b] - lastAcceptedDepth);
                                if (distance < closestVariantDistance)
                                {
                                    closestVariantDistance = distance;
                                    closeVariantIdx = b;
                                }
                            }
                        }
                        chosenIdx = closeVariantIdx;
                        tpReachedCount = 1;
                    }

                    // choose the lowest cost alternative if possible - deleted for now - need a better cost calculation
                    //if (tpReachedCount > 1)
                    //{
                    //    int minCost = int.MaxValue;
                    //    int minCostIdx = -1;
                    //    int countAtMinCost = 0;

                    //    for (int b = 0; b < 4; b++)
                    //        if (variantAlive[b] && variantExtensions[b].terminatingPrimerReached)
                    //        {
                    //            if (variantExtensions[b].cost == minCost)
                    //                countAtMinCost++;
                    //            if (variantExtensions[b].cost < minCost)
                    //            {
                    //                minCost = variantExtensions[b].cost;
                    //                minCostIdx = b;
                    //                countAtMinCost = 1;
                    //            }
                    //        }

                    //    if (countAtMinCost == 1)
                    //    {
                    //        chosenIdx = minCostIdx;
                    //        tpReachedCount = 1;

                    //        if (log != null)
                    //            reason = "cost";
                    //        contextLengthStats[statsCostIdx]++;
                    //    }
                    //}

                    // all alternatives have the same cost and all reach a TP, so just toss a coin
                    if (tpReachedCount > 1)
                    {
                        int totalTPCounts = 0;
                        int[] tpCounts = new int[tpReachedCount];
                        int[] tpBases = new int[tpReachedCount];
                        int tpcIdx = 0;
                        // only going to choose amongst TP paths
                        for (int b = 0; b < 4; b++)
                            if (variantAlive[b] && variantExtensions[b].terminatingPrimerReached)
                            {
                                totalTPCounts += kMerVariantDepths[b];
                                // cumulative counts
                                tpCounts[tpcIdx] = kMerVariantDepths[b];
                                if (tpcIdx != 0)
                                    tpCounts[tpcIdx] += tpCounts[tpcIdx - 1];
                                // track which path this TP extension came from
                                tpBases[tpcIdx] = b;
                                tpcIdx++;
                            }

                        // choose one of the paths
                        Random r = new Random();
                        int roll = r.Next(totalTPCounts);
                        Array.Sort<int, int>(tpCounts, tpBases);
                        for (int p = 0; p < tpReachedCount; p++)
                            if (roll < tpCounts[p])
                            {
                                chosenIdx = tpBases[p];
                                //if (breakOnRead && chosenIdx == 3)
                                //    Debugger.Break();
                                break;
                            }
                        tossedCoin = true;
                        contextLengthStats[statsRolledIdx]++;

                        if (logees != null)
                        {
                            reason = "toss[";
                            for (int b = 0; b < 4; b++)
                                if (variantAlive[b] && variantExtensions[b].terminatingPrimerReached)
                                    reason += kMerVariantDepths[b] + ",";
                            reason = reason.Substring(0, reason.Length - 1) + "]";
                        }
                        //if (breakOnRead)
                        //    Console.WriteLine("tossed: " + "[" + chosenIdx + "] " + trialReads[chosenIdx].ToString());
                    }

                    // finally - none of the alternatives reached a TP so choose the longest - not great but this is going to be discarded anyway most likely
                    if (tpReachedCount == 0)
                    {
                        int longestExtension = 0;
                        int longestIdx = 0;

                        for (int b = 0; b < 4; b++)
                        {
                            if (variantAlive[b] && variantExtensions[b].extendedSeq.Length > longestExtension)
                            {
                                longestExtension = variantExtensions[b].extendedSeq.Length;
                                longestIdx = b;
                            }
                        }

                        chosenIdx = longestIdx;
                        contextLengthStats[statsLongestIdx]++;         // chose longest
                        if (logees != null)
                        {
                            reason = "longest[";
                            for (int b = 0; b < 4; b++)
                                if (variantAlive[b])
                                    reason += variantExtensions[b].extendedSeq.Length + ",";
                            reason = reason.Substring(0, reason.Length - 1) + "]";
                        }
                    }

                    //if (kMer == 0x6BB26B80E6C8CDA8 && chosenIdx == 0)
                    //    Debugger.Break();

                    if (logees != null && breakOnRead)
                    {
                        logees.Add(new string(' ', level * 4) + level + "@" + seq.Length + ": " + "explored - chose (" + reason + ") " + bases[chosenIdx] +
                                    " d=" + kMerVariantDepths[chosenIdx] + " ct=" + variantExtensions[chosenIdx].coinTossed + " ab=" + variantExtensions[chosenIdx].abandoned +
                                    " tp=" + variantExtensions[chosenIdx].terminatingPrimerReached + " cost=" + variantExtensions[chosenIdx].cost + " md=" + variantExtensions[chosenIdx].meanDepth + " avg=" + variantExtensions[chosenIdx].avgDepth);
                        if (level == 1)
                            logees.Add("extended (" + variantExtensions[chosenIdx].extendedSeq.Length + ") " + variantExtensions[chosenIdx].extendedSeq.ToString());
                    }

                    // and exit here as the chosen TrialRead will also be the chosen extended read - thanks to the magic of recursion
                    seq = variantExtensions[chosenIdx].extendedSeq;
                    terminatingPrimerReached = variantExtensions[chosenIdx].terminatingPrimerReached;
                    levelAtTP = variantExtensions[chosenIdx].levelAtTP; 
                    tossedCoin = tossedCoin | variantExtensions[chosenIdx].coinTossed;
                    abandoned = variantExtensions[chosenIdx].abandoned; //anyExtensionAbandoned;
                    cost += variantExtensions[chosenIdx].cost;
                    meanDepthExtendedRead = variantExtensions[chosenIdx].meanDepth;
                    avgDepthExtendedRead = variantExtensions[chosenIdx].avgDepth;
                    break;
                }

                // just a single acceptable next base - no ambiguity to resolve - so just append this base and move on... 
                seq.Append(bases[chosenIdx]);
                kMer = kMerVariants[chosenIdx];
                lastAcceptedDepth = kMerVariantDepths[chosenIdx];
                LADs.Add(lastAcceptedDepth);
                kMersInRead++;
                basesAdded++;
                invDepthSum += 1.0 / (double)lastAcceptedDepth;
                depthSum += lastAcceptedDepth;

                // keep a loop trap to detect loops - never seen in normal use... but found in amoeba symbiont repeat!
                // loop trap will contain pairs based on context lengths as these are all known to fit within a read. 
                if (loopTrapLength < contextLengths[contextLengths.Count - 1])
                {
                    int longestPossibleContextLength = loopTrapLength;
                    for (int i = 0; i < contextLengths.Count; i++) 
                    {
                        if (contextLengths[i] < seq.Length && contextLengths[i] > longestPossibleContextLength)
                            longestPossibleContextLength = contextLengths[i];
                    }

                    if (longestPossibleContextLength > loopTrapLength)
                    {
                        loopTrapLength = longestPossibleContextLength;
                        loopTrap.Clear();
                    }
                }
                ulong predecessor;
                Sequence.CondenseMer(seq, seq.Length - loopTrapLength, kMerSize, out predecessor);
                ulong kMerPair = predecessor ^ kMer;
                if (loopTrap.Contains(kMerPair))
                {
                    if (logees != null && breakOnRead)
                        logees.Add(new string(' ', level * 4) + level + ": ER loop (r) with " + kMers.ExpandMer(kMer, kMerSize) + " in " + seq.ToString());
                    terminatingPrimerReached = false;

                    return seq;
                }
                loopTrap.Add(kMerPair);

                //if (breakOnRead && kMer == 0x8FAAE0B6C10ACD6C)
                //    Debugger.Break();

                // have we reached the end (terminal primer)??
                if (kMerIsTerminating(kMer, kMerSize, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength))
                {
                    terminatingPrimerReached = true;
                    levelAtTP = level;
                    meanDepthExtendedRead = runningMeanDepth;
                    avgDepthExtendedRead = runningAvgDepth;
                    break;
                }

            } // while (extending)

            if (logees != null && breakOnRead && failureReason != null)
                logees.Add(new string(' ', level * 4) + level + ": failed @" + seq.Length + " (" + failureReason + ") " + seq.ToString());

            return seq;
        }

        private static int CheckViabilityWithPairedReads(int level, Sequence seq, List<string> selectedHeaders, List<string> selectedReads, HashSet<int> startingReadRCed, 
                                                         int longestReadLength, int kMerSize, int[] contextLengthStats, Dictionary<int, int> readPairs, Dictionary<string, List<int>> cachedReadsContainingTarget, Dictionary<string, List<int>> cachedReadsTargetIdxs, 
                                                         Dictionary<string, List<int>> cachedReadPairs, List<int> pathToExtension, List<ulong> kMerVariants, List<int> kMerVariantDepths, 
                                                         List<bool> variantAlive, List<bool> variantPreviouslyAlive, ref int viableAlternatives, ref int chosenIdx, bool breakOnRead, List<string> logees)
        {
            int[][] coverage = new int[kMerVariants.Count][];
            char[] rcs = new char[seq.Length];

            List<int> matchesAt = null;
            List<int> matchesFor = null;
            List<int> matchesFrom = null;
            if (breakOnRead)
            {
                matchesAt = new List<int>();
                matchesFor = new List<int>();
                matchesFrom = new List<int>();
            }
            int bestCoverage = 0;
            int[] pairedReadVariantCoverage = new int[kMerVariants.Count];
            List<int> readsMatchingSeqContext = new List<int>();
            List<int> readsMatchingSeqIdx = new List<int>();
            List<int> pairedReadsForMatchingContext = new List<int>();

            for (int b = 0; b < kMerVariants.Count; b++)
            {
                if (variantAlive[b])
                {
                    // find the reads with the target (RC form for backwards check) as their pairs might go into the already assembled part of the amplicon
                    // first find all the possible pairs using a small kMer, and we'll then check which if these kMer matches should be kept
                    Span<char> seqSpan = new Span<char>(seq.Bases, 0, seq.Length+1); // Span includes the extra char for the added base being trialled
                    coverage[b] = new int[seq.Length+1];

                    seq.Bases[seq.Length] = bases[b];   // sneak in next base at end of seq
                    // pair target from (extended) seq. As we're looking for RC reads that match to the seq, the target from the seq will be RC to what we're looking for in the reads
                    ReadOnlySpan<char> pairTargetRC = new ReadOnlySpan<char>(seq.Bases, seq.Length - kMerSize+1, kMerSize);
                    string pairTarget = kMers.ReverseComplement(pairTargetRC, rcs);

                    List<int> readsForTarget = null; ;
                    List<int> pairedReadsForTarget = null; ;
                    List<int> readsTargetIdx = null;
                    bool foundInCache = false;
                    lock (cachedReadPairs)
                    {
                        if (cachedReadPairs.ContainsKey(pairTarget))
                        {
                            readsForTarget = cachedReadsContainingTarget[pairTarget];
                            pairedReadsForTarget = cachedReadPairs[pairTarget];
                            readsTargetIdx = cachedReadsTargetIdxs[pairTarget];
                            foundInCache = true;
                        }
                    }

                    if (!foundInCache)
                    {
                        readsForTarget = new List<int>();
                        pairedReadsForTarget = new List<int>();
                        readsTargetIdx = new List<int>();
                        foreach (KeyValuePair<int, int> pair in readPairs)
                        {
                            int R1Idx = pair.Key;
                            int R2Idx = pair.Value;
                            // does the R1 read contain the target
                            string read = selectedReads[R1Idx];
                            if (read != null)
                            {
                                int idx = read.IndexOf(pairTarget, StringComparison.Ordinal);
                                if (idx != -1)
                                {
                                    readsForTarget.Add(R1Idx);
                                    pairedReadsForTarget.Add(R2Idx);
                                    readsTargetIdx.Add(idx);
                                }
                            }
                            // or the R2 from the pair?
                            read = selectedReads[R2Idx];
                            if (read != null)
                            {
                                int idx = read.IndexOf(pairTarget, StringComparison.Ordinal);
                                if (idx != -1)
                                {
                                    readsForTarget.Add(R2Idx);
                                    pairedReadsForTarget.Add(R1Idx);
                                    readsTargetIdx.Add(idx);
                                }
                            }
                        }
                        lock (cachedReadPairs)
                        {
                            if (!cachedReadPairs.ContainsKey(pairTarget))
                            {
                                cachedReadsContainingTarget.Add(pairTarget, readsForTarget);
                                cachedReadPairs.Add(pairTarget, pairedReadsForTarget);
                                cachedReadsTargetIdxs.Add(pairTarget, readsTargetIdx);
                            }
                        }
                    }

                    //if (breakOnRead)
                    //    for (int i = 0; i < readsForTarget.Count; i++)
                    //        logees.Add(i + "\n" + selectedHeaders[readsForTarget[i]] + "\n" + selectedReads[readsForTarget[i]] + " (" + kMers.ReverseComplement(selectedReads[readsForTarget[i]]) + ")" +
                    //                          "\n" + selectedHeaders[pairedReadsForTarget[i]] + "\n" + (selectedReads[pairedReadsForTarget[i]] == null ? "null" : selectedReads[pairedReadsForTarget[i]] + " (" + kMers.ReverseComplement(selectedReads[pairedReadsForTarget[i]]) + ")"));

                    // may well have matched reads that don't fit into this growing contig, so check the 'matched' read actually matches the whole seq, not just the target
                    readsMatchingSeqContext.Clear();
                    pairedReadsForMatchingContext.Clear();

                    for (int i = 0; i < readsForTarget.Count; i++)
                    {
                        // target may be in the middle of the read so only check the read up to the target
                        ReadOnlySpan<char> read = selectedReads[readsForTarget[i]].AsSpan();
                        ReadOnlySpan<char> readUpToTarget = read.Slice(readsTargetIdx[i], selectedReads[readsForTarget[i]].Length - readsTargetIdx[i]);
                        ReadOnlySpan<char> readUpToTargetRC = kMers.ReverseComplementSpan(readUpToTarget, 0, readUpToTarget.Length, rcs);

                        int readMatchIdx = seqSpan.IndexOf(readUpToTargetRC);
                        if (selectedReads[pairedReadsForTarget[i]] != null && readMatchIdx >= 0)
                        {
                            readsMatchingSeqContext.Add(readsForTarget[i]);
                            readsMatchingSeqIdx.Add(readMatchIdx);
                            pairedReadsForMatchingContext.Add(pairedReadsForTarget[i]);
                        }
                    }

                    //if (breakOnRead)
                    //    for (int i = 0; i < readsMatchingSeqContext.Count; i++)
                    //        logees.Add(readsMatchingSeqContext[i] + "\n" + selectedHeaders[readsMatchingSeqContext[i]] + "\n" + selectedReads[readsMatchingSeqContext[i]] + " (" + kMers.ReverseComplement(selectedReads[readsMatchingSeqContext[i]]) + ")" +
                    //                          "\n" + pairedReadsForMatchingContext[i] + "\n" + selectedHeaders[pairedReadsForMatchingContext[i]] + "\n" + (selectedReads[pairedReadsForMatchingContext[i]] == null ? "null" : selectedReads[pairedReadsForMatchingContext[i]] + " (" + kMers.ReverseComplement(selectedReads[pairedReadsForMatchingContext[i]]) + ")"));

                    if (breakOnRead)
                    {
                        matchesAt.Clear();
                        matchesFrom.Clear();
                        matchesFor.Clear();
                    }

                    for (int p = 0; p < pairedReadsForMatchingContext.Count; p++)
                    {
                        int pairIdx = pairedReadsForMatchingContext[p];
                    
                        // get target length from start of read and see if it matches
                        string pairedRead = selectedReads[pairIdx];
                        if (pairedRead == null || pairedRead.Length < kMerSize)
                            continue;
                        // FP' reads were already RCed 
                        if (startingReadRCed.Contains(pairIdx))
                            pairedRead = kMers.ReverseComplement(pairedRead, rcs);
                        string pairedReadTargetStart = pairedRead.Substring(0, kMerSize);
                        int startTargetFoundAt = seqSpan.IndexOf(pairedReadTargetStart);

                        // we know the read matched, so see if its pair matched and record the match if they both did. 
                        if (startTargetFoundAt >= 0)
                        {
                            int basesMatched = 0;
                            int basesPossible = 0;
                            for (int i = 0; i < pairedRead.Length; i++)
                                if (startTargetFoundAt + i < seqSpan.Length)
                                {
                                    basesPossible++;
                                    if (seqSpan[startTargetFoundAt + i] == pairedRead[i])
                                        basesMatched++;
                                }

                            if (breakOnRead)
                            {
                                matchesAt.Add(startTargetFoundAt);
                                matchesFor.Add(basesMatched);
                                matchesFrom.Add(pairedRead.Length);
                            }

                            // exact match so record where they matched in the coverage
                            if (basesMatched == basesPossible)
                            {
                                string read = selectedReads[readsMatchingSeqContext[p]];
                                int readFoundAt = readsMatchingSeqIdx[p];
                                int lengthMatched = Math.Min(seq.Length + 1 - readFoundAt, read.Length);
                                for (int i = 0; i < lengthMatched; i++)
                                {
                                    int locInSeq = readFoundAt + i;
                                    if (seqSpan[locInSeq] == kMers.baseToComplement[read[read.Length-i-1] & 0x1f])
                                        coverage[b][locInSeq]++;
                                    else
                                        break;
                                }
                                for (int i = 0; i < basesPossible; i++)
                                {
                                    int locInSeq = startTargetFoundAt + i;
                                    if (locInSeq < coverage[b].Length)
                                    {
                                        if (pairedRead[i] == seqSpan[locInSeq])
                                            coverage[b][locInSeq]++;
                                        else
                                            break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // finally count the final coverage for each variant
             for (int b = 0; b < kMerVariants.Count; b++)
                if (coverage[b] != null)
                {
                    // check coverage from start to end of seq 
                    for (int i = 0; i < coverage[b].Length; i++)
                        if (coverage[b][i] >= 1)
                            pairedReadVariantCoverage[b]++;
                }

            for (int b = 0; b < kMerVariants.Count; b++)
                if (pairedReadVariantCoverage[b] > bestCoverage)
                    bestCoverage = pairedReadVariantCoverage[b];

            // did only one of the variants get paired-read matches
            int previousViableAlternatives = viableAlternatives;

            if (bestCoverage == seq.Length+1)
            {
                viableAlternatives = 0;

                for (int b = 0; b < kMerVariants.Count; b++)
                {
                    variantPreviouslyAlive[b] = variantAlive[b];
                    if (Close(pairedReadVariantCoverage[b], bestCoverage))
                    {
                        viableAlternatives++;
                        variantAlive[b] = true;
                        chosenIdx = b;
                    }
                    else
                        variantAlive[b] = false;
                }
            }

            string coverages = null;
            if (logees != null && breakOnRead)
            {
                coverages = "[ ";
                for (int b = 0; b < kMerVariants.Count; b++)
                    coverages += pairedReadVariantCoverage[b] + " ";
                coverages += "]";
            }

            if (viableAlternatives == 1)
            {
                contextLengthStats[statsPairedReadsRIdx]++;
                pathToExtension.Add(seq.Length);

                if (logees != null && breakOnRead)
                    logees.Add(new string(' ', level * 4) + level + ": " + "<-- paired reads chose '" + kMers.ExpandMer(kMerVariants[chosenIdx], kMerSize) + "' at " + seq.Length + " with " + coverages + " coverage");
            }

            if (viableAlternatives == 0)
            {
                viableAlternatives = previousViableAlternatives;
                for (int b = 0; b < variantAlive.Count; b++)
                    variantAlive[b] = variantPreviouslyAlive[b];
            }

            if (logees != null && viableAlternatives != 1 && breakOnRead)
                logees.Add(new string(' ', level * 4) + level + ": " + "<-- paired reads indecisive at " + seq.Length + " with " + coverages + " coverage ");

            return viableAlternatives;
        }

        private static int CheckSeqAgainstPairedReads(Sequence seq, List<string> selectedHeaders, List<string> selectedReads, int longestReadLength, Dictionary<int, int> readPairs, HashSet<int> startingReadRCed, ConcurrentDictionary<char[], int> cachedFullPairCoverage, bool breakOnRead, int level, List<string> logees)
        {
            int cachedPairCoverage = 0;
            if (cachedFullPairCoverage.TryGetValue(seq.Bases, out cachedPairCoverage))
                return cachedFullPairCoverage[seq.Bases];

            if (seq.Length < pairedReadKML)
                return 0;

            bool targetFound = false;
            //if (breakOnRead && seq.ToString() == "AGGCAGCAGTGGGGAATCTTGCGCAATGCGCGAAAGCGTGACGCAGCAACGCCGCGTGGGGGAAGAAGGCCTTCGGGTTGTAAACCCCTTTCAGTTGGGACGAAGCTTCGCCGGTGAATAGCCGGCTGGAGTGACGGTACCTTCACAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCTAGCGTTGTCCGGAATCATTGGGCGTAAAGCGCGTGTAGGCGGCCCGGTAAGTCCGCTGTGAAAGTCGGGGGCTCAACCCTCGGATGCCGGTGGATACTGTCGGGCTTGAGTACGGAAGAGGCGAGTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCAATTGCGAAGGCAGCTCGCTGGGACGTTACTGACGCTGAGACGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCC")
            //{
            //    targetFound = true;
            //    //Debugger.Break();
            //}

            int[] coverage = new int[seq.Length];
            char[] rcs = new char[longestReadLength];
            ReadOnlySpan<char> seqSpan = new ReadOnlySpan<char>(seq.Bases, 0, seq.Length);
            Sequence rcSeq = new Sequence(seq);
            Sequence.ReverseComplement(rcSeq);
            ReadOnlySpan<char> rcSeqSpan = new ReadOnlySpan<char>(rcSeq.Bases, 0, rcSeq.Length);
            bool firstRead = true;

            foreach (KeyValuePair<int, int> kvp in readPairs)
            {
                int R1idx = kvp.Key;
                int R2idx = kvp.Value;

                string R1 = selectedReads[R1idx];
                string R2 = selectedReads[R2idx];

                if (R1 == null || R2 == null)
                    continue;

                // short fragment with FP on both reads
                if ((startingReadRCed.Contains(R1idx) && selectedHeaders[R2idx].EndsWith(";FP")) || (startingReadRCed.Contains(R2idx) && selectedHeaders[R1idx].EndsWith(";FP")))
                    continue;

                // starting reads can be RCed so reverse them again in this case
                if (startingReadRCed.Contains(R1idx))
                    R1 = kMers.ReverseComplement(R1, 0, R1.Length, rcs);
                if (startingReadRCed.Contains(R2idx))
                    R2 = kMers.ReverseComplement(R2, 0, R2.Length, rcs);

                ReadOnlySpan<char> R1Span = R1.AsSpan();
                ReadOnlySpan<char> R2Span = R2.AsSpan();

                int minR1MatchLength = Math.Min(kMerSize, R1.Length);
                int minR2MatchLength = Math.Min(kMerSize, R2.Length);

                ReadOnlySpan<char> startR1Span = R1Span.Slice(0, minR1MatchLength);
                ReadOnlySpan<char> startR2Span = R2Span.Slice(0, minR2MatchLength);
                int R1StartAt = seqSpan.IndexOf(startR1Span);
                int R2StartAt = seqSpan.IndexOf(startR2Span);

                // both matched in forward direction - just ignore - good pairs are one forward, one reverse
                if (R1StartAt >= 0 && R2StartAt >= 0)
                    continue;

                // neither matched in forward direction - move on to next pair
                if (R1StartAt == -1 && R2StartAt == -1)
                    continue;

                string fwdRead;
                string rvsRead;
                int minFwdMatchLength = 0;
                int minRvsMatchLength = 0;
                int fwdStartAt = 0;
                int fwdIdx;
                int rvsIdx;
                ReadOnlySpan<char> rvsSpan;

                if (R1StartAt >= 0)     // R1 matches in forward, so R2 must match in reverse
                {
                    fwdRead = R1;
                    minFwdMatchLength = minR1MatchLength;
                    rvsRead = R2;
                    minRvsMatchLength = minR2MatchLength;
                    fwdStartAt = R1StartAt;
                    fwdIdx = R1idx;
                    rvsIdx = R2idx;
                    rvsSpan = R2Span;
                }
                else
                {
                    fwdRead = R2;
                    minFwdMatchLength = minR2MatchLength;
                    rvsRead = R1;
                    minRvsMatchLength = minR1MatchLength;
                    fwdStartAt = R2StartAt;
                    fwdIdx = R2idx;
                    rvsIdx = R1idx;
                    rvsSpan = R1Span;
                }

                ReadOnlySpan<char> endRvsSpan = rvsSpan.Slice(rvsRead.Length - minRvsMatchLength, minRvsMatchLength); // will be matched in reverse so this is the RC start
                int rvsRCStartAt = rcSeqSpan.IndexOf(endRvsSpan);

                // no match on the rvs read - time for next pair
                if (rvsRCStartAt == -1)
                    continue;

                int fwdMatches = 0;
                int rvsMatches = 0;
                int possibleFwdMatches = Math.Min(seq.Length - fwdStartAt, fwdRead.Length);
                int possibleRvsMatches = Math.Min(rvsRCStartAt + kMerSize, rvsRead.Length);

                for (int i = 0; i < possibleFwdMatches; i++)
                    if (seq.Bases[fwdStartAt + i] == fwdRead[i])
                        fwdMatches++;
                for (int i = 0; i < possibleRvsMatches; i++)
                    if (rcSeq.Bases[rvsRCStartAt + kMerSize - 1 - i] == rvsRead[rvsRead.Length - 1 - i])
                        rvsMatches++;

                if (targetFound && logees != null)
                {
                    if (fwdMatches == possibleFwdMatches && rvsMatches == possibleRvsMatches)
                    {
                        if (firstRead)
                        {
                            firstRead = false;
                            logees.Add(seq.ToString());
                        }
                        logees.Add(new string(' ', fwdStartAt) + fwdRead.Substring(0, possibleFwdMatches) + (fwdRead.Length > possibleFwdMatches ? fwdRead.Substring(possibleFwdMatches).ToLower() : "") + " F (" + fwdIdx + ") " + fwdMatches + "=" + (fwdMatches * 100 / possibleFwdMatches) + "%");
                        string fwdRvs = kMers.ReverseComplement(rvsRead);
                        logees.Add(new string(' ', seq.Length - rvsRCStartAt - kMerSize) + fwdRvs.Substring(0, possibleRvsMatches) + (rvsRead.Length > possibleRvsMatches ? fwdRvs.Substring(possibleRvsMatches).ToLower() : "") + " R (" + rvsIdx + ") " + rvsMatches + "=" + (rvsMatches * 100 / possibleRvsMatches) + "%");
                    }
                }

                if (fwdMatches * 100 / possibleFwdMatches == 100 && rvsMatches * 100 / possibleRvsMatches == 100)
                {
                    for (int i = 0; i < possibleFwdMatches; i++)
                        if (seq.Bases[fwdStartAt + i] == fwdRead[i])
                            coverage[fwdStartAt + i]++;
                    int rvsStartInSeq = seq.Length - (rvsRCStartAt + kMerSize);
                    for (int i = 0; i < possibleRvsMatches; i++)
                        if (rcSeq.Bases[rvsRCStartAt + kMerSize - 1 - i] == rvsRead[rvsRead.Length - 1 - i])
                            coverage[rvsStartInSeq + i]++;
                }
            }

            int basesCovered = 0;
            for (int i = 0; i < seq.Length; i++)
                if (coverage[i] > 0)
                    basesCovered++;

            cachedFullPairCoverage.TryAdd(seq.Bases, basesCovered);
            //if (logees != null)
            //    logees.Add(new string(' ', level * 4) + level + ": +covered+ " + basesCovered + " " + seq);
            return basesCovered;
        }

        private static bool TrivialPair(string R1, string R2, int longestReadLength)
        {
            // make r1 the shorter seq
            string r1 = R1;
            string r2 = R2;
            if (R2.Length < R1.Length)
            {
                r1 = R2;
                r2 = R1;
            }
            string r1RC = kMers.ReverseComplement(r1);
            int coreLength = longestReadLength * 2 / 3;
            if (r1.Length - 20 < coreLength) 
                coreLength = r1.Length - 20;
            int coreStart = r1.Length / 2 - coreLength / 2;
            string r1RCCore = r1RC.Substring(coreStart, coreLength);

            return r2.Contains(r1RCCore, StringComparison.Ordinal);
        }

        // checks coverage of (extended) seq with reads (unpaired). Coverage can come from multiple reads
        private static int CheckReadCoverage(Sequence seq, char trialBase, int kMerSize, Dictionary<ulong, List<string>> starts, int level, bool breakOnRead, List<string> logees)
        {
            ulong[] kMersInSeq = new ulong[seq.Length+1];
            bool[] kMersValid = new bool[seq.Length+1];

            //if (breakOnRead)
            //    Console.WriteLine(seq);

            int kMerCount = kMers.GenerateMersFromRead(seq, trialBase, kMerSize, ref kMersInSeq, ref kMersValid);
            int[] coverage = new int[seq.Length+1];
            int[] halves = new int[2];
            //bool firstRead = true;

            for (int i = 0; i < kMerCount; i++)
            {
                // get reads starting with this kMer (either as-is or RCed)
                if (!startsOfReads.ContainsKey(kMersInSeq[i]))
                    continue;
                
                List<string> readsToCheck = startsOfReads[kMersInSeq[i]];

                // map each of these reads onto the seq
                foreach (string read in readsToCheck)
                {
                    // only interested in fully matching reads
                    int matchesForRead = 0;
                    int possibleMatches = 0;
                    for (int b = 0; b < read.Length; b++)
                    {
                        if (b + i >= seq.Length+1)
                            break;

                        if (seq.Bases[i+b] == read[b])
                            matchesForRead++;
                        possibleMatches++;
                    }

                    // found a matching read, add it to the coverage
                    if (matchesForRead == possibleMatches)
                    {
                        //if (breakOnRead && log != null)
                        //{
                        //    if (firstRead)
                        //    {
                        //        firstRead = false;
                        //        log.WriteLine();
                        //    }
                        //    log.WriteLine(new string(' ', i) + read.Substring(0, possibleMatches) + (read.Length > possibleMatches ? read.Substring(possibleMatches).ToLower() : ""));
                        //}

                        for (int b = 0; b < read.Length; b++)
                        {
                            if (b + i >= seq.Length+1)
                                break;
                            coverage[i + b]++;
                        }

                        // remember which half the matching read started in
                        int halfIdx = i / (kMerCount / 2);
                        if (halfIdx == 2)
                            halfIdx = 1;
                        halves[halfIdx]++;
                        // and if the read completely covers the seq, allow it to be counted for both halves
                        if (matchesForRead == seq.Length+1)
                            halves[halfIdx == 0 ? 1 : 0]++;
                    }
                }
            }

            int basesCovered = 0;
            for (int i = 0; i < seq.Length+1; i++) 
                if (coverage[i] > 0)
                    basesCovered++;
            int halvesWithStarts = 0;
            for (int i = 0; i < 2; i++)
                if (halves[i] > 0)
                    halvesWithStarts++;

            return basesCovered * halvesWithStarts / 2;
        }

        // check if seq is fully covered by single reads. Trimmed matching reads have to start at the start of seq an go up to or past the end of seq 
        private static int CheckFullReadCoverage(Sequence seq, char trialBase, int kMerSize, Dictionary<ulong, List<string>> starts, int level, bool breakOnRead, List<string> logees)
        {
            // get reads starting with this initial kMer (either as-is or RCed)
            Span<char> seqSpan = new Span<char>(seq.Bases, 0, seq.Length + 1);
            seqSpan[seq.Length] = trialBase;
            ulong kMerAtStart;
            kMers.CondenseMer(seqSpan, 0, kMerSize, out kMerAtStart);
            if (!startsOfReads.ContainsKey(kMerAtStart))
                return 0; 

            List<string> readsToCheck = startsOfReads[kMerAtStart];
            int matchesForRead = 0;
            int possibleMatches = seq.Length + 1;
            int highestMatches = 0;

            // map each of these reads onto the seq
            foreach (string read in readsToCheck)
            {
                // only interested in fully matching reads
                matchesForRead = 0;
                for (int b = 0; b < read.Length; b++)
                {
                    if (b >= seq.Length + 1)
                        break;

                    if (seq.Bases[b] == read[b])
                        matchesForRead++;
                }

                if (matchesForRead > highestMatches)
                    highestMatches = matchesForRead;
                // found a fully matching read - don't bother looking any further
                if (matchesForRead == possibleMatches)
                    break;
            }

            //if (breakOnRead && matchesForRead == possibleMatches)
            //{
            //    Console.WriteLine(seqSpan.ToString());
            //    foreach (string read in readsToCheck)
            //        Console.WriteLine(read);
            //}
            
            return highestMatches;
        }

        private static Extension ExplorePossibleExtension(Sequence seq, char nextBase, List<string> selectedHeaders, List<string> selectedReads, HashSet<int> startingReadRCed, ConcurrentDictionary<string, Extension> cachedExtensions, int maxExtendedLength, int longestReadLength, double invDepthSum, int depthSum, 
                                                          int fewestRecursiveCalls, int lastAcceptedDepth, int kMerSize, Dictionary<ulong, int> kMerTable, List<int> contextLengths, int minViableContextIdx, Dictionary<int, Dictionary<ulong, int>> kMerContexts,
                                                          HashSet<ulong> loopTrap, int loopTrapLength, HashSet<ulong>[] extensionTerminatingPrimers, int rvsPrimerHeadLength, int primerCoreLength, Dictionary<int, int> readPairs,
                                                          Dictionary<string, List<int>> cachedReadsContainingTarget, Dictionary<string, List<int>> cachedReadsTargetIdxs, Dictionary<string, List<int>> cachedReadPairs, ConcurrentDictionary<char[], int> cachedFullPairCoverage,
                                                          ReadStats readStats, int level, List<int> forkPath, int readNo, int depth,  int[] contextLengthStats, bool breakOnRead, List<string> logees)
        {
            Extension extension;

            Sequence startingSeq = new Sequence(seq, Math.Min(maxExtendedLength, 5000));
            startingSeq.Append(nextBase);
            string startingRead = startingSeq.ToString();

            if (cachedExtensions.TryGetValue(startingRead, out extension) && !breakOnRead)
            {
                //extension = cachedExtensions[startingRead];
                extension.alreadyCached = true;
                contextLengthStats[statsCachedBranchIdx]++;
                if (logees != null && breakOnRead)
                    logees.Add(new string(' ', level * 4) + level + ": (cache@" + extension.cachedAtRead + ") " + "[" + nextBase + "] " + " d=" + depth + " tp=" + extension.terminatingPrimerReached +
                                                " cost=" + extension.cost + " md=" + extension.meanDepth + " avg=" + extension.avgDepth + " " + startingRead + " --> " + extension.extendedSeq.ToString());
            }
            else
            {
                extension = new Extension();
                extension.alreadyCached = false;
                extension.startingRead = startingRead;

                extension.extendedSeq = ExtendRead(level + 1, forkPath, startingSeq, readNo, selectedHeaders, selectedReads, startingReadRCed, maxExtendedLength, longestReadLength, kMerSize, kMerTable, contextLengths, minViableContextIdx, kMerContexts, contextLengthStats, invDepthSum, depthSum,
                                                   fewestRecursiveCalls, lastAcceptedDepth, cachedExtensions, loopTrap, loopTrapLength, extensionTerminatingPrimers, rvsPrimerHeadLength, primerCoreLength, readPairs, cachedReadsContainingTarget, cachedReadsTargetIdxs, cachedReadPairs, cachedFullPairCoverage,
                                                   readStats, out extension.terminatingPrimerReached, out extension.coinTossed, out extension.abandoned, out extension.cost, out extension.meanDepth, out extension.avgDepth, out extension.levelAtTP, breakOnRead, logees);
                if (logees != null && breakOnRead)
                    logees.Add(new string(' ', level * 4) + level + ": " + "extended " + "[" + nextBase + "]" + " d=" + depth +
                        " tp=" + extension.terminatingPrimerReached + " ct=" + extension.coinTossed + " ab=" + extension.abandoned + " cost=" + extension.cost + " md=" + extension.meanDepth + " avg=" + extension.avgDepth + " latTP=" + extension.levelAtTP +
                        " frc= " + fewestRecursiveCalls + " " + extension.extendedSeq.ToString());
            }

            return extension;
        }

        // get the context depths for a given context length, and returns the max depth
        private static int GetContextDepths(Sequence seq, List<ulong> kMerVariants, int kMerSize, List<bool> variantAlive, int minDepth, List<int> contextLengths, int pl, Dictionary<int, Dictionary<ulong, int>> kMerContexts, List<int> contextCounts)
        {
            int maxContextDepth = 0;

            int contextLength = contextLengths[pl];
            if (contextLength > seq.Length)
                return 0; 

            int startOfContext = seq.Length - contextLength + 1;   

            for (int b = 0; b < kMerVariants.Count; b++)
            {
                contextCounts[b] = 0;

                if (variantAlive[b])
                {
                    ulong hashedContext = HashContextVariant(seq, startOfContext, contextLength, kMerSize, kMerVariants[b]);
                    int contextDepth;
                    kMerContexts[contextLengths[pl]].TryGetValue(hashedContext, out contextDepth);
                    contextCounts[b] = contextDepth;

                    if (contextDepth >= minDepth && contextDepth > maxContextDepth)
                        maxContextDepth = contextDepth;
                }
            }

            return maxContextDepth;
        }

        private static int FindChosenVariant(List<bool> variantAlive)
        {
            int chosenIdx = 0;
            for (int v = 0; v < variantAlive.Count; v++)
                if (variantAlive[v])
                    chosenIdx = v;
            return chosenIdx;
        }

        // fast check to cull non-viable alternative kMers
        private static int CheckVariantViability(Sequence seq, int m, int kMerSize, List<ulong> kMerVariants, List<int> kMerVariantDepths, List<int> kMerContextDepths, int minDepthForRead, List<bool> variantAlive,
                                                 Dictionary<int, Dictionary<ulong, int>> kMerContexts, List<int> contextLengths, int minViableContextIdx, out string failureReason)
        {
            int viableAlternatives = 0;
            failureReason = null;

            // and stop if none of the extensions are at all viable
            int maxVariantDepth = 0;
            variantAlive.Clear();

            for (int b = 0; b < kMerVariants.Count; b++)
            {
                int kMerVariantDepth = kMerVariantDepths[b];
                if (kMerVariantDepth > maxVariantDepth)
                    maxVariantDepth = kMerVariantDepth;
            }

            for (int b = 0; b < kMerVariants.Count; b++)
            {
                int kMerVariantDepth = kMerVariantDepths[b];
                if (kMerVariantDepth >= minDepthForRead) 
                {
                    viableAlternatives++;
                    variantAlive.Add(true);
                }
                else
                    variantAlive.Add(false);

                kMerContextDepths[b] = 0;
            }
 
            if (viableAlternatives == 0)
                failureReason = "no next";

            // can we discard some variants based on a longer kMer? (context)
            if (viableAlternatives > 1 && kMerContexts != null)
            {
                // check the viable alternatives for shortest context support as well
                int contextLength = contextLengths[minViableContextIdx];
                Dictionary<ulong, int> kMerContextsTable = kMerContexts[contextLength];
                int startOfContext = m + kMerSize - contextLength;
                int minDepthForContext = minDepthForRead / 2;
                if (minDepthForContext < 1)
                    minDepthForContext = 1;
                if (viableAlternatives == 1)
                    minDepthForContext = 1;

                if (startOfContext >= 0)
                {
                    viableAlternatives = 0;

                    for (int b = 0; b < kMerVariants.Count; b++)
                    {
                        if (variantAlive[b])
                        {
                            ulong hashedContext = HashContextVariant(seq, startOfContext, contextLength, kMerSize, kMerVariants[b]);
                            int contextCount;
                            kMerContextsTable.TryGetValue(hashedContext, out contextCount);
                            bool contextPresent = contextCount >= minDepthForContext;
                            kMerContextDepths[b] = contextCount;

                            if (contextPresent)
                            {
                                viableAlternatives++;
                            }
                            else
                            {
                                kMerVariantDepths[b] = 0;
                                variantAlive[b] = false;
                            }
                        }
                    }

                    if (viableAlternatives == 0)
                        failureReason = "minContext";
                }
            }

            return viableAlternatives;
        }

        private static bool CheckExtensionViability(string seq, Dictionary<ulong, int> kMerTable, int kMerSize, int wantedReadLength, int minDepthAllowed)
        {
            ulong kMer;
            bool kMerValid = kMers.CondenseMer(seq, seq.Length - kMerSize, kMerSize, out kMer);
            bool targetGroup = false;
            //if (targetGroup)
            //    Console.WriteLine(seq);
            bool viable = kMerValid && ExploreExtensionTree(kMer, kMerTable, kMerSize, wantedReadLength, minDepthAllowed, seq.Length, 0, 0, targetGroup);
            //if (!viable && targetGroup)
            //    Console.WriteLine(seq);
            return viable;
        }

        private static bool ExploreExtensionTree(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize, int wantedReadLength, int minDepthAllowed, int extendedLength, int recursiveDepth, int indent, bool targetGroup)
        {
            bool viable = true;
            ulong lbMask = 0xffffffffffffffff << (32 - kMerSize + 1) * 2;
            ulong kMerPrefix = (kMer << 2) & lbMask;
            List<ulong> kMersFound = new List<ulong>();
            const int maxRecursiveDepth = 5;

            while (viable)
            {
                int viableCount = 0;
                ulong lastGoodVariant = 0;
                kMersFound.Clear();
                viable = false;

                for (int b = 0; b < 4; b++)
                {
                    ulong kMerVariant = kMerPrefix | ((ulong)b << (32 - kMerSize) * 2);
                    int kMerDepth = kMerTableLookup(kMerVariant, kMerTable, kMerSize);
                    
                    if (kMerDepth >= minDepthAllowed)
                    {
                        viableCount++;
                        lastGoodVariant = kMerVariant;
                        kMersFound.Add(kMerVariant);
                    }
                }

                if (viableCount == 0)
                {
                    viable = false;
                    if (targetGroup) 
                        Console.WriteLine(new string(' ', indent+1) + kMers.ExpandMer(kMerPrefix, kMerSize-1) + "_ X");
                    break;
                }
                if (viableCount > 1)
                {
                    // check if any of the candidate kMers appears in the candidate long-enough reads
                    if (recursiveDepth == maxRecursiveDepth)
                    {
                        viable = false;
                        break;
                    }
                    foreach (ulong kMerVariant in kMersFound)
                    {
                        if (targetGroup) 
                            Console.WriteLine(new string(' ', indent+1) + kMers.ExpandMer(kMerVariant, kMerSize) + " + " + kMerTableLookup(kMerVariant, kMerTable, kMerSize));
                        viable |= ExploreExtensionTree(kMerVariant, kMerTable, kMerSize, wantedReadLength, minDepthAllowed, extendedLength, recursiveDepth + 1, indent+1, targetGroup);
                    }
                    break;
                }
                // just one viable so continue extending
                extendedLength++;
                indent++;
                viable = true;
                if (targetGroup) 
                    Console.WriteLine(new string(' ', indent) + kMers.ExpandMer(lastGoodVariant, kMerSize) + "  " + kMerTableLookup(lastGoodVariant, kMerTable, kMerSize));
                kMerPrefix = (lastGoodVariant << 2) & lbMask;
                if (extendedLength >= wantedReadLength)
                {
                    if (targetGroup) 
                        Console.WriteLine(new string(' ', indent) + kMers.ExpandMer(lastGoodVariant, kMerSize) + " =");
                    break;
                }
            }

            return viable;
        }

        private static bool kMerIsTerminating(ulong kMer, int kMerSize, HashSet<ulong>[] terminatingPrimers, int primerHeadlength, int primerCoreLength)
        {
            // extract primer from RHS end of kMer (leave left-adjusted)
            kMer = kMer << (kMerSize - primerHeadlength - primerCoreLength) * 2;
            // cccccccccchhhh
            ulong coreMask = 0xFFFFFFFFFFFFFFFF << (32 - primerCoreLength) * 2;
            ulong potentialPrimerHead = kMer << primerCoreLength * 2;
            ulong potentialPrimerCore = kMer & coreMask;

            return terminatingPrimers[head].Contains(potentialPrimerHead) && terminatingPrimers[core].Contains(potentialPrimerCore);
        }

        // Find a kMer in the canonical kMer table.
        private static int kMerTableLookup(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize)
        {
            kMer = Canonical(kMer, kMerSize);

            kMerTable.TryGetValue(kMer, out int count);

            return count;
        }

        private static int GetDeepestVariantDepth(VariantsWanted vw, ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize, List<ulong> kMerVariants, List<int> kMerVariantDepths, int noiseLevel)
        {
            int deepestVariantDepth = 0;

            if (vw == VariantsWanted.allSingle)
                GetAllMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths, noiseLevel);
            if (vw == VariantsWanted.lastSingle)
                GetNextMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths);
            if (vw == VariantsWanted.allDouble)
            {
                GetAllMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerVariantDepths, -1);
                GetAllDoubleMerVariantsAndCounts(kMerVariants, kMerVariantDepths, kMerTable, kMerSize, noiseLevel);
            }

            for (int b = 0; b < kMerVariants.Count; b++)
            {
                int variantDepth = kMerVariantDepths[b];
                if (variantDepth > deepestVariantDepth)
                    deepestVariantDepth = variantDepth;
            }

            return deepestVariantDepth;
        }
    }

    public class ReadStats
    {
        public int avgDepth;            // average depth over the seq
        public int meanDepth;           // harmonic mean (reduces the effect of conserved kMers)
        public int minDepthAllowed;     // rock-bottom depth allowed. Depths lower than this are regarded as zero
        public int minDepthFound;       // actual min depth found in the seq
        public int initialGoodDepth;    // guess at an initial good depth
    }

    public class Extension
    {
        public string startingRead;
        public Sequence extendedSeq;
        public bool terminatingPrimerReached;
        public int levelAtTP;
        public bool abandoned;
        public int cost;                // not used in 2.3.0, may resurrect with better cost calculation at a later time
        public bool coinTossed;
        public bool alreadyCached;
        public int cachedAtRead;
        public int meanDepth;
        public int avgDepth;
    }
}
