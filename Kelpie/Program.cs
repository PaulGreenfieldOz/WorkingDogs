using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using WorkingDogsCore;

namespace Kelpie
{
    // Kelpie runs in two phases: 
    //      - extracting and trimming the reads that cover the inter-primer region
    //      - extending the starting reads just found, using kMers derived from the covering reads
    //
    // usage: Kelpie -f forwardPrimer -r reversePrimer readsToFilterFNP extendedReadsFN
    // 
    // Extracting/trimming 
    // -------------------
    //
    // Takes reads filtered for 16S or 18S (via FilterReads or an equivalent HMM) and finds those reads that include one of the
    // specified primers at their start or end. These reads are then tiled for kMers to build an initial filter set. The remaining reads
    // are passed over this filter to find reads that overlap with the starting/ending reads. These additional reads are tiled and their
    // kMers are added to the filter. This process is repeated until the rate of new kMer discovery has dropped sufficiently. If pairs of files
    // were provided, the kMer sets are merged and only kMers found in both sets are kept.
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
    // and the terminal primers are (F) and (R). Reads containing one of these starting primers are found, trimmed to the start of the
    // respective primer, and tiled for the initial set of filter kMers. These selected reads are added to the list of reads to be written.
    // The remaining reads are then scanned for any that start with a kMer in the starting set, and these reads are removed from the to-be-scanned list,
    // tiled and their kMers are added to the filtering set. Reads containing a terminating primer are only tiled up to this primer. This process
    // continues until there are either no new matching reads (and so no changes to the filter), or the rate of new reads looks like it is indicating
    // that we've reached the end and may be about to shoot past one or more unrecognised terminal primers.
    // 
    // Once the filter set is built, all the reads are passed by it again and matching reads are added to the to-be-written list. These
    // reads are then written. Those reads that start with a forward primer are marked with FP, and those that end with (F) are marked with 
    // FP'. The FP' reads will be reverse-complemented before being extended, ensuring all such 'starting' reads start with the forward primer 
    // (and can all be extended left-to-right). These starting reads are marked so they can be extended in the next phase.
    //
    // Extending
    // ---------
    //
    // The inter-primer reads found in the previous stage are then turned into a set of 32-mers, and this is then denoised to remove very low coverage 
    // kMers. The other longer kMers are then generated from the reads, their 32-ends are checked against the denoised 32-mer table, and then stored in
    // the appropriate kMer table. kMers are generated for each length from 32 to the read length, in steps of 8 bases. The long kMers are saved as hashes 
    // of the 32-mers found at the start and end. This allows the use of fast binary hash tables.
    //
    // Each of the starting reads is checked, trimmed and then extended, one base at a time. At each iteration of the extension loop, the ending 32-mer for 
    // each possible extension is generated and then checked forviability against the 32-mer hash table. If more than one extension base is viable, all of
    // the longer kMers are generated and checked. If only a single extension is viable at any point, this base is chosen and the extension loop moves on.
    // The extension loop terminates when either the extended read ends with a terminal primer, or when none of the possible extensions are viable.
    //
    // If there are multiple viable extensions, each of the viable extended reads is recursively explored to see of any of them end up finishing with a terminal
    // primer. If only one of the extensions can reach the end of the inter-primer region, it is chosen and the extension of the reads finishes. If more than one
    // extension can reach a terminal primer, one of them is chosen at random, in proportion to the repetition depth of the 32-mers at the end of the extended
    // read. If none of the extension reach a terminal primer, the extension that results in the longest extended read is cosen, although this will often
    // result in an incomplete extension that will be discarded anyway.
    //
    // Finally the set of extended reads are trimmed of their primers and written, and any reads that didn't reach a terminal primer are dropped. 
    // 

    class Program
    {
        const string version = "V1.0.3";

        const bool cacheThings = true;

        enum PairedReadsWanted
        {
            unspecified,
            paired,
            unpaired
        }

        const int kMerSize = 32;                // size of kMers used in filters
        const int startSize = 5;                // how close to the start of a read do we expect the first found kMer
        const int shortestPairSize = 40;        // shortest long kMer/kMer pair length
        const int maxRecursion = 20;            // limit depth of recursion during downstream tree exploration

        const int statsNextIdx = 0;             // how many times did we look for the next kMer
        const int statsSingleIdx = 1;           // only a single next kMer
        const int statsDownstreamIdx = 2;       // needed to look downstream to resolve ambiguity
        const int statsSingleDSTIdx = 3;        // only a single downstream path led to a terminal primer
        const int statsRolledIdx = 4;           // chose kMer in proportion to kMer counts
        const int statsLongestIdx = 5;          // chose kMer that gave the longest downstream extension
        const int statsCachedResult = 6;        // starting read already extended
        const int statsCachedBranch = 7;        // cached branch encountered in tree exploration
        const int statsPairsIdx = 8;            // first pair-length index

        static char[] bases = new char[] { 'A', 'C', 'G', 'T' };

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.WriteLine("usage: Kelpie [-h] -f forwardPrimer -r reversePrimer filteredReadsFNP extendedReadsFN (" + version + ")");
                return;
            }

            string forwardPrimer = null;
            string reversePrimer = null;
            string primerTag = null;
            int maxReadLength = 0;
            bool writeFullHelp = false;
            bool saveLog = false;
            int minLength = 0;
            int maxLength = Int32.MaxValue;
            int minDepth = 2;
            string extendedReadsFN = null;
            bool generateStats = false;
            bool strict = false;

            PairedReadsWanted pairedReadsArg = PairedReadsWanted.unspecified;
            bool pairedReads = false;

            List<string> FNParams = new List<string>();     // the set of file names or patterns to filter

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

                    if (lcArg == "-stats")
                    {
                        generateStats = true;
                        continue;
                    }

                    if (lcArg == "-strict")
                    {
                        strict = true;
                        continue;
                    }

                    if (lcArg == "-f" || lcArg == "-forward")
                    {
                        if (!CheckForParamValue(p, args.Length, "forward primer expected after -forward"))
                            return;
                        forwardPrimer = args[p + 1];
                        p++;
                        continue;
                    }

                    if (lcArg == "-r" || lcArg == "-reverse")
                    {
                        if (!CheckForParamValue(p, args.Length, "reverse primer expected after -reverse"))
                            return;
                        reversePrimer = args[p + 1];
                        p++;
                        continue;
                    }

                    if (lcArg == "-save")
                    {
                        if (!CheckForParamValue(p, args.Length, "primer name expected after -save"))
                            return;
                        primerTag = args[p + 1];
                        p++;
                        continue;
                    }

                    if (lcArg == "-paired")
                    {
                        pairedReadsArg = PairedReadsWanted.paired;
                        continue;
                    }

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
                            maxReadLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -maxreadlength parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (lcArg == "-log")
                    {
                        saveLog = true;
                        continue;
                    }

                    if (lcArg == "-min" || lcArg == "-minlength")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a minimum length after -min");
                            return;
                        }
                        try
                        {
                            minLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a minimum length after -min: " + args[p + 1]);
                            return;
                        }
                        if (minLength < shortestPairSize)
                        {
                            Console.WriteLine("min length must be >= " + shortestPairSize);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (lcArg == "-max" || lcArg == "-maxlength")
                    {
                        if (p + 1 == args.Length)
                        {
                            Console.WriteLine("expected a maximum length after -max");
                            return;
                        }
                        try
                        {
                            maxLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a maximum length after -max: " + args[p + 1]);
                            return;
                        }
                        if (maxLength < shortestPairSize)
                        {
                            Console.WriteLine("max length must be >= " + shortestPairSize);
                            return;
                        }
                        p++;
                        continue;
                    }

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

            if (forwardPrimer == null)
            {
                Console.WriteLine("no forward primer found");
                return;
            }

            if (reversePrimer == null)
            {
                Console.WriteLine("no reverse primer found");
                return;
            }

            if (FNParams.Count == 0)
            {
                Console.WriteLine("no file names found");
                return;
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("expected both filteredReadsFNP & extendedReadsFN. Got only 1 file name");
                return;
            }

            List<string> FNs = new List<string>();
            extendedReadsFN = FNParams[FNParams.Count - 1];

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
                pairedReads = true;
            // even number of files - assume R1/R2 pairs - R2s will be reversed for filtering
            if (FNs.Count % 2 == 0 && pairedReadsArg == PairedReadsWanted.unspecified)
                pairedReads = true;

            StreamWriter log = null;
            if (saveLog)
                log = new StreamWriter("EPRlog.txt");

            Stopwatch totalElapsedTime = new Stopwatch();
            totalElapsedTime.Start();
            Stopwatch phase1ElapsedTime = new Stopwatch();
            phase1ElapsedTime.Start();

            char[] degenerateBases;
            Dictionary<char, List<char>> degenerateBaseExpansions;
            kMers.InitialiseDegenerateBaseTables(out degenerateBases, out degenerateBaseExpansions);

            // starting primers are either the Forward primer or the Reverse primer 
            HashSet<string> startingPrimers = new HashSet<string>();
            // filter extension stops when we hit a termination primer - rc(Reverse) or rc(Forward)
            HashSet<string> terminatingPrimers = new HashSet<string>();
            // and remember just the forward primers (and their RC forms) so we can mark reads containing them
            HashSet<string> forwardPrimers = new HashSet<string>();

            // the kMers in these tables are all primerLength strings (length of shortest primer)
            // the scanning code needs a common length as the starting/termination primers contain both forward and reverse primers, 
            int primerLength = Math.Min(forwardPrimer.Length, reversePrimer.Length);
            // primer core length - balance between being long enough for distinctiveness and short enough to allow some bases for the 'head' 
            const int primerCoreLength = 15;
            int primerHeadLength = primerLength - primerCoreLength;

            // primers are split into two parts - head & core - to better emulate PCR (and primer binding) and so to allow more mismatches
            // The 3' end (RHS) binds more tighly and can tolerate the fewest differences. The 5' end binds less tighly and can tolerate more mismatches.
            // We allow 2 mismatches in the 3' (core) end of each primer, and half the bases to mismatch at the 5' (head) end.
            // Allowing 3 core mismatches gave rise to 23S matches!
            // 

            // forward primer cores
            // --------------------
            // trim forward primer to the uniform primer length (keeping the RHS)
            string trimmedForwardPrimer = forwardPrimer;
            if (forwardPrimer.Length > primerLength)
                trimmedForwardPrimer = forwardPrimer.Substring(forwardPrimer.Length - primerLength);
            string forwardPrimerCore = trimmedForwardPrimer.Substring(primerHeadLength);
            string forwardPrimerHead = trimmedForwardPrimer.Substring(0, primerHeadLength);
            // build variant tables for head and core
            HashSet<string> forwardPrimerHeadVariants = new HashSet<string>();
            HashSet<string> forwardPrimerCoreVariants = new HashSet<string>();
            int fphAllowedMismatches = forwardPrimerHead.Length / 2;
            if (fphAllowedMismatches > 2)
                fphAllowedMismatches = 2;
            kMers.GenerateSeqVariants(forwardPrimerHead, forwardPrimerHead.Length, degenerateBases, degenerateBaseExpansions, forwardPrimerHeadVariants, fphAllowedMismatches);
            kMers.GenerateSeqVariants(forwardPrimerCore, forwardPrimerCore.Length, degenerateBases, degenerateBaseExpansions, forwardPrimerCoreVariants, 2);
            // and merge them to form the starting and terminating primers derived from the forward primer
            foreach (string forwardPrimerHeadVariant in forwardPrimerHeadVariants)
                foreach (string forwardPrimerCoreVariant in forwardPrimerCoreVariants)
                {
                    string primerVariant = forwardPrimerHeadVariant + forwardPrimerCoreVariant;
                    if (!startingPrimers.Contains(primerVariant))
                    {
                        startingPrimers.Add(primerVariant);
                        forwardPrimers.Add(primerVariant);
                    }
                    string primerVariantRC = kMers.ReverseComplement(primerVariant);
                    if (!terminatingPrimers.Contains(primerVariantRC))
                        terminatingPrimers.Add(primerVariantRC);
                    if (!forwardPrimers.Contains(primerVariantRC))
                        forwardPrimers.Add(primerVariantRC);
                }
            Console.WriteLine("built kMer sets from forward primer");

            // reverse primer cores
            // --------------------
            // trim reverse primer to the uniform primer length
            string trimmedReversePrimer = reversePrimer;
            if (reversePrimer.Length > primerLength)
                trimmedReversePrimer = reversePrimer.Substring(reversePrimer.Length - primerLength);
            string reversePrimerCore = trimmedReversePrimer.Substring(primerHeadLength);
            string reversePrimerHead = trimmedReversePrimer.Substring(0, primerHeadLength);
            // build variant tables for head and core
            HashSet<string> reversePrimerHeadVariants = new HashSet<string>();
            HashSet<string> reversePrimerCoreVariants = new HashSet<string>();
            int rphAllowedMismatches = reversePrimerHead.Length / 2;
            if (rphAllowedMismatches > 2)
                rphAllowedMismatches = 2;
            kMers.GenerateSeqVariants(reversePrimerHead, reversePrimerHead.Length, degenerateBases, degenerateBaseExpansions, reversePrimerHeadVariants, rphAllowedMismatches);
            kMers.GenerateSeqVariants(reversePrimerCore, reversePrimerCore.Length, degenerateBases, degenerateBaseExpansions, reversePrimerCoreVariants, 2);
            // and merge them to form the starting primers derived from the reverse primer
            foreach (string reversePrimerHeadVariant in reversePrimerHeadVariants)
                foreach (string reversePrimerCoreVariant in reversePrimerCoreVariants)
                {
                    string primerVariant = reversePrimerHeadVariant + reversePrimerCoreVariant;
                    if (!startingPrimers.Contains(primerVariant))
                        startingPrimers.Add(primerVariant);
                    string primerVariantRC = kMers.ReverseComplement(primerVariant);
                    if (!terminatingPrimers.Contains(primerVariantRC))
                        terminatingPrimers.Add(primerVariantRC);
                }
            Console.WriteLine("built kMer sets from reverse primer");

            int totalReadsCount = 0;

            // process the pairs of files 
            // --------------------------
            //      assume that we're dealing with a single DNA sample, spread (possibly) across multiple (pairs of) files.
            //      the usual case will be a single R1/R2 pair, but some HiSeq runs produce two pairs of pairs (L1 & L2) for each sample.
            //      the R1 reads are effectively merged, as are the R2 reads.
            //
            //      the starting and ending filters are built separately for the R1 and R2 files (if 'paired' is set)
            //      if 'paired', reduce the filters by excluding any kMers that were seen in only the R1 or R2 set.

            // if we're handling paired reads file, loop over the files two at a time
            int fileStride = pairedReads ? 2 : 1;

            // A initial set of filtering kMers (starting from the base after a detected starting primer and continuing at single-base increments)
            // The filter is built by scanning the starts of both the R1 and R2 reads (if paired) and only those
            // kMers present in both read sets will be kept for the final filtering pass.
            // The kMers in these tables are kMerLength strings (longer than minimum primer length)
            HashSet<string>[] startingFilter = new HashSet<string>[fileStride];
            // And the set of kMers that follow (and include the terminal primers). Used to stop building the starting filter when we get to the end of the region
            HashSet<string>[] endingFilter = new HashSet<string>[fileStride];

            // save all of the reads in a (paired) set of files (read and kept to make it faster to do multiple passes over them)
            List<string>[] headersToFilter = new List<string>[fileStride];
            List<string>[] readsToFilter = new List<string>[fileStride];
            // track which reads have already been tiled for the filter (don't want to change the actual reads set)
            HashSet<int>[] readsAlreadyProcessed = new HashSet<int>[fileStride];

            for (int fip = 0; fip < fileStride; fip++)
            {
                headersToFilter[fip] = new List<string>();
                readsToFilter[fip] = new List<string>();
                startingFilter[fip] = new HashSet<string>();
                endingFilter[fip] = new HashSet<string>();
                readsAlreadyProcessed[fip] = new HashSet<int>();
            }

            // read all the reads from all the files (possibly in pairs)
            for (int f = 0; f < FNs.Count; f += fileStride)
            {
                // read the headers and reads for this (pair) of files
                for (int fip = 0; fip < fileStride; fip++)
                {
                    string FN = FNs[f + fip];
                    int fileFormat = SeqFiles.DetermineFileFormat(FN);
                    StreamReader reads = new StreamReader(FN);
                    Console.WriteLine("filtering " + FN);

                    bool EOF = false;
                    while (!EOF)
                    {
                        string header;
                        string read = SeqFiles.ReadRead(reads, fileFormat, out header);
                        if (read == null)
                            break;

                        totalReadsCount++;
                        headersToFilter[fip].Add(header);
                        if (maxReadLength > 0 && read.Length > maxReadLength)
                            read = read.Substring(0, maxReadLength);
                        readsToFilter[fip].Add(read);
                    }
                }
            }

            // remember which reads were found to start with one of the primers. These reads are saved, removed and later written separately as they are not filtered again. 
            List<string>[] startingHeaders = new List<string>[fileStride];
            List<string>[] startingReads = new List<string>[fileStride];
            int startingReadsCount = 0;
            for (int fip = 0; fip < fileStride; fip++)
            {
                startingHeaders[fip] = new List<string>();
                startingReads[fip] = new List<string>();

                // Get an initial set of reads with primers near their starts. 
                // These will be trimmed and used to generate the initial starting filter set.
                // The terminal primers are also found and used to populate the ending kMer set (along with the terminal primers themselves)
                FindReadsWithPrimers(headersToFilter[fip], readsToFilter[fip], primerLength, startingPrimers, terminatingPrimers, forwardPrimers, endingFilter[fip], startingHeaders[fip], startingReads[fip], readsAlreadyProcessed[fip], kMerSize);

                startingReadsCount += startingReads[fip].Count;
            }
            Console.WriteLine("found " + startingReadsCount + " reads matching starting/terminating primers");

            // build region-between-primers filter, starting from these starting-primer reads and extending until we reach terminating primers
            for (int fip = 0; fip < fileStride; fip++)
            {
                // add the starts of the starting-primer-containing reads to the filter
                ExtendFilterSet(startingReads[fip], startingFilter[fip], kMerSize, endingFilter[fip]);

                int readsAddedToFilterLastTime = 1000000;
                bool stopExtending = false;
                int filterIterations = 0;
                while (!stopExtending)
                {
                    filterIterations++;

                    List<string> matchingReads = FindMatchingReads(readsToFilter[fip], readsAlreadyProcessed[fip], kMerSize, startingFilter[fip]);

                    int readsAddingToFilter = ExtendFilterSet(matchingReads, startingFilter[fip], kMerSize, endingFilter[fip]);

                    // Stop when we either found no more reads to tile - or if we appear to be finding an increasing number of reads to tile
                    // The number of reads to tile should be approximately flat while traversing the inter-primer region and then drop rapidly
                    // as we get to the end of the region. An increase in this number indicates that we've overshot the terminating primer.
                    stopExtending = readsAddingToFilter == 0 || readsAddingToFilter > 2 * readsAddedToFilterLastTime;
                    readsAddedToFilterLastTime = readsAddingToFilter;
                }
                //Console.WriteLine("filter contains " + startingFilter[fip].Count + " kMers");
            }
            startingHeaders = null;
            startingReads = null;
            readsAlreadyProcessed = null;

            // Have now built a rough filter but it may include kMers derived from sequencing errors or overshoots.
            // We will now improve this filter by culling any kMers that were not found in forward and reverse forms.
            // And if we're running with -strict, we'll also demand that all kMers are found in both files in a file pair

            HashSet<string> finalStartingFilter = new HashSet<string>(startingFilter[0]);
            HashSet<string> finalEndingFilter = new HashSet<string>(endingFilter[0]);
            if (fileStride > 1)
            {
                finalStartingFilter.UnionWith(startingFilter[1]);
                finalEndingFilter.UnionWith(endingFilter[1]);
            }

            Console.WriteLine("initial starting filter contains " + finalStartingFilter.Count + " kMers");
            Console.WriteLine("initial ending filter contains " + finalEndingFilter.Count + " kMers");

            List<string> kMersToDelete = new List<string>();
            foreach (string kMer in finalStartingFilter)
            {
                string kMerRC = kMers.ReverseComplement(kMer);
                if (!finalStartingFilter.Contains(kMerRC))
                    kMersToDelete.Add(kMer);
            }
            finalStartingFilter.ExceptWith(kMersToDelete);

            // strict pairs - kMers must also exist (in either form) in both files in the pair
            if (strict && fileStride > 1)
            {
                kMersToDelete.Clear();
                foreach (string kMer in finalStartingFilter)
                {
                    bool presentInBoth = (startingFilter[0].Contains(kMer) || startingFilter[0].Contains(kMers.ReverseComplement(kMer))) &&
                                         (startingFilter[1].Contains(kMer) || startingFilter[1].Contains(kMers.ReverseComplement(kMer)));
                    if (!presentInBoth)
                        kMersToDelete.Add(kMer);
                }
                finalStartingFilter.ExceptWith(kMersToDelete);
            }

            startingFilter = null;
            endingFilter = null;
            kMersToDelete = null;

            Console.WriteLine("final starting filter contains " + finalStartingFilter.Count + " kMers");
            Console.WriteLine("final ending filter contains " + finalEndingFilter.Count + " kMers");

            // now do a final pass though the starting reads, matching the starts of the reads against the starting filter
            // and writing out the selected reads if requested. The selected reads are saved for the extension phase..
            List<string> selectedHeaders = new List<string>();
            List<string> selectedReads = new List<string>();
            int readsAlreadySelected = 0;

            for (int fip = 0; fip < fileStride; fip++)
            {
                // go through all of the reads and filter out those that match the inter-primer kMers
                FinalFilterReads(headersToFilter[fip], readsToFilter[fip], kMerSize, finalStartingFilter, selectedHeaders, selectedReads);
                
                // and write them out if requested
                if (primerTag != null)
                    WriteFilteredReads(FNs[fip], primerTag, readsAlreadySelected, selectedHeaders, selectedReads, kMerSize);

                readsAlreadySelected = selectedReads.Count;
            }

            headersToFilter = null;
            readsToFilter = null;
            finalStartingFilter = null;
            finalEndingFilter = null;

            phase1ElapsedTime.Stop();
            Console.WriteLine("filtered and kept " + selectedReads.Count + "/" + totalReadsCount + " primer-region reads from " + FNs.Count + " files in " + phase1ElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s");

            // ====================
            // Read extension phase
            // ====================

            Stopwatch phase2ElapsedTime = new Stopwatch();
            phase2ElapsedTime.Start();

            // build the (packed) terminating primers (derived from the reverse primer only for the extension phase)
            HashSet<ulong> extensionTerminatingPrimers = new HashSet<ulong>();
            foreach (string reversePrimerHeadVariant in reversePrimerHeadVariants)
                foreach (string reversePrimerCoreVariant in reversePrimerCoreVariants)
                {
                    string primerVariant = reversePrimerHeadVariant + reversePrimerCoreVariant;
                    string primerVariantRC = kMers.ReverseComplement(primerVariant);
                    ulong primerVariantRCPacked = kMers.CondenseMer(primerVariantRC);
                    if (!extensionTerminatingPrimers.Contains(primerVariantRCPacked))
                        extensionTerminatingPrimers.Add(primerVariantRCPacked);
                }

            Console.WriteLine("built terminating primer set for extension phase");

            // process the filtered reads and remember all of the reads containing a forward primer (;FPx) as 'starting' reads 
            // RC the ;FP' starting reads so that they also have their forward primer at their start
            // tile each read on the way through and add its kMers to the kMer table

            // just the reads with forward primers
            List<string> readsWithStartingPrimers = new List<string>();
            // and all the other reads
            List<string> nonStartingReads = new List<string>();
            // remember the forward primers that were found so they can be trimmed before writing
            Dictionary<string, int> forwardPrimersFound = new Dictionary<string, int>();

            // a canonical kMers table
            Dictionary<ulong, int> kMerTable = new Dictionary<ulong, int>(100000);
            // save kMer pairs for a range of lengths. Shorter ones have a better chance of being found, but are less likely to be distinct
            List<int> pairSizes = new List<int>();
            Dictionary<int, Dictionary<ulong, int>> kMerPairs = new Dictionary<int, Dictionary<ulong, int>>();
            // remember the start of the inter-primer regions so we can scavenge reads that don't have complete primers at their start
            HashSet<string> startsOfRegion = new HashSet<string>();

            // tile the filtered reads (both R1 and R2)
            Console.WriteLine("building kMer table for extension phase");
            ulong[] kMersFromRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int[] depthsForRead = new int[1000];

            for (int r = 0; r < selectedReads.Count; r++)
            {
                string header = selectedHeaders[r];
                string read = selectedReads[r];

                //if (read == "GGATTTCACCCCCGACTTAACAATCCGCCTACGCGCGCTTTACGCCCAGTAATTCCGAACAACGCTAGCCCCCTCCGTATTACCGCGGCTGCTGGCAC")
                //    Debugger.Break();

                if (header.EndsWith(";FP"))
                {
                    readsWithStartingPrimers.Add(read);
                }
                else if (header.EndsWith(";FP'"))
                {
                    read = kMers.ReverseComplement(read);
                    readsWithStartingPrimers.Add(read);
                }
                else
                    nonStartingReads.Add(read);

                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);

                //if (read == "GTGCCAGCAGCCGCGGTAATACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTGTTAAGTCGGGGGTGAAATCC")
                //    Debugger.Break();

                // initial set of as-read and canonical kMers
                for (int k = 0; k < kMerCount; k++)
                    if (kMersValid[k])
                    {
                        //if (kMersFromRead[k] == 0x667F1952C3D60419 || kMersFromRead[k] == 0x9BEF683C7A9B0266 || kMersFromRead[k] == 0x5C6667F1952C3D60 || kMersFromRead[k] == 0xF683C7A9B02666CA)
                        //{
                        //    Console.WriteLine(kMersFromRead[k].ToString("X16") + " " + kMers.ExpandMer(kMersFromRead[k], 32) + " from " + header);
                        //    Debugger.Break();
                        //}

                        AddMerToTable(Canonical(kMersFromRead[k], kMerSize), kMerTable, kMerSize);
                    }
            }
            
            // remove kMers that appear to be sequencing errors
            int noiseLevel;
            int meanDepth;
            int kMersCulled = DeNoiseMerTable(selectedReads, kMerSize, kMerTable, minDepth, out noiseLevel, out meanDepth, log);
            Console.WriteLine("removed " + kMersCulled + " kMers from kMerTable");

            int totalPairs = 0;
            int pairsDropped = 0;

            // re-process the reads again and build the pairs tables
            for (int r = 0; r < selectedReads.Count; r++)
            {
                string header = selectedHeaders[r];
                string read = selectedReads[r];

                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);
                if (kMerCount > depthsForRead.Length)
                    Array.Resize<int>(ref depthsForRead, kMerCount + 1000);
                for (int k = 0; k < kMerCount; k++)
                {
                    ulong kMer = kMersFromRead[k];
                    depthsForRead[k] = kMerTableLookup(kMer, kMerTable, kMerSize);
                }

                for (int pl = shortestPairSize; pl < read.Length; pl += 8)
                {
                    // initialise when each pair length is first encountered
                    if (!kMerPairs.ContainsKey(pl))
                    {
                        pairSizes.Add(pl);
                        kMerPairs.Add(pl, new Dictionary<ulong, int>(10000));
                    }

                    int pairCount = read.Length - pl + 1;
                    for (int p = 0; p < pairCount; p++)
                    {
                        int secondMerIdx = p + pl - kMerSize;
                        if (depthsForRead[p] > 0 && depthsForRead[secondMerIdx] > 0)
                        {
                            AddMerToTable(HashPair(kMersFromRead[p], kMersFromRead[secondMerIdx], kMerSize), kMerPairs[pl], kMerSize);
                            //trace.Add("added " + kMers.ExpandMer(kMersFromRead[p], kMerSize) + "+" + kMers.ExpandMer(kMersFromRead[secondMerIdx], kMerSize) + " " +
                            //              depthsForRead[p] + " " + depthsForRead[secondMerIdx] + " " + HashPair(kMersFromRead[p], kMersFromRead[secondMerIdx], kMerSize).ToString("X16"));
                        }
                        else
                        {
                            pairsDropped++;
                            //trace.Add("dropped " + kMers.ExpandMer(kMersFromRead[p], kMerSize) + "+" + kMers.ExpandMer(kMersFromRead[secondMerIdx], kMerSize) + " " + 
                            //              depthsForRead[p] + " " + depthsForRead[secondMerIdx]);
                        }
                    }
                }
            }

            foreach (int ps in pairSizes)
                totalPairs += kMerPairs[ps].Count;

            // stats on pairSize usage
            int[] pairSizeStats = new int[pairSizes.Count + statsPairsIdx];

            if (generateStats)
            {
                foreach (int ps in pairSizes)
                    Console.WriteLine(ps + ": " + kMerPairs[ps].Count);
            }

            Console.WriteLine("finished building kMer tables: " + (kMerTable.Count - kMersCulled) + " " + kMerSize + "-mers; " +
                               totalPairs + " kMer pairs over " + pairSizes.Count + " pair lengths");

            // build final set of starting reads - after cleaning, trimming, extending, dropping
            List<Sequence> finalStartingReads = new List<Sequence>();

            int uncleanStartingReads = 0;
            int cleanStartingReads = 0;
            int shortStartingReads = 0;
            int extendedShortStartingReads = 0;

            for (int i = 0; i < readsWithStartingPrimers.Count; i++)
            {
                string startingRead = readsWithStartingPrimers[i];
                Sequence startingSeq = new Sequence(startingRead);

                // drop any reads containing Ns
                int firstNIdx = startingRead.IndexOf('N');
                if (firstNIdx >= 0)
                    continue;

                bool cleanRead = CleanStartingRead(kMerSize, kMerTable, pairSizes, kMerPairs, noiseLevel, meanDepth, startingSeq, log);

                // starting primer/kMer is OK
                if (cleanRead)
                {
                    cleanStartingReads++;

                    string fp = startingSeq.SubSeq(0, forwardPrimer.Length).ToString();
                    if (forwardPrimersFound.ContainsKey(fp))
                        forwardPrimersFound[fp]++;
                    else
                        forwardPrimersFound.Add(fp, 1);

                    string shortRead = null;
                    if (log != null)
                        shortRead = startingSeq.ToString();

                    if (startingSeq.Length < (shortestPairSize + primerLength))
                    {
                        shortStartingReads++;

                        // try extending the read
                        bool extended = ExtendShortStartingRead(startingSeq, kMerTable, kMerSize, (shortestPairSize + primerLength), pairSizes, kMerPairs, noiseLevel, meanDepth, log);

                        if (extended)
                        {
                            if (log != null)
                                log.WriteLine("extended(SSR): " + shortRead + " --> " + startingSeq.ToString());
                            extendedShortStartingReads++;
                        }
                    }
                    
                    // read long enough (primer + short kMer pair)
                    if (startingSeq.Length >= (shortestPairSize + primerLength))
                        finalStartingReads.Add(startingSeq);
                }
                else
                // start of potential starting read was poor so add it the non-starting set and it may be rescued later
                {
                    nonStartingReads.Add(startingRead);
                    uncleanStartingReads++;
                }
            }
            Console.WriteLine("cleaned " + cleanStartingReads + "/" + readsWithStartingPrimers.Count + " starting reads");
            Console.WriteLine("extended " + extendedShortStartingReads + "/" + shortStartingReads + " short starting reads");
            Console.WriteLine(finalStartingReads.Count + " reads starting with F primer after extension");
            readsWithStartingPrimers = null;           // finished with the initial set of starting reads

            // remember what follows the primers in the starting reads
            for (int i = 0; i < finalStartingReads.Count; i++)
            {
                string startingRead = finalStartingReads[i].ToString();

                // trim the forward primer and save the rest of the read (at least short kMer pair in length)
                startingRead = startingRead.Substring(primerLength, startingRead.Length - primerLength);

                string startOfRegion = startingRead.Substring(0, shortestPairSize);
                if (!startsOfRegion.Contains(startOfRegion))
                {
                    startsOfRegion.Add(startOfRegion);
                    //if (startOfRegion == "AGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGC")
                    //    Debugger.Break();
                }
            }

            // find any other reads that start with a partial forward primer followed by one of the starting seqs
            List<string> additionalStartingReads = new List<string>();
            foreach (string read in nonStartingReads)
            {
                //if (read == "GTGCCAGCAGCCGCGGTAATCGGGACGGTTCAAGCGTTATCCGGATTTACTAGGTTTAAAGGGTG")
                //    Debugger.Break();

                int firstNIdx = read.IndexOf('N');
                if (firstNIdx >= 0)
                    continue;

                for (int i = 0; i < forwardPrimer.Length; i++)
                {
                    if (i + shortestPairSize > read.Length)
                        break;

                    string startFragment = read.Substring(i, shortestPairSize);
                    if (startsOfRegion.Contains(startFragment))
                    {
                        additionalStartingReads.Add(read.Substring(i));
                        break;
                    }
                }
            }
            // and add these to the starting reads (after cleaning/extending etc)
            int addedStartingReads = 0;
            if (log != null)
                log.WriteLine("rescuing starting reads with poor first kMers");

            foreach (string additionalStartingRead in additionalStartingReads)
            {
                Sequence additionalStartingReadSeq = new Sequence(additionalStartingRead);
                bool cleanRead = CleanStartingRead(kMerSize, kMerTable, pairSizes, kMerPairs, noiseLevel, meanDepth, additionalStartingReadSeq, log);

                //if (additionalStartingRead.ToString() == "AGAGTTTGATCCTGGCTCAGGATGAACGCTAGCTACAGGCTTAACACATGCAAGTGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT")
                //    Debugger.Break();

                if (cleanRead)
                {
                    string shortRead = null;
                    if (log != null)
                        shortRead = additionalStartingRead;

                    if (additionalStartingRead.Length < shortestPairSize)
                    {
                        shortStartingReads++;

                        // try extending the read
                        bool extended = ExtendShortStartingRead(additionalStartingReadSeq, kMerTable, kMerSize, shortestPairSize, pairSizes, kMerPairs, noiseLevel, meanDepth, log);

                        if (extended)
                        {
                            if (log != null)
                                log.WriteLine("extended(SRR): " + shortRead + " --> " + additionalStartingReadSeq.ToString());
                            extendedShortStartingReads++;
                        }
                    }
                    // read long enough so add it to the starting set
                    if (additionalStartingReadSeq.Length >= shortestPairSize)
                    {
                        finalStartingReads.Add(additionalStartingReadSeq);
                        addedStartingReads++;
                    }
                }
            }
            Console.WriteLine("added " + addedStartingReads + " starting reads");

            // Starting from a set of starting reads and incrementally growing them from the available kMers until they can be grown no further
            Dictionary<Sequence, int> extendedReadsSet = ExtendFromStartingReads(finalStartingReads, kMerSize, kMerTable, pairSizes, kMerPairs, pairSizeStats,
                                                                                 extensionTerminatingPrimers, primerLength, noiseLevel, meanDepth, log);

            // now have a set of final extended reads - almost ready for writing
            // but first trim away the starting and terminating primers (and discard any reads that are too short or too long)
            Dictionary<string, int> trimmedExtendedReads = TrimExtendedReads(extendedReadsSet, forwardPrimer.Length, primerLength, forwardPrimersFound, extensionTerminatingPrimers, minLength, maxLength, "F-->R", log);

            // and write out all the extended reads
            StreamWriter extendedReads = new StreamWriter(extendedReadsFN);

            int readNo = 1;
            foreach (KeyValuePair<string, int> kvp in trimmedExtendedReads)
            {
                string extendedRead = kvp.Key;
                int count = kvp.Value;
                for (int r = 0; r < count; r++)
                {
                    SeqFiles.WriteRead(extendedReads, (">R" + readNo), extendedRead, SeqFiles.formatFNA);
                    readNo++;
                }
            }

            extendedReads.Close();
            if (log != null)
                log.Close();

            phase2ElapsedTime.Stop();
            totalElapsedTime.Stop();

            Console.WriteLine("extended and wrote " + (readNo - 1) + " extended F->R reads (" + (minLength > 0 ? (">=" + minLength + "-mers") : "full length") +
                              ") in " + phase2ElapsedTime.Elapsed.TotalSeconds.ToString("F1") + "s");

            if (generateStats)
            {
                Console.WriteLine("read extensions:\t" + pairSizeStats[statsNextIdx]);
                Console.WriteLine("read cached:\t" + pairSizeStats[statsCachedResult]);
                Console.WriteLine("branch cached:\t" + pairSizeStats[statsCachedBranch]);
                Console.WriteLine("single next kMer:\t" + pairSizeStats[statsSingleIdx]);
                Console.WriteLine("looked downstream:\t" + pairSizeStats[statsDownstreamIdx]);
                Console.WriteLine("single good downstream:\t" + pairSizeStats[statsSingleDSTIdx]);
                Console.WriteLine("chose in proportion by depth:\t" + pairSizeStats[statsRolledIdx]);
                Console.WriteLine("chose longest downstream:\t" + pairSizeStats[statsLongestIdx]);
                for (int pl = 0; pl < pairSizes.Count; pl++)
                    Console.WriteLine("single choice at length " + pairSizes[pl] + ":\t" + pairSizeStats[statsPairsIdx + pl]);
            }
        }

        private static void WriteFullUsage()
        {
            Console.WriteLine("usage: Kelpie [-h] -f forwardPrimer -r reversePrimer filteredReadsFNP extendedReadsFN (" + version + ")");
            Console.WriteLine("       -h                - Write out this extended help and exit");
            Console.WriteLine("       -f forwardPrimer e.g. GTGYCAGCMGCCGCGGTAA");
            Console.WriteLine("       -r reversePrimer e.g. GGACTACNVGGGTWTCTAAT");
            Console.WriteLine("       filteredReadsFNP  - List of file names/patterns of filtered reads to be processed. e.g. S1_270_reads_anonymous_R?_16S.fasta");
            Console.WriteLine("       extendedReadsFN   - File name for extended between-primer reads. e.g. S1_270_reads_anonymous_16S_EMB.fasta");
            Console.WriteLine();
            Console.WriteLine("uncommonly needed options that can usually be ignored");
            Console.WriteLine("       -paired|unpaired  - Use -paired if you have paired-end reads, and -unpaired if you want to force pairs of reads to be processed individually.");
            Console.WriteLine("                           Default is to assume pairs of reads are paired. Used when cleaning final kMer filter table.");
            Console.WriteLine("       -min nn           - Minimum length after assembly/extension. Default is that extended 'reads' must be full-length (from a starting primer to a terminating primer).");
            Console.WriteLine("       -max nn           - Drop extended reads longer than this length. A rarely useful option, but added to better handle improperly pair-merged reads found in some test datasets.");
            Console.WriteLine("       -mindepth nn      - kMers found fewer than this many times are dropped from the filtering table. Added to help pull out full length 'amplicons' from a very deep dataset.");
            Console.WriteLine("       -save primerTag   - Save the filtered/trimmed between-primer reads, and add this tag when building the file names. e.g. EMB --> S1_270_reads_anonymous_R?_16S_EMB.fasta");
            Console.WriteLine("       -log              - Debugging log from read extension phase.");
        }

        private static void FindReadsWithPrimers(List<string> headersToFilter, List<string> readsToFilter, int primerLength, HashSet<string> startingPrimers, HashSet<string> terminatingPrimers, HashSet<string> forwardPrimers,
                                                 HashSet<string> endingFilter, List<string> selectedHeaders, List<string> selectedReads, HashSet<int> readsAlreadyProcessed, int kMerSize)
        {
            // add reads to this list as they are selected - but can be dropped again if they contain Ns
            List<int> readsToSelect = new List<int>();

            HashSet<ulong> startingPrimersPacked = new HashSet<ulong>();
            HashSet<ulong> terminatingPrimersPacked = new HashSet<ulong>();
            foreach (string startingPrimer in startingPrimers)
                startingPrimersPacked.Add(kMers.CondenseMer(startingPrimer, primerLength));
            foreach (string terminatingPrimer in terminatingPrimers)
                terminatingPrimersPacked.Add(kMers.CondenseMer(terminatingPrimer, primerLength));

            for (int r = 0; r < readsToFilter.Count; r++)
            {
                string read = readsToFilter[r];

                //if (read == "TTATGCCTGAGGGTTCCACATCTGACTTGCCAAGCCGCCTACGAGCTCTTTACGCCCAATGAATCCGGACAACGCTTGCCCCCTACGTCTTACCGCGGCTGCTGGCACGTAGTTAGCCGGGGCTTCTTCGCCAGGTACAGTCATCTTCGT")
                //    Debugger.Break();

                int startingPrimerIdx = -1;
                int terminalPrimerIdx = -1;
                int lastPossiblePrimerIdx = read.Length - primerLength;
                bool previousPossiblePrimerKnown = false;
                ulong previousPossiblePrimer = 0;

                for (int i = 0; i <= lastPossiblePrimerIdx; i++)
                {
                    ulong possiblePrimerPacked;
                    if (previousPossiblePrimerKnown)
                        previousPossiblePrimerKnown = kMers.CondenseMerIncremental(primerLength, previousPossiblePrimer, read, i, out possiblePrimerPacked);
                    else
                        previousPossiblePrimerKnown = kMers.CondenseMer(read, i, primerLength, out possiblePrimerPacked);

                    previousPossiblePrimer = possiblePrimerPacked;

                    if (startingPrimersPacked.Contains(possiblePrimerPacked))
                    {
                        startingPrimerIdx = i;
                        break;
                    }
                    if (terminatingPrimersPacked.Contains(possiblePrimerPacked))
                    {
                        terminalPrimerIdx = i;
                        break;
                    }
                }

                // (+) ....F.....................(R).....
                // R1  ....F.....................(R).....
                // R2  ....R.....................(F).....
                //
                // (-) ....R.....................(F)..... (after reversal)
                // R1  ....R.....................(F).....
                // R2  ....F.....................(R).....

                // if we found a starting primer (F or R) in a read we'll trim the read and save it 
                if (startingPrimerIdx != -1)
                {
                    string possiblePrimer = read.Substring(startingPrimerIdx, primerLength);
                    // if the starting primer was a forward primer, mark the read (forwards direction)
                    if (forwardPrimers.Contains(possiblePrimer))
                        headersToFilter[r] += ";FP";
                    else
                        headersToFilter[r] += ";RP";

                    // trim away anything before the start of the primer (and leave the primer)
                    readsToFilter[r] = read.Substring(startingPrimerIdx);

                    if (read.Length >= kMerSize)
                    {
                        // save the selected read and remember that we've processed it 
                        readsToSelect.Add(r);
                        readsAlreadyProcessed.Add(r);
                    }

                    continue;
                }

                // if we found a terminal primer (F' or R') in a read, remember it and trim the read
                if (terminalPrimerIdx != -1)
                {
                    // add the kMer after the terminal primer to the endingFilter set
                    int startingIdx = terminalPrimerIdx + primerLength - kMerSize;
                    if (startingIdx > 0)
                    {
                        string terminalMer = read.Substring(startingIdx, kMerSize);
                        if (!endingFilter.Contains(terminalMer))
                            endingFilter.Add(terminalMer);
                    }

                    string possiblePrimer = read.Substring(terminalPrimerIdx, primerLength);
                    // mark the read if it contains a forward primer (RC form in this case as we know it's a terminating primer)
                    if (forwardPrimers.Contains(possiblePrimer))
                        headersToFilter[r] += ";FP'";
                    else
                        headersToFilter[r] += ";RP'";

                    // trim away anything after the primer (and leave the primer)
                    readsToFilter[r] = read.Substring(0, terminalPrimerIdx + primerLength);

                    if (read.Length >= kMerSize)
                    {
                        // save the selected read and remember that we've processed it 
                        readsToSelect.Add(r);
                        readsAlreadyProcessed.Add(r);
                    }
                }
            }

            // save the selected reads (unless they contain N bases)
            for (int r = 0; r < readsToSelect.Count; r++)
            {
                string header = headersToFilter[readsToSelect[r]];
                string read = readsToFilter[readsToSelect[r]];

                // add starting reads to the 'selected' list unless they contain Ns
                if (read.IndexOf('N') < 0)
                {
                    selectedHeaders.Add(header);
                    selectedReads.Add(read);
                }
            }
        }

        private static int ExtendFilterSet(List<string> selectedReads, HashSet<string> startingFilter, int kMerSize, HashSet<string> endingFilter)
        {
            int initialFilterCount = startingFilter.Count;

            int readsAddingToFilter = 0;

            for (int r = 0; r < selectedReads.Count; r++)
            {
                string read = selectedReads[r];
                //bool readWritten = false;
                bool readAddedToFilter = false;

                // tile for kMers starting anywhere in the read (as long as there's room for at least one kMer)
                int lastPossibleMerIdx = read.Length - kMerSize;
                int addedFromRead = 0;

                //if (read == "TTATGCCTGAGGGTTCCACATCTGACTTGCCAAGCCGCCTACGAGCTCTTTACGCCCAATGAATCCGGACAACGCTTGCCCCCTACGTCTTACCGCGGCTGCTGGCAC")
                //    Debugger.Break();

                // starting from the beginning of the read
                for (int i = 0; i <= lastPossibleMerIdx; i++)
                {
                    // get next kMer tile from this read
                    string kMer = read.Substring(i, kMerSize);

                    // or if the kMer contains an N (these can block terminating primers and don't want to include them anyway)
                    if (kMer.IndexOf('N') != -1)
                        break;

                    if (HomopolymerStart(kMer))
                        break;

                    // and if we haven't see it before, add it to the filter
                    if (!startingFilter.Contains(kMer))
                    {
                        startingFilter.Add(kMer);
                        addedFromRead++;

                        //if (kMer == "CCATCTTGTATTAATCTTCCTTTCAGAAGGC")
                        //    Debugger.Break();

                        if (!readAddedToFilter)
                        {
                            readAddedToFilter = true;
                            readsAddingToFilter++;
                        }
                    }

                    // stop tiling if we hit a terminating kMer
                    if (endingFilter.Contains(kMer))
                        break;
                }
            }

            int kMersAddedToFilter = startingFilter.Count - initialFilterCount;
            //Console.WriteLine(readsAddingToFilter + " reads added " + kMersAddedToFilter + " kMers to starting filter");
            return readsAddingToFilter;
        }

        private static bool HomopolymerStart(string kMer)
        {
            string startOfMer = kMer.Substring(0, 10);

            return startOfMer == "AAAAAAAAAA" || startOfMer == "CCCCCCCCCC" || startOfMer == "GGGGGGGGGG" || startOfMer == "TTTTTTTTTT";
        }

        // filter reads against the agreed, final set of kMers
        private static void FinalFilterReads(List<string> headersToFilter, List<string> readsToFilter, int kMerSize, HashSet<string> kMerFilter, List<string> selectedHeaders, List<string> selectedReads)
        {
            for (int r = 0; r < readsToFilter.Count; r++)
            {
                //if (headersToFilter[r].StartsWith(">RM2|S1|R41203241/1"))
                //    Debugger.Break();

                string header = headersToFilter[r];
                string read = readsToFilter[r];
                //int lastSemiIdx = header.LastIndexOf(';');
                //if (lastSemiIdx > 0)
                //{
                //    string tag = header.Substring(lastSemiIdx);
                //    if (tag == ";FP" || tag == ";FP'")              //  || tag == ";RP" || tag == ";RP'")
                //    {
                //        // save all potential starting reads
                //        selectedHeaders.Add(headersToFilter[r]);
                //        selectedReads.Add(read);
                //        continue;
                //    }
                //}

                // look for a match towards the start of a read
                int lastStartIdx = Math.Min(startSize, read.Length - kMerSize);
                for (int i = 0; i < lastStartIdx; i++)
                {
                    string startingMer = read.Substring(i, kMerSize);
                    if (kMerFilter.Contains(startingMer))
                    {
                        // save the selected read 
                        selectedHeaders.Add(headersToFilter[r]);
                        selectedReads.Add(read);
                        break;
                    }
                }
            }
        }

        // find reads that match against the (newly-updated) filter kMers, these reads will then be tiled for more kMers to add to the filter
        private static List<string> FindMatchingReads(List<string> readsToFilter, HashSet<int> readsAlreadyProcessed, int kMerSize, HashSet<string> kMerFilter)
        {
            List<string> matchingReads = new List<string>();

            for (int r = 0; r < readsToFilter.Count; r++)
            {
                // skip over reads that have already been processed
                if (readsAlreadyProcessed.Contains(r))
                    continue;

                bool breakOnRead = false;
                //if (readsToFilter[r] == "GAACCCCGATTGCGAAGGCAGCATACTGGGCTATAACTGACGCTGAAGCACGAAAGCGTGGGTATCGAACAGGATTAGATACCCTGGTAGTC")
                //    //breakOnRead = true;
                //    Debugger.Break();

                // look for a match at the start of a read
                string read = readsToFilter[r];
                if (read.Length < kMerSize)
                    continue;

                string startingMer = read.Substring(0, kMerSize);
                if (kMerFilter.Contains(startingMer))
                {
                    // save the selected read and remember not to process it again
                    matchingReads.Add(readsToFilter[r]);
                    readsAlreadyProcessed.Add(r);
                    //if (breakOnRead)
                    //    Debugger.Break();
                }
            }

            //int selectedReadsLengthSum = 0;
            //foreach (string selectedRead in matchingReads)
            //    selectedReadsLengthSum += selectedRead.Length;
            //Console.WriteLine("selected " + matchingReads.Count + " reads (" + selectedReadsLengthSum + " bases) from " + (readsToFilter.Count - readsAlreadyProcessed.Count));

            return matchingReads;
        }

        private static int WriteFilteredReads(string readsFN, string suffix, int firstReadToWrite, List<string> selectedHeaders, List<string> selectedReads, int minLength)
        {
            int readsWritten = 0;
            string FNWithoutSuffix = readsFN;
            string fileType = ".fa";
            int lastDot = readsFN.LastIndexOf('.');
            if (lastDot > 0)
            {
                FNWithoutSuffix = readsFN.Substring(0, lastDot);
                fileType = readsFN.Substring(lastDot);
            }

            string filteredReadsFN = FNWithoutSuffix + "_" + suffix + fileType;
            StreamWriter filteredReads = new StreamWriter(filteredReadsFN);

            for (int i = firstReadToWrite; i < selectedHeaders.Count; i++)
            {
                string read = selectedReads[i];

                if (read.Length < minLength)
                    continue;

                string header = selectedHeaders[i];
                if (header[0] == '@')
                    header = ">" + header.Substring(1);
                SeqFiles.WriteRead(filteredReads, header, read, SeqFiles.formatFNA);
                readsWritten++;
            }

            filteredReads.Close();
            Console.WriteLine("wrote " + readsWritten + " reads to " + filteredReadsFN);
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

        private static void AddMerToTable(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize)
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

        private static ulong HashPair(ulong kMer1, ulong kMer2, int kMerSize)
        {
            return Canonical(Canonical(kMer1, kMerSize) ^ Canonical(kMer2, kMerSize), kMerSize);
        }

        private static int DeNoiseMerTable(List<string> selectedReads, int kMerSize, Dictionary<ulong, int> kMerTable, int minDepth, out int noiseLevel, out int meanDepth, StreamWriter log)
        {
            ulong[] kMerVariants = new ulong[4];
            int[] kMerCounts = new int[4];

            // compute a mean depth (excluding singletons etc)
            long nonSingletons = 0;
            Dictionary<int, int> depthHistogram = new Dictionary<int, int>();

            // use harmonic mean rather than average to reduce the impact of conserved regions and dominant organisms
            double invDepthSum = 0.0;
            long totalCount = 0;
            foreach (int count in kMerTable.Values)
            {
                if (count >= minDepth)
                {
                    invDepthSum += 1.0f / (double)count;
                    nonSingletons++;
                    totalCount += count;
                }
                if (depthHistogram.ContainsKey(count))
                    depthHistogram[count]++;
                else
                    depthHistogram.Add(count, 1);
            }

            meanDepth = (int)((double)nonSingletons / invDepthSum);
            long avgDepth = totalCount / nonSingletons;
            noiseLevel = Math.Max(minDepth, meanDepth / 10);    // 10% of mean
            meanDepth = Math.Max(meanDepth, noiseLevel * 5);    // in case mean is too low

            int kMersCulled = 0;

            // scan through the reads again looking for possible sequencng errors
            // triggers are sudden drops down to noise level. Low depth reads are OK
            int previousDepth = -1;

            foreach (string read in selectedReads)
            {
                //if (read == "CCCTCACGGGGGTGAGCCCCGCACTTTTAAGAAAGACTTACCGTCCCACCTGCGCTCCCTTTACGCCCAATGATTCCGGACAACGCTTGCCGCTTACGTATTACCGCGGCTGCTGGCAC")
                //    Debugger.Break();

                ulong[] kMersFromRead = new ulong[1000];
                bool[] kMersValid = new bool[1000];
                int kMerCount = kMers.GenerateMersFromRead(read, kMerSize, ref kMersFromRead, ref kMersValid);

                for (int m = 0; m < kMerCount; m++)
                {
                    ulong kMer = kMersFromRead[m];
                    if (!kMersValid[m])
                        continue;

                    int currentDepth = kMerTableLookup(kMer, kMerTable, kMerSize);
                    if (previousDepth < 0)
                        previousDepth = currentDepth;

                    // skip if kMer has been culled already (propagation)
                    if (currentDepth == 0)
                        continue;
                    // nothing much has changed... 
                    if (Close(previousDepth, currentDepth))
                        continue;
                    // and we haven't gone up rather than down
                    if (currentDepth > previousDepth)
                        continue;

                    // check the following kMer to see if the depth remains unchanged (indicating that we have a low-frequency variant rather than an error)
                    ulong kMerFollowing;
                    int idxOfFollowingMer = m + kMerSize;
                    int followingDepth = -1;
                    if (idxOfFollowingMer  < kMerCount)
                    {
                        kMers.CondenseMer(read, idxOfFollowingMer, kMerSize, out kMerFollowing);
                        followingDepth = kMerTableLookup(kMerFollowing, kMerTable, kMerSize);
                    }
                    if (Close(currentDepth, followingDepth))
                        continue;

                    //if (kMer == 0xAEC9AE03992334A3)
                    //    Debugger.Break();

                    // get counts for all variants of this kMer (last base only), including the kMer itself
                    GetMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerCounts);

                    int deepestDepth = 0;       // depth of deepest variant
                    ulong deepestMer = 0;       // the deepest variant kMer (could be the current kMer)
                    int nonZeroCount = 0;       // how many non-zero depth variants were found

                    for (int b = 0; b < kMerVariants.Length; b++)
                    {
                        int variantDepth = kMerCounts[b];
                        if (variantDepth > 0)
                            nonZeroCount++;

                        if (variantDepth > deepestDepth)
                        {
                            deepestDepth = variantDepth;
                            deepestMer = kMerVariants[b];
                        }
                    }

                    // erase the count of erroneous kMers and increment the depth of the highest depth alternative
                    if (currentDepth < noiseLevel ||
                        (deepestDepth > currentDepth &&
                         currentDepth < deepestDepth / 100 &&
                         currentDepth < meanDepth / 2))                                  // small (1%) and deeper alternative present
                    {
                        // shift the count to the deepest alternative
                        bool countTransferred = false;
                        if (deepestMer != kMer && deepestDepth > currentDepth)
                        {
                            kMerTable[Canonical(deepestMer, kMerSize)] += currentDepth;
                            countTransferred = true;
                        }

                        // leave the kMer in the table but zero its count
                        kMerTable[Canonical(kMer, kMerSize)] = 0;
                        kMersCulled++;

                        if (log != null)
                        {
                            log.Write("kdr: " + kMers.ExpandMer(kMer, kMerSize) + " [" + kMers.ExpandMer(Canonical(kMer, kMerSize), kMerSize) + "] (" + currentDepth + ") ");
                            for (int b = 0; b < kMerCounts.Length; b++)
                                log.Write(kMerCounts[b] + "  ");
                            if (countTransferred)
                                log.WriteLine("to " + kMers.ExpandMer(deepestMer, kMerSize) + " " + kMerTable[Canonical(deepestMer, kMerSize)]);
                            else
                                log.WriteLine();
                        }

                        // propagate cull downstream
                        for (m++; m < kMerCount; m++)
                        {
                            ulong kMerDownstream;
                            kMers.CondenseMer(read, m, kMerSize, out kMerDownstream);
                            int downstreamDepth = kMerTableLookup(kMerDownstream, kMerTable, kMerSize);

                            if (downstreamDepth == 0)
                                continue;

                            if (downstreamDepth <= currentDepth)
                            {
                                GetMerVariantsAndCounts(kMerDownstream, kMerTable, kMerSize, kMerVariants, kMerCounts);
                                deepestDepth = 0;
                                nonZeroCount = 0;

                                for (int b = 0; b < kMerVariants.Length; b++)
                                {
                                    if (kMerCounts[b] > deepestDepth)
                                    {
                                        deepestDepth = kMerCounts[b];
                                        deepestMer = kMerVariants[b];
                                        nonZeroCount++;
                                    }
                                }

                                kMerTable[Canonical(kMerDownstream, kMerSize)] = 0;
                                kMersCulled++;

                                if (log != null)
                                    log.WriteLine("kdp: " + kMers.ExpandMer(kMerDownstream, kMerSize) + " [" + kMers.ExpandMer(Canonical(kMerDownstream, kMerSize), kMerSize) + "] (" + downstreamDepth + ") ");

                                // deeper alternative so transfer the count and stop propagating
                                if (deepestDepth > downstreamDepth)
                                {
                                    kMerTable[Canonical(deepestMer, kMerSize)] += downstreamDepth;
                                    break;
                                }
                            }
                            else
                                // no longer propagating as we've hit something higher
                                break;
                        }

                        // if we terminated the propagation early, we'll restart the outer scanning loop at the point of termination
                        if (m < kMerCount)
                            m--;
                    }
                }
            }

            // recalculate the noise level now that we've cleaned up the kMer table
            invDepthSum = 0.0;
            nonSingletons = 0;
            totalCount = 0;
            depthHistogram.Clear();
            foreach (int count in kMerTable.Values)
            {
                if (count > minDepth)
                {
                    invDepthSum += 1.0f / (double)count;
                    nonSingletons++;
                    totalCount += count;
                }
                if (depthHistogram.ContainsKey(count))
                    depthHistogram[count]++;
                else
                    depthHistogram.Add(count, 1);
            }

            avgDepth = totalCount / nonSingletons;
            meanDepth = (int)((double)nonSingletons / invDepthSum);
            noiseLevel = Math.Max(minDepth, meanDepth / 10);    // 10% of mean
            meanDepth = Math.Max(meanDepth, noiseLevel * 5);    // in case mean is too low

            return kMersCulled;
        }

        private static void GetMerVariantsAndCounts(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize, ulong[] kMerVariants, int[] kMerCounts)
        {
            ulong mask = 0xffffffffffffffff << (32 - kMerSize + 1) * 2;
            ulong kMerPrefix = kMer & mask;

            for (int b = 0; b < 4; b++)
            {
                ulong kMerVariant = kMerPrefix | ((ulong)b << (32 - kMerSize) * 2);
                kMerCounts[b] = kMerTableLookup(kMerVariant, kMerTable, kMerSize);
                kMerVariants[b] = kMerVariant;
            }
        }

        private static bool CleanStartingRead(int kMerSize, Dictionary<ulong, int> kMerTable, List<int> pairSizes, Dictionary<int, Dictionary<ulong, int>> kMerPairs, int noiseLevel, int meanDepth, Sequence seq, StreamWriter log)
        {
            int consecutiveChanges = 0;
            bool readClean = true;
            bool readChanged = false;
            ulong[] kMerVariants = new ulong[4];
            int[] kMerCounts = new int[4];
            List<ulong> kMerVariantsList = new List<ulong>(500);
            List<int> kMerCountsList = new List<int>(500);
            ulong[] kMersFromRead = new ulong[1000];
            bool[] kMersValid = new bool[1000];
            int previousDepth = -1;

            string originalRead = (log == null) ? null : seq.ToString();

            //if (seq.ToString() == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTTGCGTCGGTCGTGAAAACCTGGGGCTTAACCCTGGGCTTGCGGTCGATACGGGCAGACTTGAG")
            //    Debugger.Break();

            // scan the read looking for aberrant kMers
            for (int i = 0; i < seq.Length - kMerSize + 1; i++)
            {
                kMerVariantsList.Clear();
                kMerCountsList.Clear();
                
                // look for kMers that appear to be suspect - poor kMer depth or no pair (when possible)
                ulong kMer;
                Sequence.CondenseMer(seq, i, kMerSize, out kMer);
                int currentDepth = kMerTableLookup(kMer, kMerTable, kMerSize);
                if (previousDepth < 0)
                    previousDepth = currentDepth;

                // does this kMer look OK? (above noise level, not too much lower than previous)
                if (currentDepth >= noiseLevel && Close(currentDepth,previousDepth))
                {
                    consecutiveChanges = 0;
                    previousDepth = currentDepth;
                    continue;
                }

                // possibly suspect kMer found... so generate variants and get their counts
                if (i == 0)
                {
                    // check all double-base variants for the first kMer
                    HashSet<string> doubleSubVariants = new HashSet<string>();
                    kMers.GenerateSeqVariants(kMers.ExpandMer(kMer, kMerSize), kMerSize, doubleSubVariants, 2);
                    foreach (string kMerVariant in doubleSubVariants)
                    {
                        ulong kmv = kMers.CondenseMer(kMerVariant, kMerSize);
                        kMerVariantsList.Add(kmv);
                        kMerCountsList.Add(kMerTableLookup(kmv, kMerTable, kMerSize));
                    }
                }
                else
                {
                    // just vary the last base for all other kMers
                    GetMerVariantsAndCounts(kMer, kMerTable, kMerSize, kMerVariants, kMerCounts);

                    foreach (ulong kMerVariant in kMerVariants)
                        kMerVariantsList.Add(kMerVariant);
                    foreach (int kMerCount in kMerCounts)
                        kMerCountsList.Add(kMerCount);
                }

                // scan the just-created variants to find the deepest and remember a few things about the variants
                int deepestDepth = 0;
                int deepestMerIdx = 0;
                ulong deepestMer = 0;
                int currentMerIdx = 0;
                for (int v = 0; v < kMerCountsList.Count; v++)
                {
                    int variantDepth = kMerCountsList[v];
                    if (variantDepth > deepestDepth)
                    {
                        deepestDepth = variantDepth;
                        deepestMer = kMerVariantsList[v];
                        deepestMerIdx = v;
                    }
                    if (kMerVariantsList[v] == kMer)
                        currentMerIdx = v;
                }

                // check all the still-viable variants for valid backwards pairs
                for (int v = 0; v < kMerVariantsList.Count; v++)
                {
                    bool pairsOKForVariant = false;
                    if (kMerCountsList[v] > 0)
                        pairsOKForVariant = CheckAllPairs(kMerVariantsList[v], seq, i, kMerSize, pairSizes, kMerPairs);

                    // cull any variants where a pair was not found
                    if (!pairsOKForVariant)
                        kMerCountsList[v] = 0;
                }

                // update current depth in case we have set it to zero
                currentDepth = kMerCountsList[currentMerIdx];

                // check following kMer from this potential repair to see if the depth remains unchanged (indicating that the repair is possibly breaking a low-frequency variant)
                ulong kMerFollowing;
                int idxOfFollowingMer = i + kMerSize;
                int followingDepth = -1;
                if (idxOfFollowingMer + kMerSize < seq.Length)
                {
                    Sequence.CondenseMer(seq, idxOfFollowingMer, kMerSize, out kMerFollowing);
                    followingDepth = kMerTableLookup(kMerFollowing, kMerTable, kMerSize);
                }

                // correct the kMer if it looks like it's likely to be an error coming from a dominant organism
                // never do this though if the count is anywhere near average, or close to the previous count
                if (kMer != deepestMer && ((currentDepth < 5 * deepestDepth / 100) || currentDepth == 0) 
                                       && !Close(currentDepth, previousDepth)
                                       && !Close(currentDepth, followingDepth)
                                       && currentDepth < meanDepth / 2
                                       && deepestDepth >= noiseLevel)
                {
                    seq.Replace(i, deepestMer, kMerSize);
                    currentDepth = deepestDepth;
                    consecutiveChanges++;
                    readChanged = true;
                }

                previousDepth = currentDepth;

                if (currentDepth == 0 || consecutiveChanges > 2)
                {
                    // don't truncate read if we can't fix the first kMer (containing the primer), just report it and it will be dropped from the starting reads later
                    if (i == 0)
                        readClean = false;
                    else
                    {
                        // truncate the read at this first irrecoverable poor kMer
                        seq.Length = i + kMerSize - 1 - consecutiveChanges;
                        readChanged = true;
                    }
                    break;
                }
            }

            if (log != null)
            {
                if (readChanged)
                {
                    string trimmed = (originalRead.Length > seq.Length) ? "/trimmed" : "";
                    log.WriteLine("cleaned" + trimmed + ": " + originalRead + " --> " + seq.ToString());
                }
                else
                { 
                    if (readClean)
                        log.WriteLine("clean: " + originalRead + " --> " + seq.ToString());
                    else
                        log.WriteLine("unclean: " + originalRead + " --> " + seq.ToString());
                }
            }

            return readClean;
        }

        private static bool Close(int currentDepth, int previousDepth)
        {
            int lowestDepth = Math.Min(currentDepth, previousDepth);
            int highestDepth = Math.Max(currentDepth, previousDepth);
            if (lowestDepth <= 0)
                return false;
            return (highestDepth - lowestDepth) <= highestDepth / 2;
        }

        private static bool CheckAllPairs(ulong kMer, Sequence seq, int ki, int kMerSize, List<int> pairSizes, Dictionary<int, Dictionary<ulong, int>> kMerPairs)
        {
            bool pairsOK = true;
            int pli = 0;
            while (pli < pairSizes.Count && pairsOK)
            {
                int pairSize = pairSizes[pli];
                int pairStart = ki + kMerSize - pairSize;
                if (pairStart < 0)
                    break;

                Dictionary<ulong, int> kMerPairsTable = kMerPairs[pairSize];

                ulong firstMerInPair;
                Sequence.CondenseMer(seq, pairStart, kMerSize, out firstMerInPair);

                int pairDepth = kMerTableLookup(HashPair(firstMerInPair, kMer, kMerSize), kMerPairsTable, kMerSize);

                pairsOK = pairDepth > 0;
                pli++;
            }

            return pairsOK;
        }

        private static bool ExtendShortStartingRead(Sequence read, Dictionary<ulong, int> kMerTable, int kMerSize, int minLengthRequired, List<int> pairSizes, Dictionary<int, Dictionary<ulong, int>> kMerPairs, int noiseLevel, int meanDepth, StreamWriter log)
        {
            // extend this read until we reach a fork, can't find a following kMer or have extended the read enough
            bool extending = true;
            int initialLength = read.Length;

            //if (read.ToString() == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTA")
            //    Debugger.Break();

            TrimTrailingPoorBases(read, kMerTable, kMerSize, log);

            if (read.Length < kMerSize)
                return false;

            ulong finalMerInRead;
            Sequence.CondenseMer(read, read.Length - kMerSize, kMerSize, out finalMerInRead);
            int previousCount = kMerTableLookup(finalMerInRead, kMerTable, kMerSize);

            ulong[] kMerVariants = new ulong[4];
            int[] kMerCounts = new int[4];

            while (extending)
            {
                if (read.Length >= minLengthRequired)
                    break;

                // look up each of the possible next kMers and get their counts
                GetMerVariantsAndCounts(finalMerInRead << 2, kMerTable, kMerSize, kMerVariants, kMerCounts);

                // validate the just-created variants
                int viableMers = 0;
                int deepestViableBaseIdx = 0;
                int deepestCount = 0;
                ulong deepestMer = 0;
                int totalCount = 0;

                // find the deepest (before cleaning)
                for (int v = 0; v < kMerCounts.Length; v++)
                {
                    int count = kMerCounts[v];
                    if (count > deepestCount)
                        deepestCount = count;
                }

                // clean variants
                for (int v = 0; v < kMerCounts.Length; v++)
                {
                    int count = kMerCounts[v];
                    if (count < noiseLevel &&
                        deepestCount > count &&
                        count < 5 * deepestCount / 100 &&
                        count < meanDepth / 2)
                    {
                        kMerCounts[v] = 0;
                    }
                    bool pairsOKForVariant = false;
                    if (kMerCounts[v] > 0)
                        pairsOKForVariant = CheckAllPairs(kMerVariants[v], read, (read.Length-kMerSize+1), kMerSize, pairSizes, kMerPairs);

                    // cull any variants where a pair was not found
                    if (!pairsOKForVariant)
                        kMerCounts[v] = 0;
                }

                // remember the deepest from the still-viable variants
                deepestCount = 0;
                for (int v = 0; v < kMerCounts.Length; v++)
                {
                    int count = kMerCounts[v];
                    if (count > 0)
                    {
                        viableMers++;
                        totalCount += count;
                        if (count > deepestCount)
                        {
                            deepestCount = count;
                            deepestViableBaseIdx = v;
                            deepestMer = kMerVariants[v];
                        }
                    }
                }

                // extend read if one base represents 98%+ of all of variants or we're just continuing
                if ((totalCount > 0 && deepestCount * 1000 / totalCount > 980) || (deepestCount >= noiseLevel && Close(deepestCount, previousCount)))
                {
                    read.Append(bases[deepestViableBaseIdx]);
                    finalMerInRead = deepestMer;
                    previousCount = deepestCount;
                }
                else
                    extending = false;
            }

            return read.Length > initialLength;
        }

        private static void TrimTrailingPoorBases(Sequence read, Dictionary<ulong, int> kMerTable, int kMerSize, StreamWriter log)
        {
            bool stillTrimming = true;
            bool readWasTrimmed = false;
            int basesTrimmed = 0;
            string startingRead = log == null ? null : read.ToString();

            while (stillTrimming)
            {
                if (read.Length < kMerSize)
                    break;

                ulong finalMerInRead;
                Sequence.CondenseMer(read, read.Length - kMerSize, kMerSize, out finalMerInRead);

                int count = kMerTableLookup(finalMerInRead, kMerTable, kMerSize);

                if (count > 1)
                    break;

                read.Length--;
                readWasTrimmed = true;
                basesTrimmed++;
            }

            if (readWasTrimmed && log != null)
            {
                log.WriteLine("trimmed: " + startingRead + " --> " + read.ToString());
            }
        }

        private static Dictionary<string, int> TrimExtendedReads(Dictionary<Sequence, int> extendedReads, int forwardPrimerLength, int reversePrimerLength, Dictionary<string, int> forwardPrimersFound, HashSet<ulong> terminatingPrimers, int minLength, int maxLength, string direction, StreamWriter log)
        {
            Dictionary<string, int> trimmedReads = new Dictionary<string, int>(extendedReads.Count);
            int droppedReadCount = 0;

            foreach (KeyValuePair<Sequence, int> kvp in extendedReads)
            {
                string extendedRead = kvp.Key.ToString();
                int extendedReadCount = kvp.Value;
                bool foundTerminatingPrimer = false;

                //if (extendedRead == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTTGCGTCGGTCGTGAAAACCTGGGGCTTAACCCTGGGCTTGCGGTCGATACGGGCAGACTTGAGTTCGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGTGGCGAAGGCGGGTCTCTGGGAAACAACTGACGCTGAGGAACGAAAGCGTGGGTAGCAAACAGGATTAGATACCCTGGTAGTC")
                //    Debugger.Break();
                //if (extendedRead.Contains("TACGGAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGACCGGAAAGTTGGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCTTCAAAACTATCGGTCTGGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAAATCGTAGATATTCGGAAGAACACCAGTGGCGAAGGCGGCTCACTGGACGACTGTTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"))
                //    Debugger.Break();

                // trim the forward primer at the start
                string possibleFP = extendedRead.Substring(0, forwardPrimerLength);
                if (forwardPrimersFound.ContainsKey(possibleFP))
                    extendedRead = extendedRead.Substring(forwardPrimerLength);

                // check for a terminal primer at the end
                int lastPossiblePrimerIdx = extendedRead.Length - reversePrimerLength;
                string possiblePrimer = extendedRead.Substring(lastPossiblePrimerIdx, reversePrimerLength);
                ulong possiblePrimerPacked = kMers.CondenseMer(possiblePrimer, reversePrimerLength);
                if (terminatingPrimers.Contains(possiblePrimerPacked))
                {
                    extendedRead = extendedRead.Substring(0, lastPossiblePrimerIdx);
                    foundTerminatingPrimer = true;
                }

                if (extendedRead.Length <= maxLength && (foundTerminatingPrimer || (minLength > 0 && extendedRead.Length >= minLength)))
                {
                    if (trimmedReads.ContainsKey(extendedRead))
                        trimmedReads[extendedRead] += extendedReadCount;
                    else
                        trimmedReads.Add(extendedRead, extendedReadCount);
                }
                else
                {
                    droppedReadCount += extendedReadCount;
                    if (log != null)
                        log.WriteLine("dropped: " + extendedRead + " " + (foundTerminatingPrimer ? "TP" : "NTP") + ", " + extendedRead.Length + "b");
                }
            }

            Console.WriteLine("dropped " + droppedReadCount + " short " + direction + " reads");
            return trimmedReads;
        }

        private static Dictionary<Sequence, int> ExtendFromStartingReads(List<Sequence> startingReads, int kMerSize, Dictionary<ulong, int> kMerTable,
                                                              List<int> pairSizes, Dictionary<int, Dictionary<ulong, int>> kMerPairs, int[] pairSizeStats,
                                                              HashSet<ulong> terminatingPrimers, int primerSize, int noiseLevel, int meanDepth, StreamWriter log)
        {
            Dictionary<Sequence, int> fullyExtendedReads = new Dictionary<Sequence, int>();         // reads that have been extended as far as possible (and will be returned from this method) (reads and counts)
            Dictionary<string, Sequence> cachedExtensions = new Dictionary<string, Sequence>();     // read prefixes that have been extended previously - no need to repeat
            Dictionary<string, bool> cachedReachedTerminalPrimer = new Dictionary<string, bool>();  // and remember if this extension reached a terminal primer

            HashSet<ulong> loopTrap = new HashSet<ulong>();

            Stopwatch sw = new Stopwatch();

            // extending the 'starting' reads
            foreach (Sequence read in startingReads)
            {
                string startingRead = read.ToString();

                if (log != null)
                {
                    log.WriteLine("extending: " + startingRead);
                    sw.Restart();
                }

                bool breakOnRead = false;
                //if (startingRead == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAA")
                //{
                //    Debugger.Break();
                //    breakOnRead = true;
                //}

                if (cachedExtensions.ContainsKey(startingRead))
                {
                    Sequence cachedExtendedRead = cachedExtensions[startingRead];
                    if (fullyExtendedReads.ContainsKey(cachedExtendedRead))
                        fullyExtendedReads[cachedExtendedRead]++;
                    else
                        fullyExtendedReads.Add(cachedExtendedRead, 1);
                    pairSizeStats[statsCachedResult]++;
                    if (log != null)
                    {
                        sw.Stop();
                        log.WriteLine("cached:    " + cachedExtensions[startingRead].ToString()); // + " " + sw.Elapsed.TotalSeconds.ToString("f3"));
                    }
                    continue;
                }

                loopTrap.Clear();
                bool terminatingPrimerReached;
                bool coinTossed = false;
                //double averageDepth = AverageMerDepth(read, kMerSize, kMerTable);
                int previousDepth = -1;
                Sequence extendedRead = ExtendRead(read, kMerSize, kMerTable, pairSizes, kMerPairs, pairSizeStats, previousDepth, /*averageDepth,*/
                                                   cachedExtensions, cachedReachedTerminalPrimer, loopTrap, terminatingPrimers,
                                                   primerSize, noiseLevel, meanDepth, out terminatingPrimerReached, out coinTossed, log, 1, breakOnRead);
                if (fullyExtendedReads.ContainsKey(extendedRead))
                    fullyExtendedReads[extendedRead]++;
                else
                    fullyExtendedReads.Add(extendedRead, 1);
                //if (extendedRead.ToString() == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTTGCGTCGGTCGTGAAAACCTGGGGCTTAACCCTGGGCTTGCGGTCGATACGGGCAGACTTGAGTTCGGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCCACTACTGACGCTGAGAAGCGAAAGCATGGGGAGCGAACAGGATTAGATACCCTGGTAGTC")
                //    Debugger.Break();

                if (cacheThings && !coinTossed)
                {
                    cachedExtensions.Add(startingRead, extendedRead);
                    cachedReachedTerminalPrimer.Add(startingRead, terminatingPrimerReached);
                }

                if (log != null)
                {
                    sw.Stop();
                    log.WriteLine("extended:  " + extendedRead.ToString() + " (" /*+ sw.Elapsed.TotalSeconds.ToString("f3")*/ + ") ct=" + coinTossed);
                }
            }

            return fullyExtendedReads;
        }

        private static double AverageMerDepth(Sequence read, int kMerSize, Dictionary<ulong, int> kMerTable)
        {
            ulong[] kMersInRead = new ulong[read.Length - kMerSize + 1];
            bool[] kMersValid = new bool[kMersInRead.Length];

            Sequence.GenerateMersFromRead(read, kMerSize, ref kMersInRead, ref kMersValid);

            int depthSum = 0;
            int validCount = 0;
            for (int i = 0; i < kMersInRead.Length; i++)
                if (kMersValid[i])
                {
                    depthSum += kMerTableLookup(kMersInRead[i], kMerTable, kMerSize);
                    validCount++;
                }

            if (validCount > 0)
                return (double)depthSum / (double)validCount;
            else
                return 0;
        }

        private static Sequence ExtendRead(Sequence read, int kMerSize, Dictionary<ulong, int> kMerTable,
                                           List<int> pairSizes, Dictionary<int, Dictionary<ulong, int>> kMerPairs, int[] pairSizeStats, int previousDepth,  /*double averageDepth,*/
                                           Dictionary<string, Sequence> cachedExtensions, Dictionary<string, bool> cachedReachedTerminalPrimer,
                                           HashSet<ulong> loopTrap, HashSet<ulong> terminatingPrimers, int primerSize, int noiseLevel, int meanDepth,
                                           out bool terminatingPrimerReached, out bool tossedCoin, StreamWriter log, int level, bool breakOnRead)
        {
            bool traceRead = false;

            // keep extending until there are no following kMers 
            bool extending = true;
            terminatingPrimerReached = false;
            tossedCoin = false;
            int kMersInRead = read.Length - kMerSize + 1;
            int basesAdded = 0;

            ulong kMer;
            Sequence.CondenseMer(read, read.Length - kMerSize, kMerSize, out kMer);

            //if (read.ToString() == "GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTGTCCGGAATTATTGGGCGTAA")
            //{
            //    breakOnRead = true;
            //    Debugger.Break();
            //}

            if (loopTrap.Contains(kMer))
            {
                //if (log != null)
                //    log.WriteLine(new string(' ', level * 4) + level + ": ER loop (r) with " + kMers.ExpandMer(kMer, kMerSize) + " in " + read.ToString());
                terminatingPrimerReached = false;
                return read;
            }

            if (kMerIsTerminating(kMer, kMerSize, terminatingPrimers, primerSize))
            {
                terminatingPrimerReached = true;
                return read;
            }

            if (previousDepth == -1)
                previousDepth = kMerTableLookup(kMer, kMerTable, kMerSize);

            ulong[] kMerVariants = new ulong[4];
            int[] kMerCounts = new int[4];
            bool[] variantAlive = new bool[4];
            bool[] variantPreviouslyAlive = new bool[4];

            //if (log != null)
            //    log.WriteLine(new string(' ', level * 4) + level + ": " + "extending " + read.ToString());

            while (extending)
            {
                int viableAlternatives = 0;
                int chosenIdx = -1;
                pairSizeStats[statsNextIdx]++;

                //if (kMer == 0x2B8BA0F58BB22B80 /*&& breakOnRead*/)
                //    Debugger.Break();

                // generate all 4 possible following kMers 
                GetMerVariantsAndCounts(kMer << 2, kMerTable, kMerSize, kMerVariants, kMerCounts);

                // and stop if none of the extensions are at all viable
                for (int b = 0; b < 4; b++)
                    if (kMerCounts[b] > 0)
                    {
                        viableAlternatives++;
                        chosenIdx = b;
                        variantAlive[b] = true;
                    }
                    else
                        variantAlive[b] = false;
                if (viableAlternatives == 0)
                    break;

                // denoise alternatives (initial denoise only checked as-read kMers, and we could now be going in the opposite direction)
                // this is done primarily for performance, stopping the search of a large tree generated through these not-denoised kMer variants.
                if (viableAlternatives > 1)
                {
                    int deepestCount = 0;

                    // find the deepest variant
                    for (int b = 0; b < 4; b++)
                    {
                        int variantCount = kMerCounts[b];
                        if (variantCount > deepestCount)
                        {
                            deepestCount = variantCount;
                            chosenIdx = b;
                        }
                    }

                    // check whether any other of the other variants seem to be error variants from this deepest kMer
                    for (int b = 0; b < 4; b++)
                    {
                        if (!variantAlive[b])
                            continue;
                        int currentDepth = kMerCounts[b];

                        if (currentDepth > 0 &&
                            deepestCount > currentDepth &&
                            currentDepth < 5 * deepestCount / 100 &&
                            currentDepth < meanDepth / 2)                                  // small (5%) and deeper alternative present
                        {
                            if (variantAlive[b])
                            {
                                variantAlive[b] = false;
                                viableAlternatives--;
                            }
                        }
                    }
                }

                if (viableAlternatives == 1)
                    pairSizeStats[statsSingleIdx]++;             // only one next kMer 

                // if there were multiple next kMers, iterate through the pair lengths, until we run out of pair support, get just one survivor, or get to the maximum pair length
                int pli = 0;
                while (viableAlternatives > 1 && pli < pairSizes.Count)
                {
                    int pairSize = pairSizes[pli];
                    if (read.Length < pairSize)
                        break;

                    Dictionary<ulong, int> kMerPairsTable = kMerPairs[pairSize];
                    int previousViableAlternatives = viableAlternatives;
                    viableAlternatives = 0;

                    ulong initialMerInPair;
                    // start of pair allows for the base that is being added (hence the +1)
                    Sequence.CondenseMer(read, read.Length - pairSize + 1, kMerSize, out initialMerInPair);
                    int highestPairDepth = 0;

                    for (int b = 0; b < 4; b++)
                    {
                        variantPreviouslyAlive[b] = variantAlive[b];

                        if (variantAlive[b])
                        {
                            ulong hashedPair = HashPair(initialMerInPair, kMerVariants[b], kMerSize);
                            int pairDepth = kMerTableLookup(hashedPair, kMerPairsTable, kMerSize);

                            if (pairDepth > 0)
                            {
                                viableAlternatives++;
                                chosenIdx = b;
                                variantAlive[b] = true;
                                if (pairDepth > highestPairDepth)
                                    highestPairDepth = pairDepth;
                            }
                            else
                                variantAlive[b] = false;
                        }
                    }

                    // stop looking for pair matches if we've run out, or if we're down to just one and it could be from a noisy kMer or we've run into the noise level for pairs as well
                    if (viableAlternatives == 0 || (viableAlternatives == 1 && (kMerCounts[chosenIdx] < noiseLevel) || (highestPairDepth < noiseLevel && pli > pairSizes.Count /2)))
                    {
                        for (int b = 0; b < 4; b++)
                            variantAlive[b] = variantPreviouslyAlive[b];
                        viableAlternatives = previousViableAlternatives;
                        break;
                    }

                    if (viableAlternatives == 1)
                    {
                        //if (breakOnRead)
                        //{
                        //    Console.WriteLine("pair decided @" + pairSize + " (" + highestPairDepth + ") [" + chosenIdx + "]");
                        //}
                        pairSizeStats[pli + statsPairsIdx]++;   // stopped at this pair-length with a single good next kMer
                        break;
                    }

                    pli++;
                }

                // don't do a further recursive call if we're already at the maximum call depth
                if (viableAlternatives > 1 && level == maxRecursion)
                {
                    viableAlternatives = 0;
                    if (log != null)
                        log.WriteLine("recursive limit: " + read.ToString());
                    break;
                }

                // still have multiple alternatives, so see how far downstream each of them can get
                if (viableAlternatives > 1)
                {
                    //if (kMer == 0x2B8BA0F58BB22B80)
                    //    Debugger.Break();

                    Sequence[] trialReads = new Sequence[4];
                    string[] trialStartingReads = new string[4];
                    int[] trialLengths = new int[4];
                    HashSet<ulong>[] trialLoopTraps = new HashSet<ulong>[4];
                    bool[] reachedTerminalPrimer = new bool[4];
                    bool[] coinTossed = new bool[4];
                    int tpReachedCount = 0;

                    pairSizeStats[statsDownstreamIdx]++;     // count how many times we needed to look downstream to resolve ambiguity

                    //if (log != null)
                    //    log.WriteLine(new string(' ', level * 4) + level + ": " + "exploring " + read.ToString());

                    // recursively explore each of the potential paths and choose one of them
                    // if multiple paths can reach a terminal primer, choose amongst these alternatives in proportion to their abundance
                    // all of these paths are equally 'correct'... 

                    for (int b = 0; b < 4; b++)
                    {
                        if (variantAlive[b])
                        {
                            //if (log != null)
                            //    log.WriteLine(new string(' ', level * 4) + level + ": " + "append  " + bases[b] + "@" + read.Length);

                            trialReads[b] = new Sequence(read);
                            trialReads[b].Append(bases[b]);
                            trialStartingReads[b] = trialReads[b].ToString();
                            trialLoopTraps[b] = new HashSet<ulong>(loopTrap);

                            if (cachedExtensions.ContainsKey(trialStartingReads[b]))
                            {
                                trialReads[b] = cachedExtensions[trialStartingReads[b]];
                                reachedTerminalPrimer[b] = cachedReachedTerminalPrimer[trialStartingReads[b]];
                                pairSizeStats[statsCachedBranch]++;
                                //if (breakOnRead)
                                //{
                                //    Console.WriteLine("cached: " + "[" + b + "] " + trialStartingReads[b] + " --> " + trialReads[b].ToString());
                                //}
                            }
                            else
                            {
                                trialReads[b] = ExtendRead(trialReads[b], kMerSize, kMerTable, pairSizes, kMerPairs, pairSizeStats, previousDepth, /*averageDepth,*/
                                                           cachedExtensions, cachedReachedTerminalPrimer, trialLoopTraps[b], terminatingPrimers, primerSize, noiseLevel, meanDepth,
                                                           out reachedTerminalPrimer[b], out coinTossed[b], log, level + 1, breakOnRead);
                                //if (breakOnRead)
                                //    Console.WriteLine("extended: " + "[" + b + "] " + trialReads[b].ToString());
                                if (cacheThings && !coinTossed[b])
                                {
                                    cachedExtensions.Add(trialStartingReads[b], trialReads[b]);
                                    cachedReachedTerminalPrimer.Add(trialStartingReads[b], reachedTerminalPrimer[b]);
                                }
                            }
                            trialLengths[b] = trialReads[b].Length;
                            if (reachedTerminalPrimer[b])
                            {
                                tpReachedCount++;
                                chosenIdx = b;
                            }
                        }
                    }

                    // chose amongst the alternatives
                    //      ones that lead to a terminal primer are always chosen if any exist
                    //      if there was only one we would have chosen its index on the way
                    //      if more than one such TP path exists, choose among the paths in proportion to their kMer counts
                    //      if no TP paths were found, just choose the longest path - which will probably be rejected as too short later anyway

                    if (tpReachedCount == 1)
                        pairSizeStats[statsSingleDSTIdx]++;

                    if (tpReachedCount > 1)
                    {
                        int totalTPCounts = 0;
                        int[] tpCounts = new int[tpReachedCount];
                        int[] tpBases = new int[tpReachedCount];
                        int tpcIdx = 0;
                        // only going to choose amongst TP paths
                        for (int b = 0; b < 4; b++)
                            if (reachedTerminalPrimer[b])
                            {
                                totalTPCounts += kMerCounts[b];
                                // cumulative counts
                                tpCounts[tpcIdx] = kMerCounts[b];
                                if (tpcIdx != 0)
                                    tpCounts[tpcIdx] += tpCounts[tpcIdx - 1];
                                // track which path this TP extension came from
                                tpBases[tpcIdx] = b;
                                tpcIdx++;
                            }

                        // chose one of the paths
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
                        pairSizeStats[statsRolledIdx]++;
                        //if (breakOnRead)
                        //    Console.WriteLine("tossed: " + "[" + chosenIdx + "] " + trialReads[chosenIdx].ToString());
                    }

                    if (tpReachedCount == 0)
                    {
                        int longestExtension = 0;
                        int longestIdx = 0;

                        for (int b = 0; b < 4; b++)
                        {
                            if (variantAlive[b] && trialLengths[b] > longestExtension)
                            {
                                longestExtension = trialLengths[b];
                                longestIdx = b;
                            }
                        }

                        chosenIdx = longestIdx;
                        pairSizeStats[statsLongestIdx]++;         // chose longest
                    }

                    //if (kMer == 0x6BB26B80E6C8CDA8 && chosenIdx == 0)
                    //    Debugger.Break();

                    //if (log != null)
                    //{
                    //    log.WriteLine(new string(' ', level * 4) + level + ": " + "explored - chose " + bases[chosenIdx] + " " + longestAlternative);
                    //    if (level == 1)
                    //        log.WriteLine("extended (" + trialReads[chosenIdx].Length + ") " + trialReads[chosenIdx].ToString());
                    //}

                    // and exit here as the chosen TrialRead will be the chosen extended read
                    read = trialReads[chosenIdx];
                    terminatingPrimerReached = reachedTerminalPrimer[chosenIdx];
                    tossedCoin = tossedCoin | coinTossed[chosenIdx]; 
                    break;
                }

                // at least one viable following kMer so extend the read and continue
                //if (log != null && traceRead)
                //    log.WriteLine("added " + bases[chosenIdx] + " " + kMerCounts[0] + "+" + savedMerCounts[0] + " " + kMerCounts[1] + "+" + savedMerCounts[1] + " " +
                //                                                      kMerCounts[2] + "+" + savedMerCounts[2] + " " + kMerCounts[3] + "+" + savedMerCounts[3]);

                previousDepth = kMerCounts[chosenIdx];
                read.Append(bases[chosenIdx]);
                kMer = kMerVariants[chosenIdx];
                kMersInRead++;
                basesAdded++;

                if (loopTrap.Contains(kMer))
                {
                    //if (log != null)
                    //    log.WriteLine(new string(' ', level * 4) + level + ": ER loop (i) with " + kMers.ExpandMer(kMer, kMerSize) + " in " + read.ToString());
                    break;
                }
                else
                    loopTrap.Add(kMer);

                if (kMerIsTerminating(kMer, kMerSize, terminatingPrimers, primerSize))
                {
                    terminatingPrimerReached = true;
                    break;
                }

            } // while (extending)

            //if (read.ToString() == "GTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTGTCCGGAATCATTGGGCGTAAAGAGCACGTAGGCGGCCTGTCAAGTCGGATGTGAAAACCCGTGGCTCAACTGCGGGACTGCATTCGAAACTGGCAGGCTAGAGTCCAGTAGAGGAAAGTGGAATTCCCAGTGTAG")
            //    Debugger.Break();

            return read;
        }

        private static bool kMerIsTerminating(ulong kMer, int kMerSize, HashSet<ulong> terminatingPrimers, int primerSize)
        {
            ulong potentialPrimer = kMer << (kMerSize - primerSize) * 2;
            return terminatingPrimers.Contains(potentialPrimer);
        }

        // Find a kMer in the canonical kMer table.
        private static int kMerTableLookup(ulong kMer, Dictionary<ulong, int> kMerTable, int kMerSize)
        {
            int count = 0;
            kMer = Canonical(kMer, kMerSize);

            if (kMerTable.ContainsKey(kMer))
                count = kMerTable[kMer];

            return count;
        }
    }
}


