#undef TRACE
#undef DEBUG
#undef HLTRACE

using System.Diagnostics;
using WorkingDogsCore;

namespace Clancy
{
    // Compares (at a k-mer level) the set of sequences from a set of genes from one genome onto some other (whole) genomes.
    // 
    // Clancy is an offspring of Pup - but Pup is a mostly-symmetrical GxG comparer, for generating gene/genome identity matrices - and Clancy maps just one file of called genes (.ffn) onto other (usually whole genomes).
    // MapGenomesOntoGenomes does much the same thing without the Blue-inspired match-with-correction goodness.
    //
    // Clancy takes a file of called genes, and any number of target genomes or genes. All files are assumed to be in FASTA format.
    //
    // The -targets option says how to treat the target.
    //      '-target genes' says that the target files contain genes seqs, and each of these seqs will be a target.
    //      '-target genomes' says each target file is a complete genomes and the matches will be done against the entire genome. Default is '-target genomes'
    //
    // The comparisons can be exact or single-base variants (fuzzy), and either just subs or subs/dels/ins.
    //
    // The result is a table with each row being a called gene/contig followed by how many matches it got from each of the genomes
    //
    // usage: Clancy [-k kmerSize] [-exact|-inexact|-subs] [-t #threads] [-s] [-genes|-genomes] -out countsFN sourceFN targetFNs
    //
    // NM-KH31_C3.ffn Priestia_aryabhattai_genomes\*.fna Priestia_megaterium_genomes\*.fna -k 25 -s -subs -t 32 -out Clancy_counts_PaPm_subs.txt
    // NM-KH31_C3.ffn Close_relatives\*.fna -k 25 -s -subs -t 32 -out Clancy_counts_CR_subs.txt
    // NM-KH31_C3.ffn "Priestia plasmids.fna" -genes -k 25 -s -subs -t 32 -out Clancy_counts_plasmid_subs.txt

    public enum VariantType
    {
        none,
        sub,
        ins,
        del,
        N,
        invalid
    }

    public enum VariantsAllowed
    {
        subOnly,
        subInsDel
    }

    class Program
    {

        public const string clancyVersion = "v1.0.0";

        enum VariantsWanted
        {
            lastVariantOnly,
            firstVariantOnly,
            allSubvariants
        }

        enum ScanType
        {
            forwardScan,    // normal left-to-right scan/match of the source sequence
            reverseScan     // recursive right-to-left (RC) scan used to recover matches after a gap
        }

        enum MatchType
        {
            full,           // single base sub, ins & del gaps when looking for next match
            scanning,       // single/multiple subs, or grow seed matches when trying to find the next anchoring match
            spanGap         // try to span a long gap
        }

        enum MatchMode
        {
            continuing,     // continuing a match
            scanning        // scanning for the start of a matching region
        }

        public enum Strand
        {
            plus = 0,
            minus = 1,
            either = 2
        }

        // costs of various changes - feed into the cost function
        public enum Costs
        {
            costMatch = 0,
            costSub = 1,
            costGapOpen = 2,
            costGapExtend = 3
        }

        const int reducedMatch = 0;
        const int reducedX1Sub = 1;
        const int reducedX1InDel = 2;
        const int reducedX2Sub = 1;
        const int reducedX2Indel = 2;
        const int reducedGapIndel = 3;
        const int reducedGapInversion = 0;

        const int plusStrand = 0;
        const int minusStrand = 1;

        static int largeVariantSetSize = 0;
        static int smallVariantSetSize = 0;

        static int maxMismatches = 5;
        static int maxCost = 10;
        static int maxDownstream = 20;
        static int smallGap = 2;
        static int biggerGap = 5;
        static int longGap = 12;
        static int maxConsecutiveFixes = 5;
        static int trailingCostsLength = 10;
        static int scanCheckInterval = 10;
        static bool skipping = true;
        static bool doubleSubs = false;

        static bool useSeeds = true;
        static int seedLength = 0;
        static ulong seedMask = 0;
        static int maxSeedMismatches = 0;
        const int smallestSeedLength = 8;
        const int longestSeedLength = 14;
        const double maxSeedUsage = 0.05;
        const int noViableVariants = -1;
        const int maxDegenerate = 5;

        static TiledGenome[]? targets = null;

#if PERF
        static long[][] perfCounts;
        const int noOfCounts = 25;
        const string perfHeaders = "Exact\tScanTests\tScanSingleMatches\tScanMultMatches\tExt1Tests\tExt1SingleMatches\tExt1MultMatches\tExt2Tests\tExt2Matches\tExt3Tests\tExt3Matches\tMatched\tCostly\tNoMatch\tSeedTest\tSeedViable\tSeedGrown\tScanMTests\tScanMMatches\tScanMMatchedM\tScanBackTracks\tSkipbacks\tRvsExact\tRvsFuzzy\tSkipbackProbs";
        //                          0     1          2                  3                4          5                  6                7          8            9          10           11       12      13       14        15          16         17          18            19             20              21         22        23
        enum PC //  ptGrowSeed
        {
            pcExact, pcScanTests, pcScanSingleMatches, pcScanMultMatches, pcExt1Tests, pcExt1SingleMatches, pcExt1MultMatches, pcExt2Tests, pcExt2Matches, pcExt3Tests, pcExt3Matches, pcMatched, pcCostly, pcNoMatch, pcSeedTest, pcSeedViable, pcSeedGrown, pcScanMTests, pcScanMMatches, pcScanMMatchedM, pcScanBackTracks, pcSkipbacks, pcRvsExact, pcRvsFuzzy, pcSkipbackProbs
        }
        static Stopwatch[][] perfTimers;
        const int noOfPerfTimers = 9;
        const string perfTimersHeaders = "Scan\tTVMScan\tTVM1\tTVM2\tTVM3\tRVSE\tRVSF\tGVMV\tGrowSeed";
        //                                0     1        2     3    4     5     6     7     8
        enum PT
        {
            ptScan, ptTVMScan, ptTVM1, ptTVM2, ptTVM3, ptRVSE, ptRVSF, ptGVMV, ptGrowSeed
        }
#endif

#if (DEBUG)
        static int foundInMM = 0;
        static int foundIn2xMM = 0;
#endif
#if (TRACE || HLTRACE)
        static StreamWriter trace = new StreamWriter("trace.txt");
#endif
#if (PERF)
        static StreamWriter stats;
#endif

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                WriteUsage();
                return;
            }

            string? sourceFlag = null;
            SearchOption so = SearchOption.TopDirectoryOnly;
            bool exactMatch = false;
            VariantsAllowed allowedVariants = VariantsAllowed.subOnly;
            int kMerSize = 21;
            int noThreads = Environment.ProcessorCount/2;
            List<string> FNParams = new List<string>();
            string? countsFN = null;
            // each target file is treated as a single target, or the single target file is treated as a file of target seqs (genes/contigs)
            bool mapToGenome = true;
            // min length for target - useful when matching to metagenome
            int minLength = 0;

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-k")
                    {
                        if (!CheckForParamValue(p, args.Length, "kMer length number expected after -k"))
                            return;
                        try
                        {
                            kMerSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        if (kMerSize > 31)
                        {
                            Console.WriteLine("kMer length must be <= 31");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-s")
                    {
                        so = SearchOption.AllDirectories;
                        continue;
                    }

                    if (args[p] == "-targets")
                    {
                        string targetsArg = args[p+1].ToLower();
                        if (targetsArg == "genes")
                        {
                            mapToGenome = false;
                            p++;
                            continue;
                        }
                        if (targetsArg == "genomes")
                        {
                            mapToGenome = true;
                            p++;
                            continue;
                        }
                        Console.WriteLine("expected genes or genomes after -targets: " + args[p]);
                        return;
                    }

                    if (args[p] == "-noskip")
                    {
                        skipping = false;
                        continue;
                    }

                    if (args[p] == "-subx2")
                    {
                        doubleSubs = true;
                        continue;
                    }

                    if (args[p] == "-subx1")
                    {
                        doubleSubs = false;
                        continue;
                    }

                    if (args[p] == "-skip")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -skip"))
                            return;
                        try
                        {
                            scanCheckInterval = Convert.ToInt32(args[p + 1]);
                            skipping = true;
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -skip parameter: " + args[p + 1]);
                            return;
                        }

                        p++;
                        continue;
                    }

                    if (args[p] == "-t" || args[p] == "-threads")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -t|-threads"))
                            return;
                        try
                        {
                            noThreads = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -t|-threads parameter: " + args[p + 1]);
                            return;
                        }

                        p++;
                        continue;
                    }

                    if (args[p] == "-minlength")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -minlength"))
                            return;
                        try
                        {
                            minLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -minlength parameter: " + args[p + 1]);
                            return;
                        }

                        p++;
                        continue;
                    }

                    if (args[p] == "-mm" || args[p] == "-maxmismatch")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -mm|-maxMismatch"))
                            return;
                        try
                        {
                            maxMismatches = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -mm|-maxMismatch parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-mc" || args[p] == "-maxcost")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -mc|-maxCost"))
                            return;
                        try
                        {
                            maxCost = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -mc|-maxCost parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-mf" || args[p] == "-maxfollowers")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -mf|-maxFollowers"))
                            return;
                        try
                        {
                            maxDownstream = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -mf|-maxFollowers parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-mg" || args[p] == "-maxgap")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -mg|-maxGap"))
                            return;
                        try
                        {
                            biggerGap = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -mg|-maxGap parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-sl" || args[p] == "-seedlength")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -sl|-seedLength"))
                            return;
                        try
                        {
                            seedLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -sl|-seedLength parameter: " + args[p + 1]);
                            return;
                        }
                        seedMask = 0xffffffffffffffff << (64 - seedLength * 2);
                        if (seedLength == 0)
                            useSeeds = false;
                        p++;
                        continue;
                    }

                    if (args[p] == "-sm" || args[p] == "-seedmismatch")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -sm|-seedMismatch"))
                            return;
                        try
                        {
                            maxSeedMismatches = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -sm|-seedMismatch parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-out" || args[p] == "-output")
                    {
                        if (!CheckForParamValue(p, args.Length, "file name expected after -out"))
                            return;

                        countsFN = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-inexact")
                    {
                        exactMatch = false;
                        allowedVariants = VariantsAllowed.subInsDel;
                        continue;
                    }

                    if (args[p] == "-subs")
                    {
                        exactMatch = false;
                        allowedVariants = VariantsAllowed.subOnly;
                        continue;
                    }

                    if (args[p] == "-exact")
                    {
                        exactMatch = true;
                        useSeeds = false;
                        continue;
                    }

                    Console.WriteLine("Unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                FNParams.Add(args[p]);
            }

            smallVariantSetSize = 3 + 1 + 4;
            largeVariantSetSize = 3 * kMerSize + kMerSize + 4 * kMerSize;

            if (seedLength != 0 && maxSeedMismatches == 0)
            {
                Console.WriteLine("Must specify seed mismatches if seed length is specified");
                return;
            }

            if (FNParams.Count < 2)
            {
                Console.WriteLine("need both source and at least one target file)");
                return;
            }

            string sourceFN = FNParams[0];
            if (!File.Exists(sourceFN))
            {
                Console.WriteLine(sourceFN + " not found");
                return;
            }
            string cleanSourceFN = sourceFN;
            int lastSlosh = sourceFN.LastIndexOf(Path.PathSeparator);
            if (lastSlosh != -1)
                cleanSourceFN = sourceFN.Substring(lastSlosh + 1);
            List<string> targetFNs = new List<string>();

            for (int i = 1; i < FNParams.Count; i++)
            {
                string FN = FNParams[i];
                string path;
                string FNP;
                GetPathFN(FN, out path, out FNP);

                SearchOption searchOption = so;
                if (sourceFlag == "-d")
                    searchOption = SearchOption.AllDirectories;
                string[] foundFNs = Directory.GetFiles(path, FNP, searchOption);

                if (foundFNs.Length == 0)
                    Console.WriteLine("No files found matching " + FN);

                foreach (string foundFN in foundFNs)
                    if (!foundFN.Contains(sourceFN))
                        targetFNs.Add(foundFN);
            }

            int noOfTargets = targetFNs.Count;
            if (noOfTargets == 0)
            {
                Console.WriteLine("no genome target files found to compare against");
                return;
            }

            targetFNs.Sort();

            // check that output file name has been specified
            if (countsFN == null)
            {
                Console.WriteLine("no results file specified (-out)");
                return;
            }

            // mapping to genes? open all target files and read their contigs/genes (and keep these as our targets)
            // mapping to genomes? defer this tiling until needed during matching
            List<string> targetHeaders = new List<string>();
            List<string> targetSeqs = new List<string>();
            if (!mapToGenome)
            {
                foreach (string FN in targetFNs)
                {
                    List<string> th = new List<string>(50000);
                    List<string> ts = new List<string>(50000);
                    GetSequences(FN, kMerSize, th, ts, minLength);
                    foreach (string h in th)
                        targetHeaders.Add(h);
                    foreach (string s in ts)
                        targetSeqs.Add(s);
                }
                noOfTargets = targetHeaders.Count;
            }

#if PERF
            stats = new StreamWriter(countsPrefix + "_stats.txt");
            stats.Write("pup (" + pupVersion + ") ");
            foreach (string arg in args)
                stats.Write(" " + arg);
            stats.WriteLine();
            stats.WriteLine("\t" + perfHeaders + "\t" + perfTimersHeaders);
#endif

            DateTime comparisonStart = DateTime.Now;

            char[] FNSeparator = new char[] { Path.DirectorySeparatorChar };

            //Console.WriteLine("SL=" + seedLength + "; SM=" + maxSeedMismatches + "; K=" + kmerSize);

            // get all the source (gene/contig) seqs
            List<string> geneHeaders = new List<string>(50000);
            List<string> geneSeqs = new List<string>(50000);
            GetSequences(sourceFN, kMerSize, geneHeaders, geneSeqs, 0);

            Console.WriteLine("matching " + geneSeqs.Count + " genes from " + cleanSourceFN + " to " + noOfTargets + " target genes/genomes");

            // counts matrix (source gene#, matches to target genome#)
            long[,] matches = new long[geneSeqs.Count, noOfTargets];
            // timings matrix
            long[,] timesForMatches = new long[geneSeqs.Count, noOfTargets];
            // and a nice form of the target file names
            string[] targetNames = new string[noOfTargets];
            // target sizes (in kMers)
            long[] targetSizes = new long[noOfTargets];

            if (mapToGenome)
                for (int t = 0; t < noOfTargets; t++)
                {
                    string targetName = targetFNs[t].Substring(targetFNs[t].LastIndexOf(Path.DirectorySeparatorChar) + 1);
                    if (targetName.Contains('.'))
                        targetName = targetName.Substring(0, targetName.LastIndexOf('.'));
                    targetNames[t] = targetName;
                }
            else
            {
                for (int t = 0; t < noOfTargets; t++)
                {
                    // >CP120616.1 Priestia megaterium strain DSM 1804 plasmid unnamed1, complete sequence
                    // >NC_014103.1 Priestia megaterium DSM 319, complete sequence
                    string th = targetHeaders[t].Substring(1);
                    int commaIdx = th.IndexOf(",");
                    if (commaIdx > 0)
                        th = th.Substring(0, commaIdx);
                    targetNames[t] = th;
                }
            }
#if PERF
            perfCounts = new long[noThreads][];
            perfTimers = new Stopwatch[noThreads][];
            for (int i = 0; i < noThreads; i++)
            {
                perfCounts[i] = new long[noOfCounts];
                perfTimers[i] = new Stopwatch[noOfPerfTimers];
                for (int t = 0; t < noOfPerfTimers; t++)
                    perfTimers[i][t] = new Stopwatch();
            }
#endif
            // calculate required seed size (if not explicitly specified as a parameter)
            if (useSeeds && seedLength == 0)
            {
                long longestGenome = 0;
                for (int g = 0; g < noOfTargets; g++)
                {
                    long genomeFileSize = targetSeqs.Count == 0 ? new FileInfo(targetFNs[g]).Length : targetSeqs[g].Length;

                    if (genomeFileSize > longestGenome)
                        longestGenome = genomeFileSize;
                }
                double seedCountGuess = longestGenome * 2;
                double seedUsage = 0.0;
                int neededSeedLength = 0;
                for (int sk = 8; sk <= (kMerSize * 3) / 4; sk++)
                {
                    seedUsage = seedCountGuess / Math.Pow(4.0, (double)sk);
                    //Console.WriteLine("seed=" + sk + " usage=" + (seedUsage * 100.0).ToString("F1") + "%");

                    if (Math.Pow(4.0, (double)sk) / 10.0 > seedCountGuess)
                    {
                        neededSeedLength = sk;
                        seedUsage = seedCountGuess / Math.Pow(4.0, (double)sk);
                        break;
                    }
                }
                seedLength = neededSeedLength;
                seedMask = 0xffffffffffffffff << (64 - seedLength * 2);
                maxSeedMismatches = Math.Min((kMerSize - seedLength) / 2, 4);
                if (seedLength == 0)
                {
                    Console.WriteLine("Not using seed kMers: not sufficiently distinct");
                    useSeeds = false;
                }
                else
                    Console.WriteLine("Seed length " + seedLength + " " + (seedUsage * 100.0).ToString("F1") + "%");
            }

            // make sure results file can be opened - better than waiting until the end of a long comparison run
            StreamWriter counts;
            try
            {
                counts = new StreamWriter(countsFN);
            }
            catch
            {
                Console.WriteLine("Failed to open counts file: " + countsFN);
                return;
            }
            counts.Close();

            // cache of previously tiled target genomes (empty until needed, removed after all genomes have been matched against it)
            targets = new TiledGenome[noOfTargets];
            for (int t = 0; t < noOfTargets; t++)
                targets[t] = new TiledGenome(geneSeqs.Count);

            // queue matches all of our source contigs/genes against all of the target genome hash sets
            //   in target order - S1T1, S2T1, S1T2, S2T2
            Queue<MatchingRequest> queuedRequests = new Queue<MatchingRequest>();
            for (int t = 0; t < noOfTargets; t++)
            {
                for (int s = 0; s < geneSeqs.Count; s++)
                {
                    MatchingRequest mr = new MatchingRequest(s, t);
                    queuedRequests.Enqueue(mr);
                }
                MatchingRequest mre = new MatchingRequest(-1, t);
                queuedRequests.Enqueue(mre);
            }

#if (DEBUG || TRACE || HLTRACE)
            noThreads = 1;
#endif
#if (TRACE || HLTRACE)                        
            trace.AutoFlush = true;
#endif
            MatchingThreadParams[] matchingParams = new MatchingThreadParams[noThreads];
            Thread[] matchingThreads = new Thread[noThreads];

            // ready a new thread for each parallel healer
            for (int b = 0; b < noThreads; b++)
            {
                matchingParams[b] = new MatchingThreadParams();
                matchingParams[b].threadNo = b;
                matchingParams[b].exactMatch = exactMatch;
                matchingParams[b].allowedVariants = allowedVariants;
                matchingParams[b].kmerSize = kMerSize;
                matchingParams[b].queuedRequests = queuedRequests;
                matchingParams[b].sourceSeqs = geneSeqs;
                matchingParams[b].targetFNs = targetFNs;
                matchingParams[b].targetSeqs = targetSeqs;
                matchingParams[b].targetNames = targetNames;
                matchingParams[b].targetSizes = targetSizes;
                matchingParams[b].minLength = minLength;
                matchingParams[b].matches = matches;
                matchingParams[b].timesForMatches = timesForMatches;

                matchingThreads[b] = new Thread(new ParameterizedThreadStart(Program.MatchingThread));
                matchingThreads[b].Priority = ThreadPriority.BelowNormal;
                matchingThreads[b].Name = "MatchingWorker#" + b.ToString();
                matchingThreads[b].Start(matchingParams[b]);
                Thread.Sleep(100);
            }

            // and wait for all threads to finish
            for (int b = 0; b < noThreads; b++)
            {
                matchingThreads[b].Join();
                matchingThreads[b] = null;
            }

            // and finally write out the results (as counts, symmetric similarity PIM, sparse similarity PIM

            WriteResults(countsFN, geneHeaders, geneSeqs, kMerSize, targetNames, targetSizes, matches);
#if PERF
            WriteSquare(countsPrefix + "_timings_" + countsSuffix, noOfTargets, targetNames, genomeSizes, timesForMatches, null, 0.0f);
            stats.Close();
#endif

#if (TRACE || HLTRACE)
            trace.Close();
#endif

            Console.WriteLine("Finished cross-comparing " + noOfTargets + " genes/genomes in " + (DateTime.Now - comparisonStart).TotalSeconds.ToString("#.0") + "s");
            if (!useSeeds)
                Console.WriteLine((skipping ? "skip " + scanCheckInterval : "no skip") + ", " + (doubleSubs ? "double" : "single") + " subs");
#if DEBUG
            Console.WriteLine("Found in 1xMM: " + foundInMM + " found in 2xMM: " + foundIn2xMM);
#endif
        }

        private static void WriteResults(string resultsFN, List<string> headers, List<string> seqs, int kMerSize, string[] targetNames, long[] targetSizes, long[,] matches)
        {
            StreamWriter results = new StreamWriter(resultsFN);
            int noOfTargets = targetNames.Length;

            results.Write("\t\t");
            for (int i = 0; i < noOfTargets; i++)
            {
                results.Write(targetNames[i]);
                results.Write("\t");
            }
            results.WriteLine();
            results.Write("\t\t");
            for (int i = 0; i < noOfTargets; i++)
            {
                results.Write(targetSizes[i]);
                results.Write("\t");
            }
            results.WriteLine();

            for (int s = 0; s < headers.Count; s++)
            {
                results.Write(headers[s]);
                results.Write("\t");
                results.Write(seqs[s].Length - kMerSize +1);
                results.Write("\t");

                for (int t = 0; t < noOfTargets; t++)
                {
                    results.Write(matches[s, t]);
                    results.Write("\t");
                }
                results.WriteLine();
            }

            results.Close();
        }

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Clancy [-k kmerSize] [-exact|-inexact|-subs] [-t #threads] [-targets genes|genomes] -out countsFN [-s] sourceFN targetFNs " + clancyVersion);
        }

        private static bool CheckForParamValue(int p, int argsLength, string msg)
        {
            if (p == argsLength)
            {
                Console.WriteLine(msg);
                return false;
            }
            return true;
        }


        private static void MatchingThread(object threadParams)
        {
            MatchingThreadParams theseParams = (MatchingThreadParams)threadParams;
            int threadNo = theseParams.threadNo;
            bool exactMatch = theseParams.exactMatch;
            VariantsAllowed allowedVariants = theseParams.allowedVariants;
            int kmerSize = theseParams.kmerSize;
            Queue<MatchingRequest> queuedRequests = theseParams.queuedRequests;
            List<string> sourceSeqs = theseParams.sourceSeqs;
            List<string> targetFNs = theseParams.targetFNs;
            List<string> targetSeqs = theseParams.targetSeqs;
            string[] targetNames = theseParams.targetNames;
            long[] targetSizes = theseParams.targetSizes;
            int minLength = theseParams.minLength;
            long[,] matches = theseParams.matches;
            long[,] timesForMatches = theseParams.timesForMatches;

            bool finishThread = false;
            MatchingRequest matchingRequest = null;

            List<VariantMerSet> kmerSetPool = new List<VariantMerSet>();
            VariantMerSet outerMerSet = new VariantMerSet(largeVariantSetSize);
            List<GrowingMerSet> growingMerPool = new List<GrowingMerSet>();

            //Console.WriteLine("starting thread " + threadNo);

            while (!finishThread)
            {
                lock (queuedRequests)
                {
                    if (queuedRequests.Count > 0)
                        matchingRequest = queuedRequests.Dequeue();
                    else
                    {
                        finishThread = true;
                        matchingRequest = null;
                    }
                }

                if (matchingRequest != null)
                {
                    int s = matchingRequest.source;
                    int t = matchingRequest.target;

                    if (s == -1)
                    {
                        continue;
                    }

                    //Console.WriteLine("starting match of " + targetNames[s] + " to " + targetNames[t] + " on #" + threadNo);
                    //trace.WriteLine("starting match of " + targetNames[s] + " to " + targetNames[t] + " on #" + threadNo);

                    string? targetSeq = null;
                    string? targetFN = null;
                    if (targetSeqs.Count > 0)
                        targetSeq = targetSeqs[t];
                    else
                        targetFN = targetFNs[t];

                    targetSizes[t] = PrepTargetMerSet(t, targetFN, targetSeq, kmerSize, targetNames[t], minLength);
                    string sourceSeq = sourceSeqs[s];

                    DateTime matchingStart = DateTime.Now;
                    int timeTaken;

                    matches[s, t] = MatchSourceToTarget(exactMatch, allowedVariants, kmerSize, sourceSeqs[s], targets[t], kmerSetPool, growingMerPool, outerMerSet, threadNo, out timeTaken);
                    timesForMatches[s, t] = timeTaken;
                    //Console.WriteLine("matched " + s + ", " + t);

#if PERF
                    lock (stats)
                    {
                        stats.Write(s + "," + t + "\t");
                        for (int c = 0; c < noOfCounts; c++)
                        {
                            stats.Write(perfCounts[threadNo][c] + "\t");
                            perfCounts[threadNo][c] = 0;
                        }
                        for (int c = 0; c < noOfPerfTimers; c++)
                        {
                            stats.Write(perfTimers[threadNo][c].Elapsed.TotalSeconds.ToString("F1") + "\t");
                            perfTimers[threadNo][c].Reset();
                        }
                        stats.WriteLine();
                    }
#endif
                    ReleaseTargetMers(t);
                }
            }

            //Console.WriteLine("finished matching thread " + threadNo);
        }

        private static long PrepTargetMerSet(int targetNo, string targetFN, string? targetSeq, int kmerSize, string targetName, int minLength)
        {
            // first matcher to want this target genome? then load it up its hashtables
            if (targets[targetNo].mersFromGenomePlus == null)
            {
                lock (targets)
                {
                    // some other thread may have loaded this table while we were waiting
                    if (targets[targetNo].mersFromGenomePlus != null)
                        return targets[targetNo].mersFromGenomePlus.Count;

                    //if (targetNo == 10)
                    //    Debugger.Break();

                    //Console.WriteLine("loading targets[" + targetNo + "]");
                    long totalLength = 0;
                    List<string> headers = new List<string>(10000);
                    List<string> genomeSequences = new List<string>(10000);
                    bool readTarget = false;

                    if (targetSeq == null)
                    {
                        GetSequences(targetFN, kmerSize, headers, genomeSequences, minLength);
                        foreach (string genomeSequence in genomeSequences)
                            totalLength += genomeSequence.Length;
                        readTarget = true;
                    }
                    else
                    {
                        genomeSequences.Add(targetSeq);
                        totalLength = targetSeq.Length;
                    }

                    MerHashSet plusMers = new MerHashSet(totalLength, kmerSize, false);
                    MerHashSet minusMers;
                    MerHashSet allMers;
                    MerHashSet targetSeeds = null;

                    Stopwatch tilingSW = new Stopwatch();
                    tilingSW.Start();

                    // generate the distinct k-mers for the next target genome (not canonical) and build the 'plus' kMer set
                    Stopwatch plusSW = new Stopwatch();
                    plusSW.Start();
                    long mersInTarget = GetMersFromSequences(genomeSequences, kmerSize, plusMers);
                    plusSW.Stop();

                    if (readTarget)
                        Console.WriteLine("read & tiled target (" + targetName + ") in " + plusSW.Elapsed.TotalSeconds.ToString("F1") + "s" + " " + mersInTarget + " kMers");

                    // build the set of the RC forms of the distinct 'plus' kMers, and the plus+minus set used for scanning
                    minusMers = new MerHashSet(plusMers.Count, kmerSize, false);
                    allMers = new MerHashSet(plusMers.Count * 2, kmerSize, false);
                    if (seedLength != 0)
                    {
                        long targetSeedsLength = Math.Min((long)(plusMers.Count * 2), (long)Math.Pow(4, seedLength));
                        targetSeeds = new MerHashSet(targetSeedsLength, seedLength, false);
                    }

                    Stopwatch restSW = new Stopwatch();
                    restSW.Start();
                    foreach (ulong mer in plusMers)
                    {
                        ulong merRC = kMers.ReverseComplement(mer, kmerSize);
                        minusMers.AddNoCheck(merRC);
                        allMers.AddNoCheck(mer);
                        allMers.AddNoCheck(merRC);
                        if (useSeeds)
                        {
                            targetSeeds.AddIfNotPresent(mer & seedMask);
                            targetSeeds.AddIfNotPresent(merRC & seedMask);
                        }
                    }
                    //minusMers.Optimise();
                    //allMers.Optimise();
                    //Dictionary<int, int>[] plusChecks = plusMers.CheckHashTable();
                    //Dictionary<int, int>[] allChecks = allMers.CheckHashTable();
                    //Dictionary<int, int>[] rcChecks = minusMers.CheckHashTable();
                    //if (seedLength != 0)
                    //{
                    //    //Dictionary<int, int>[]  seedChecks = targetSeeds.CheckHashTable();
                    //    targetSeedsUsed = (double)targetSeeds.Count / Math.Pow(4, seedLength);
                    //}

                    restSW.Stop();
                    //Console.WriteLine("rc/all kMers in " + restSW.Elapsed.TotalSeconds.ToString("F1") + "s");

                    tilingSW.Stop();

                    //if (readTarget)
                    //{
                    //    string seedMsg = useSeeds ? (targetSeeds.Count + " " + seedLength + "-mer seeds") : "no seeds";
                    //    Console.WriteLine("finished loading target #" + targetNo + " (" + targetName + ")" + " plus=" + plusMers.Count + " both=" + allMers.Count + " (" + seedMsg + ")" + " in " + tilingSW.Elapsed.TotalSeconds.ToString("F1") + "s");
                    //}
                    targets[targetNo].seedsFromGenome = targetSeeds;
                    targets[targetNo].mersFromGenomeMinus = minusMers;
                    targets[targetNo].mersFromGenomeEither = allMers;
                    targets[targetNo].mersFromGenomePlus = plusMers;        // last act before lock is released
                }
            }

            return targets[targetNo].mersFromGenomePlus.Count;
        }

        private static void ReleaseTargetMers(int targetNo)
        {
            lock (targets)
            {
                targets[targetNo].refCount--;
                if (targets[targetNo].refCount == 0)
                {
                    targets[targetNo].mersFromGenomePlus = null;
                    targets[targetNo].mersFromGenomeMinus = null;
                    targets[targetNo].mersFromGenomeEither = null;
                    targets[targetNo].seedsFromGenome = null;
                    //GC.Collect(2);
                    //Console.WriteLine("released target set " + targetNo);
                }
            }
        }

        private static long MatchSourceToTarget(bool exactMatch, VariantsAllowed allowedVariants, int kmerSize, string sourceSeq, TiledGenome target, 
                                               List<VariantMerSet> kmerSetPool, List<GrowingMerSet> growingMerPool, VariantMerSet outerMerSet, int threadNo, out int timeTaken)
        {
            long matchingMers = 0;
            timeTaken = 0;
            Stopwatch sw = new Stopwatch();
            sw.Start();

            //Console.WriteLine("matching " + sourceName + " to " + targetName);

            // tile and match this sequence against the current target
            Sequence sequence = new Sequence(sourceSeq);
            if (sequence.Length < kmerSize)
                return 0;

                int mersMatchedForSeq = MatchSequenceToTarget(exactMatch, allowedVariants, ScanType.forwardScan, kmerSize, sequence, 0, target, Strand.either, kmerSetPool, outerMerSet, growingMerPool, threadNo);
                matchingMers += mersMatchedForSeq;

            sw.Stop();
            timeTaken = (int)sw.ElapsedMilliseconds;
            return matchingMers;
        }

        private static int MatchSequenceToTarget(bool exactMatchesOnly, VariantsAllowed allowedVariants, ScanType scanDirection, int kmerSize, Sequence sequence, int start, TiledGenome target, Strand strand,
                                                 List<VariantMerSet> kmerSetPool, VariantMerSet variantMerSet, List<GrowingMerSet> growingMerPool, int threadNo)
        {
            int matches = 0;
            ulong kMer = 0;                 // valid/matched kmer 
            int kMerNo = 0;                 // kmer number for tracing etc
            int lastMatchingMerNo = 0;      // last time we had a match of any form (used to control all/last variants)
            int lastMatchingFBI = 0;        // FBI value for the last match (E or V). This is the index of the first unmatched base of a possible unmatched run. 
            bool inMismatchRegion = false;  // in a mismatch region (1 or more consecutive 'no match' k-mers)
            MatchMode matchingMode = MatchMode.scanning;
            int consecutiveNoMatches = 0;   // number of consecutive non-matching kMers (gapped scans)
            int consecutiveExactMatches = 0;// and the number of consecutive matches
            int consecutiveFuzzyMatches = 0;// and the number of consecutive variant matches
            List<int> trailingCosts = new List<int>(trailingCostsLength);     // costs for the last <=10 matches (exact or inexact)      
            List<ulong> trailingMers = new List<ulong>(trailingCostsLength);  // and the last (consecutive) exact-matching kMers (used to check for should-have-been changed kMers on the first mismatch after a run of matches)
            int cumulativeCost = 0;         // sum of trailing costs
            int skipbackFBI = -1;           // last point where we skipped back to check for a variant to an exact matching kMer (stops loops)
            int lastScanbackFBI = 0;
            HashSet<int> skipbackTrap = new HashSet<int>();
            VariantType previousVariantType = VariantType.none;
            int previousIndelMerNo = -1;
            int indelBalance = 0;

            MerHashSet currentTargetMers = null;    // current set of target kMers (plus, minus or all)
            Strand currentStrand = strand;          // current preferred target kMer set (either, plus or minus)
            if (strand == Strand.either)
                currentTargetMers = target.mersFromGenomeEither;
            if (strand == Strand.plus)
                currentTargetMers = target.mersFromGenomePlus;
            if (strand == Strand.minus)
                currentTargetMers = target.mersFromGenomeMinus;

            // get the starting k-mer - skipping over any initial k-mers containing N or other non-ACGT bases - and we'll assume there's at least one valid k-mer in the sequence
            bool gotGoodMer = false;
            while (!gotGoodMer)
            {
                gotGoodMer = Sequence.CondenseMer(sequence, start, kmerSize, out kMer);
                if (gotGoodMer)
                    break;
                start++;
                kMerNo++;
                // ran off the end before finding a valid k-mer
                if (start >= sequence.Length - kmerSize)
                    return 0;
            }
            // adjust first kmer so following CondenseMerIncremental will work 
            kMer = kMer >> 2;
            int fbi = start + kmerSize - 1;              // initial followingBaseIndex is last base in first k-mer
            int kMerNoTrap = sequence.Length + Math.Max(sequence.Length / 10, 1000);

            // skim along the source sequence 

            // tile and compare all of the k-mers in this sequence (from the beginning if forward pass, and from 'start' for reverse pass)
            // ---------------------- main scanning loop ---------------------------
            while (fbi < sequence.Length)
            {
                ulong previousMer = kMer;
                bool merValid = true;

                if (kMerNo > kMerNoTrap)
                {
                    // never expected to happen - here to trap del loops and other such bugs
                    Console.WriteLine("kmerNo off end of sequence:" + "kmerNo=" + kMerNo + " " + scanDirection);
                    //Debugger.Break();
                    return matches;
                }

                // generate the next kMer
                merValid = kMers.CondenseMerIncremental(kmerSize, kMer, sequence.Bases[fbi], out kMer);
                fbi++;

                // next base is invalid so replace it with a nonviable one if possible, or skip until we have another valid kmer
                if (!merValid)
                {
                    // count consecutive Ns (until the next non-N base)
                    int consecutiveNs = CountConsecutiveNs(sequence, fbi - 1);
                    // if the number of Ns is small enough, replace each N with a non-viable base (so it will be corrected)
                    if (consecutiveNs < maxDegenerate)
                    {
                        char chosenBase;
                        kMer = GenerateNonViableNVariant(previousMer, kmerSize, currentTargetMers, out chosenBase);
                        sequence.Bases[fbi - 1] = chosenBase;
                    }
                    else
                    // too many consecutive mismatches, so just skip until we get a good kMer
                    {
                        while (!merValid)
                        {
                            merValid = Sequence.CondenseMer(sequence, (fbi - kmerSize + 1), kmerSize, out kMer);
                            fbi++;
                            kMerNo++;
                            if (fbi == sequence.Length)
                                return matches;
                        }
                    }
                }

                //if (kmer == 0x5B4737A586000000)
                //    Debugger.Break();

                // first try for an exact match
                // ----------------------------
                if (currentTargetMers.Contains(kMer))
                {
                    matches++;
#if PERF
                    perfCounts[threadNo][(int)PC.pcExact]++;
#endif
                    lastMatchingMerNo = kMerNo;
                    matchingMode = MatchMode.continuing;
                    consecutiveNoMatches = 0;
                    consecutiveFuzzyMatches = 0;
                    consecutiveExactMatches++;
                    int droppedCost = 0;
                    previousVariantType = VariantType.none;

                    if (trailingCosts.Count == trailingCostsLength)
                    {
                        droppedCost = trailingCosts[0];
                        trailingCosts.RemoveAt(0);
                    }
                    cumulativeCost = cumulativeCost - droppedCost + (int)Costs.costMatch;
                    trailingCosts.Add((int)Costs.costMatch);

                    if (trailingMers.Count == trailingCostsLength)
                        trailingMers.RemoveAt(0);
                    trailingMers.Add(kMer);
#if ( TRACE || HLTRACE )
                    //if (kmer == 0xD07F904816000000)
                    //    Debugger.Break();
                    trace.WriteLine(kMerNo + ": E\t" + kMers.ExpandMer(kMer, kmerSize) + " " + fbi +
                        " +" + (fbi == sequence.Length ? "end" : sequence.Bases[fbi].ToString()) + " " + currentStrand.ToString() + " " +
                        (scanDirection == ScanType.reverseScan ? (" [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                        " cc=" + cumulativeCost + " cm=" + consecutiveExactMatches);
#endif
                    // set current strand for the matching kMer if we're scanning, rather than continuing
                    if (currentStrand == Strand.either)
                        SetCurrentStrand(kMer, target, out currentStrand, out currentTargetMers);

                    kMerNo++;       // deferred until after trace

                    // found an exact match while scanning so it's possible that inexact matches preceding this kMer could also have succeeded (but were too expensive to try)
                    // so we'll scan/match backwards from here, allowing mismatches and continue only until the first mismatch is encountered
                    if (inMismatchRegion && scanDirection == ScanType.forwardScan && !exactMatchesOnly)
                    {
#if ( TRACE || HLTRACE )
                        trace.WriteLine(kMerNo + ": rescan after exact match");
#endif
#if PERF
                        perfCounts[threadNo][(int)PC.pcRvsExact]++;
                        perfTimers[threadNo][(int)PT.ptRVSE].Start();
#endif
                        matches += RescanMismatchRegion(exactMatchesOnly, allowedVariants, kmerSize, sequence, target, currentStrand, kMer, lastMatchingFBI, fbi, kmerSetPool, growingMerPool, variantMerSet, threadNo);
#if PERF
                        perfTimers[threadNo][(int)PT.ptRVSE].Stop();
#endif
                    }

                    inMismatchRegion = false;
                    lastMatchingFBI = fbi;
                    continue;
                } // exact match

                // didn't get an exact match so try for a match to a variant k-mer. 
                // Any matches here won't be counted to the total but the source sequence will be 'corrected' (in principle only) to avoid mismatch craters.
                if (!exactMatchesOnly)
                {
                    bool foundReplacement = false;                              // was a replacement found at all?
                    bool foundStrongReplacement = false;                        // was a good replacement found? (otherwise we'll try a bit harder in some circumstances)
                    Mer bestReplacement = null;                                 // best replacement kMer found (and its cost)
                    int replacementIdx = -1;                                    // and its index in the variant list
                    int lengthLeft = Math.Min(sequence.Length - fbi, maxDownstream);    // how many more bases downstream to we have to look for matches
                    int costOfVariant = 0;                                      // cost of a replacement kMer
                    bool checkForVariant = true;                                // check for a variant during scanning (every N bases is checked)
                    bool retryCurrentMer = false;                               // should the current kmer be included in the list of variants?

                    // if we have hit a mismatch after a run of matches, it's possible that we missed an earlier 'correction' when
                    // a kMer from our source matched the 'wrong' one from the target (good kMer but different context)
                    // we'll try going backwards for a little while to see if we come across such an ambiguous kMer and 
                    // try variants of this kmer in case they give better results than the exact match.
                    if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && consecutiveExactMatches >= trailingCostsLength && fbi > skipbackFBI)
                    {
#if (TRACE)
                        trace.WriteLine(kMerNo + ": looking backwards from " + kMerNo + " on mismatch " + kMers.ExpandMer(kMer, kmerSize) + " after " + consecutiveExactMatches + " exact matches");
                        //for (int i = 0; i < trailingMers.Count; i++)
                        //    trace.WriteLine(i + " " + kMers.ExpandMer(trailingMers[i], kmerSize));
                        //if (fbi == 3181508)
                        //    Debugger.Break();
#endif
                        // remember we've done a skipback from here (stops loops)
                        skipbackFBI = fbi;

                        int skippedExactMatches = 0;
                        for (int i = trailingMers.Count - 1; i >= 0; i--)
                        {
                            // trailingMers holds the last 'few' exact matching kMers. These must be consecutive as trailingMers is cleared on a mismatch.tra
                            ulong trailingMer = trailingMers[i];
                            variantMerSet.Reset();
                            variantMerSet.AddViable(VariantType.none, trailingMer, 0, 0, (int)Costs.costMatch);
                            skippedExactMatches++;

                            // generate all viable alternatives to this exact-matching kMer (including the kMer itself)
                            if (CheckForViableVariants(trailingMer, currentTargetMers, kmerSize, variantMerSet))
                            {
                                // found an exact match with alternatives, so generate the metrics on all of them and choose the best
                                int bestAlternativeIdx = TryVariantMers(sequence, fbi - skippedExactMatches, kmerSize, MatchType.full, allowedVariants, currentTargetMers, target.seedsFromGenome, lengthLeft + skippedExactMatches, 0, 1, kMerNo - skippedExactMatches, maxMismatches,
                                                         smallGap, variantMerSet, kmerSetPool, growingMerPool, previousVariantType, previousIndelMerNo, indelBalance, threadNo);

                                // if one of the alternatives is better than the initial match ([0]), switch over to this new kMer and fix up the scan state to match
                                if (bestAlternativeIdx > 0 && variantMerSet.merSet[bestAlternativeIdx].exactMatchesFollowing > variantMerSet.merSet[0].exactMatchesFollowing)
                                {
#if (TRACE || HLTRACE)
                                    trace.WriteLine(kMerNo + ": found alternative to exact match. " + kMers.ExpandMer(trailingMer, kmerSize) + "-->" + kMers.ExpandMer(variantMerSet.merSet[bestAlternativeIdx].mer, kmerSize) + " bem=" + variantMerSet.merSet[bestAlternativeIdx].exactMatchesFollowing); ;
#endif
                                    for (int t = i; t < trailingMers.Count; t++)
                                    {
                                        trailingMers.RemoveAt(t);
                                        cumulativeCost -= trailingCosts[t];
                                        trailingCosts.RemoveAt(t);
                                    }

                                    fbi = fbi - skippedExactMatches;
                                    matches = matches - skippedExactMatches;
                                    consecutiveExactMatches = consecutiveExactMatches - skippedExactMatches;
                                    kMerNo -= skippedExactMatches;
                                    kMer = trailingMer;
                                    lengthLeft += skippedExactMatches;
                                    previousVariantType = VariantType.none;
                                    previousIndelMerNo = -1;
                                    indelBalance = 0;
                                    retryCurrentMer = true;
#if (TRACE || HLTRACE)
                                    trace.WriteLine(kMerNo + ": skipping back over " + skippedExactMatches + " matched variants at " + kMerNo + " " + kMers.ExpandMer(trailingMer, kmerSize));
#endif
#if PERF
                                    perfCounts[threadNo][(int)PC.pcSkipbacks]++;
#endif
                                }
                                // stop after first exact kmer with viable alternatives
                                break;
                            }
                        }
                        // fall through to look for variants of the ambiguous kMer (matchingMode = continuing still)
#if (TRACE || HLTRACE)
                        if (!retryCurrentMer)
                            trace.WriteLine(kMerNo + ": found no alternatives to exact matches");
#endif
                    }

                    // looking for an alternative on a mismatch, so stop counting matches
                    if (!retryCurrentMer)
                        consecutiveExactMatches = 0;

                    // non-seed scanning only checks the occasional mismatching kMer for matching variants (and then backtracks to check intervening kMers if a match is found)
                    if (skipping && matchingMode == MatchMode.scanning && !useSeeds && consecutiveNoMatches > 0)
                        checkForVariant = consecutiveNoMatches % scanCheckInterval == 0;
                    else
                        checkForVariant = true;

                    //if (kMer == 0xEA5A124E65000000)
                    //    Debugger.Break();
                    //if (kMerNo == 1532)
                    //    Debugger.Break();
#if (DEBUG || TRACE || HLTRACE)
                    string whenMatched = "";
#endif
                    // look for the best variant of the current (usually non-matching) kMer (if one exists)
                    if (checkForVariant)
                    {
                        foundReplacement = false;
                        foundStrongReplacement = false;
                        bestReplacement = null;
                        MatchType matchType = matchingMode == MatchMode.scanning ? MatchType.scanning : MatchType.full;
#if PERF
                        perfCounts[threadNo][matchingMode == MatchMode.scanning ? (int)PC.pcScanTests : (int)PC.pcExt1Tests]++;
                        perfTimers[threadNo][(int)PT.ptGVMV].Start();
#endif
                        // generate an initial set of variant kMers (variant type depends on whether we're scanning or extending) 
                        int fbiAdjustment = 0;
                        ulong overlookedMer;
                        variantMerSet.Reset();
                        bool checkFor2xSubs = doubleSubs ? (consecutiveNoMatches % scanCheckInterval == 0) : false;

                        int noOfVariantMers = GenerateViableMerVariants(matchType, allowedVariants, checkFor2xSubs, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, growingMerPool, sequence, fbi, lastMatchingFBI, (kMerNo - previousIndelMerNo),
                                                                        currentTargetMers, target.seedsFromGenome, smallGap, lastScanbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out lastScanbackFBI);
#if PERF
                        perfTimers[threadNo][(int)PT.ptGVMV].Stop();
#endif
                        // backtrack after a successful (but occasional scanning match) and fall into the checking-variants code
                        if (fbiAdjustment > 0)
                        {
#if TRACE || HLTRACE
                            trace.WriteLine(kMerNo + ": skipping back to " + (kMerNo - fbiAdjustment) + " (" + kMers.ExpandMer(overlookedMer, kmerSize) + ")");
#endif
                            consecutiveNoMatches -= fbiAdjustment;
                            fbi -= fbiAdjustment;
                            kMerNo -= fbiAdjustment;
                        }

                        if (noOfVariantMers == 1)
                        {
                            foundReplacement = true;
                            foundStrongReplacement = true;
                            bestReplacement = variantMerSet.merSet[0];
                            costOfVariant = (int)Costs.costSub;
#if PERF
                            perfCounts[threadNo][(int)PC.pcScanSingleMatches]++;
#endif
#if TRACE
                            trace.WriteLine(kMerNo + ": " + kMers.ExpandMer(kMer, kmerSize) + " (1)--> " + kMers.ExpandMer(bestReplacement.mer, kmerSize));
#endif
                        }

                        // if there's more than one viable variant, choose between them by considering their followers
                        if (noOfVariantMers > 1)
                        {
#if PERF
                            perfTimers[threadNo][matchingMode == MatchMode.scanning ? (int)PT.ptTVMScan : (int)PT.ptTVM1].Start();
#endif
#if TRACE
                            string variants = "";
                            for (int i = 0; i < variantMerSet.setSize; i++)
                                variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                            trace.WriteLine(kMerNo + ": TVM1 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif
                            // first try for an fairly close match (small number of mismatches, small gaps)
                            replacementIdx = TryVariantMers(sequence, fbi, kmerSize, (matchingMode == MatchMode.scanning ? MatchType.scanning : MatchType.full),
                                                           allowedVariants, currentTargetMers, target.seedsFromGenome, lengthLeft, 0, 1, kMerNo, maxMismatches,
                                                           smallGap, variantMerSet, kmerSetPool, growingMerPool, previousVariantType, previousIndelMerNo, indelBalance, threadNo);
#if PERF
                            perfTimers[threadNo][matchingMode == MatchMode.scanning ? (int)PT.ptTVMScan : (int)PT.ptTVM1].Stop();
#endif

                            // was a replacement found at all?
                            if (replacementIdx >= 0)
                            {
#if PERF
                                perfCounts[threadNo][matchingMode == MatchMode.scanning ? (int)PC.pcScanMultMatches : (int)PC.pcExt1MultMatches]++;
#endif
                                bestReplacement = variantMerSet.merSet[replacementIdx];
                                //if (bestReplacement.variantType != VariantType.sub)
                                //    Debugger.Break();
                                foundReplacement = true;
                                foundStrongReplacement = (bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing) >= lengthLeft && bestReplacement.costOfFollowing < 0;

                                costOfVariant = bestReplacement.costOfVariant;
#if DEBUG
                                if (foundStrongReplacement)
                                {
                                    foundInMM++;
                                    whenMatched = "x1";
                                }
#endif
                            } // first attempt at choosing between viable variants

                            // if we didn't find a strong replacement - and we're continuing a forward matching region - we'll try a more expensive exploration with more mismatches allowed
                            if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && foundReplacement && !foundStrongReplacement && trailingMers.Count > (trailingCostsLength / 2) && cumulativeCost < 1)
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcExt2Tests]++;
                                perfTimers[threadNo][(int)PT.ptTVM2].Start();
#endif
                                // generate the variants again, with 'full' and a larger allowable gap
                                variantMerSet.Reset();
                                noOfVariantMers = GenerateViableMerVariants(MatchType.full, allowedVariants, false, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, growingMerPool, sequence, fbi, lastMatchingFBI, (kMerNo - previousIndelMerNo),
                                                                                currentTargetMers, target.seedsFromGenome, biggerGap, skipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out skipbackFBI);
#if TRACE
                                variants = "";
                                for (int i = 0; i < variantMerSet.setSize; i++)
                                    variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                                trace.WriteLine(kMerNo + ": TVM2 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif

                                replacementIdx = TryVariantMers(sequence, fbi, kmerSize, MatchType.full, allowedVariants,
                                                                currentTargetMers, target.seedsFromGenome, lengthLeft, 0, 1, kMerNo, maxMismatches + 1,
                                                                biggerGap, variantMerSet, kmerSetPool, growingMerPool, previousVariantType, previousIndelMerNo, indelBalance, threadNo);
#if PERF
                                perfTimers[threadNo][(int)PT.ptTVM2].Stop();
                                if (replacementIdx > 0)
                                    perfCounts[threadNo][(int)PC.pcExt2Matches]++;
#endif
                                foundReplacement = replacementIdx >= 0;
                                if (foundReplacement)
                                    bestReplacement = variantMerSet.merSet[replacementIdx];

                                foundStrongReplacement = foundReplacement && (bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing) >= lengthLeft;
#if DEBUG
                                if (foundStrongReplacement)
                                {
                                    foundIn2xMM++;
                                    whenMatched = "x2";
                                }
#endif
                            } // second attempt

                            // if this is the first mismatch after a reasonable looking run of matches, try for a long gap or even reversal
                            //if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && foundReplacement && !foundStrongReplacement && activeTrailers > trailingCostsLength / 2 && cumulativeCost / activeTrailers < 5)
                            //if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && foundReplacement && !foundStrongReplacement && activeTrailers > trailingCostsLength / 2 && cumulativeCost < 1)
                            if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && foundReplacement && !foundStrongReplacement && trailingMers.Count == trailingCostsLength && cumulativeCost / trailingCosts.Count < 5)
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcExt3Tests]++;
                                perfTimers[threadNo][(int)PT.ptTVM3].Start();
#endif
                                variantMerSet.Reset();
                                noOfVariantMers = GenerateViableMerVariants(MatchType.spanGap, allowedVariants, false, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, growingMerPool, sequence, fbi, lastMatchingFBI, (kMerNo - previousIndelMerNo),
                                                                            currentTargetMers, target.seedsFromGenome, longGap, skipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out skipbackFBI);
#if TRACE
                                variants = "";
                                for (int i = 0; i < variantMerSet.setSize; i++)
                                    variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                                trace.WriteLine(kMerNo + ": TVM3 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif
                                replacementIdx = TryVariantMers(sequence, fbi, kmerSize, MatchType.spanGap, allowedVariants,
                                                                target.mersFromGenomeEither, target.seedsFromGenome, lengthLeft, 0, 1, kMerNo, maxMismatches + 1,
                                                                longGap, variantMerSet, kmerSetPool, growingMerPool, previousVariantType, previousIndelMerNo, indelBalance, threadNo);
#if PERF
                                perfTimers[threadNo][(int)PT.ptTVM3].Stop();
                                if (replacementIdx > 0)
                                    perfCounts[threadNo][(int)PC.pcExt3Matches]++;
#endif

                                foundReplacement = replacementIdx >= 0;
                                bool reversal = false;
                                if (foundReplacement)
                                {
                                    bestReplacement = variantMerSet.merSet[replacementIdx];
                                    MerHashSet priorTargetMers = currentTargetMers;
                                    SetCurrentStrand(bestReplacement.mer, target, out currentStrand, out currentTargetMers);
                                    reversal = priorTargetMers != currentTargetMers;
                                }

                                foundStrongReplacement = foundReplacement && (bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing) >= lengthLeft;
#if DEBUG
                                if (foundStrongReplacement)
                                {
                                    foundIn2xMM++;
                                    if (reversal)
                                        whenMatched = "reversed";
                                    else
                                        whenMatched = "gap";
                                }
#endif
                            } // third attempt
                        } // choose between multiple viable variants

                        // found a sensible looking replacement (not abandoned early due to too many mismatches)
                        if (foundReplacement)
                        {
#if PERF
                            perfCounts[threadNo][(int)PC.pcMatched]++;
#endif
                            // penalise consecutive subs - getting close to rewriting
                            if (previousVariantType == VariantType.sub && bestReplacement.variantType == VariantType.sub)
                                costOfVariant = 2 * (int)Costs.costSub;
                            int droppedCost = 0;
                            if (trailingCosts.Count == trailingCostsLength)
                            {
                                droppedCost = trailingCosts[0];
                                trailingCosts.RemoveAt(0);
                            }
                            cumulativeCost = cumulativeCost - droppedCost + costOfVariant;
                            trailingCosts.Add(costOfVariant);

                            if (cumulativeCost <= maxCost || matchingMode == MatchMode.scanning)
                            {
                                kMer = bestReplacement.mer;
                                lastMatchingMerNo = kMerNo;
                                // adjust next base pointer appropriately (should not move for del, and may have to be moved by multiple bases for ins)
                                fbi += bestReplacement.followingBaseIncrement;
                                previousVariantType = bestReplacement.variantType;
                                if (previousVariantType == VariantType.ins || previousVariantType == VariantType.del)
                                {
                                    previousIndelMerNo = kMerNo;
                                    indelBalance += bestReplacement.indelBalance;
                                }
#if (TRACE || HLTRACE)
                                Mer traceMer = new Mer();
                                traceMer.CopyFrom(bestReplacement);
#endif
                                if (currentStrand == Strand.either)
                                    SetCurrentStrand(kMer, target, out currentStrand, out currentTargetMers);

                                // ending a mismatch region. Extension is much more accomodating to mismatches than scanning, so look back through the mismatch region to pick up any other matches
                                if (inMismatchRegion && scanDirection == ScanType.forwardScan)
                                {
                                    bool tryRescan = (fbi - lastMatchingFBI) >= kmerSize;
                                    if (tryRescan)
                                    {
#if (TRACE || HLTRACE)
                                        trace.WriteLine(kMerNo + ": rescan after fuzzy match");
#endif
#if PERF
                                        perfCounts[threadNo][(int)PC.pcRvsFuzzy]++;
                                        perfTimers[threadNo][(int)PT.ptRVSF].Start();
#endif
                                        matches += RescanMismatchRegion(exactMatchesOnly, allowedVariants, kmerSize, sequence, target, currentStrand, kMer, lastMatchingFBI, fbi, kmerSetPool, growingMerPool, variantMerSet, threadNo);
#if PERF
                                        perfTimers[threadNo][(int)PT.ptRVSF].Stop();
#endif
                                    }
                                }

                                inMismatchRegion = false;
                                lastMatchingFBI = fbi;
                                consecutiveNoMatches = 0;
                                consecutiveFuzzyMatches++;
                                trailingMers.Clear();
#if (TRACE || HLTRACE)
                                trace.WriteLine(kMerNo + ": V\t" + kMers.ExpandMer(kMer, kmerSize) +
                                    " +" + (fbi == sequence.Length ? "end" : sequence.Bases[fbi].ToString()) + " " + currentStrand.ToString() + " " +
                                    traceMer.variantType + " (" + traceMer.gapLength + ")" + " bem=" + traceMer.exactMatchesFollowing + " bmm=" + traceMer.mismatchesFollowing +
                                    " " + matchingMode.ToString() + " " + whenMatched + (scanDirection == ScanType.reverseScan ? ("  [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                                    " cost=" + traceMer.costOfVariant + "/" + traceMer.costOfFollowing + " cc=" + cumulativeCost + " cnm=" + consecutiveNoMatches);
#endif

                                matchingMode = MatchMode.continuing;
                            }
                            else
                            // cumulative cost getting too high, so revert to source sequence again
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcCostly]++;
#endif
#if (TRACE || HLTRACE)
                                trace.WriteLine(kMerNo + ": A\t" + kMers.ExpandMer(kMer, kmerSize) + " " + currentStrand.ToString() + " " +
                                    bestReplacement.variantType + " (" + bestReplacement.gapLength + ")" + " bem=" + bestReplacement.exactMatchesFollowing + " bmm=" + bestReplacement.mismatchesFollowing +
                                    " " + matchingMode.ToString() + " " + whenMatched + (scanDirection == ScanType.reverseScan ? ("  [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                                    " cost=" + bestReplacement.costOfVariant + "/" + bestReplacement.costOfFollowing + " cc=" + cumulativeCost + " cnm=" + consecutiveNoMatches);
#endif
                                foundReplacement = false;
                                // restore the kmer from the sequence
                                fbi = fbi - consecutiveFuzzyMatches;
                                kMerNo = kMerNo - consecutiveFuzzyMatches;
                                Sequence.CondenseMer(sequence, fbi - kmerSize, kmerSize, out kMer);

                            } // cumulative cost too high - abandon replacement
                        } // found a suitable replacement
                    } // need to check for variant

                    // either couldn't get either seed match or an acceptable variant
                    if (!foundReplacement)
                    {
                        //if (kMer == 0x81236186D8220000)
                        //    Debugger.Break();
#if (TRACE || HLTRACE)
                        trace.WriteLine(kMerNo + (checkForVariant ? ": no match " : ": skipped ") + kMers.ExpandMer(kMer, kmerSize) + " " + currentStrand.ToString() + " " +
                                                    (scanDirection == ScanType.reverseScan ? (" [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : ""));
#endif
#if PERF
                        perfCounts[threadNo][(int)PC.pcNoMatch]++;
#endif
                        // if we can't get a match on a reverse scan, just abandon the scan as we've already checked the
                        // rest of the reversed sequence on the forward pass (and nothing matched then)
                        if (scanDirection == ScanType.reverseScan)
                        {
#if (TRACE || HLTRACE)
                            trace.WriteLine(kMerNo + ": stopping reverse scan");
#endif
                            break;
                        }

                        consecutiveNoMatches++;
                        consecutiveFuzzyMatches = 0;
                        if (checkForVariant)
                        {
                            inMismatchRegion = true;
                            matchingMode = MatchMode.scanning;
                            trailingCosts.Clear();
                            trailingMers.Clear();
                            cumulativeCost = 0;
                            currentStrand = Strand.either;
                            currentTargetMers = target.mersFromGenomeEither;
                        }
                    }

                    kMerNo++;

                } // fuzzy matching

            } // for each k-mer in the source sequence

            return matches;
        }

        private static bool CheckForViableVariants(ulong kMer, MerHashSet targetMers, int kmerSize, VariantMerSet variantMerSet)
        {
            // always having the previous exact matching kMer in the list
            int variantsFound = 1;

            ulong baseMask = 0xc000000000000000;
            ulong merWithHole = kMer & ~(baseMask >> (kmerSize * 2 - 2));
            for (ulong b = 0; b <= 3; b++)
            {
                ulong newBase = b << (64 - kmerSize * 2);
                ulong merVariant = merWithHole | newBase;

                if (merVariant == kMer)
                    continue;

                if (targetMers.Contains(merVariant))
                {
                    variantMerSet.AddViable(VariantType.sub, merVariant, 0, 0, (int)Costs.costMatch);
                    variantsFound++;
                }
            }

            return variantsFound > 1;
        }

        private static void SetCurrentStrand(ulong kmer, TiledGenome target, out Strand currentStrand, out MerHashSet currentTargetMers)
        {
            if (target.mersFromGenomePlus.Contains(kmer))
            {
                currentStrand = Strand.plus;
                currentTargetMers = target.mersFromGenomePlus;
            }
            else
            {
                currentStrand = Strand.minus;
                currentTargetMers = target.mersFromGenomeMinus;
            }
#if TRACE
            trace.WriteLine("changing target set to " + currentStrand.ToString());
#endif
        }

        private static bool KeepOnlyViableVariants(VariantMerSet variantMerSet)
        {
            int highestExactMatches = -1;
            int mersInStartingSet = variantMerSet.setSize;
            bool noInsVariant = true;

            for (int m = 0; m < mersInStartingSet; m++)
            {
                if (variantMerSet.merSet[m].exactMatchesFollowing > highestExactMatches)
                    highestExactMatches = variantMerSet.merSet[m].exactMatchesFollowing;
            }

            int nextIdx = 0;

            // condense list of variants down to the still-viable subset
            for (int m = 0; m < mersInStartingSet; m++)
            {
                if (variantMerSet.merSet[m].exactMatchesFollowing == highestExactMatches)
                {
                    variantMerSet.CopyfromMer(m, nextIdx);
                    if (variantMerSet.merSet[m].variantType == VariantType.ins)
                        noInsVariant = false;
                    nextIdx++;
                }
            }

            variantMerSet.setSize = nextIdx;
            return noInsVariant;
        }

        private static int RescanMismatchRegion(bool exactMatching, VariantsAllowed allowedVariants, int kmerSize, Sequence sequence, TiledGenome target, Strand strand, ulong kmer, int lastMatchingFBI, int fbi,
                                                List<VariantMerSet> kmerSetPool, List<GrowingMerSet> growingMerPool, VariantMerSet outerMerSet, int threadNo)
        {
            int matches = 0;

            //// first check that the preceding kmer is viable - regardless of how it can be generated
            //int precedingViableVariants = GenerateViableMerSubVariants(kmer >> 2, kmerSize, outerMerSet, VariantsWanted.firstVariantOnly, target.mersFromGenomeEither);
            //if (precedingViableVariants == 0)
            //    return 0;

            int rescanStart = lastMatchingFBI - kmerSize + 1;
            if (rescanStart < 0)
                rescanStart = 0;

            int rescanLength = fbi - lastMatchingFBI - 1;
            if (rescanLength < kmerSize + 1)
                return 0;

            Sequence mismatchRegion = sequence.SubSeq(rescanStart, rescanLength);
#if (TRACE)
            trace.WriteLine(mismatchRegion.ToString());
#endif
            mismatchRegion.Replace(mismatchRegion.Length - kmerSize + 1, kmer, kmerSize - 1);
            mismatchRegion.ReverseComplement();

#if (TRACE || HLTRACE)
            trace.WriteLine("reverse scan of mismatch region." + " lfbi=" + lastMatchingFBI + " fbi=" + fbi + " length=" + mismatchRegion.Length);
            trace.WriteLine(mismatchRegion.ToString());
#endif
            Strand currentStrand;
            if (strand == Strand.plus)
                currentStrand = Strand.minus;
            else
                currentStrand = Strand.plus;

            matches = MatchSequenceToTarget(exactMatching, allowedVariants, ScanType.reverseScan, kmerSize, mismatchRegion, 0, target, currentStrand, kmerSetPool, outerMerSet, growingMerPool, threadNo);

#if (TRACE || HLTRACE)
            trace.WriteLine("finished reverse mismatch scan - found " + matches + " additional exact matches");
#endif
            return matches;
        }

        private static void ReplaceMerInSeq(ulong changedMer, VariantType variantType, int m, int kmerSize, int basesToDelete, Sequence sequence)
        {
            if (variantType == VariantType.del)
            {
                // the replacement 'del' mer has had a base inserted somewhere. This means that we have to insert 
                // a placeholder base in the read string, prior to replacing the healed mer  

                sequence.Insert(m + kmerSize - 1, 'X');
                sequence.Replace(m, changedMer, kmerSize);
            }

            if (variantType == VariantType.sub)
            {
                // for substitutions, the structure of the healed read doesn't change, 
                // so just overwrite the part of the read corresponding to the mer being healed
                sequence.Replace(m, changedMer, kmerSize);
            }

            if (variantType == VariantType.ins)
            {
                // the replacement 'ins' mer has had a base deleted somewhere. This means that we have to delete 
                // a base from the read string, prior to replacing the healed mer (which includes the following base from the read)
                sequence.Remove(m, basesToDelete);
                sequence.Replace(m, changedMer, kmerSize);
            }
        }

        // 
        private static int TryVariantMers(Sequence sequence, int fbi, int kmerSize, MatchType matchType, VariantsAllowed allowedVariants, MerHashSet targetMers, MerHashSet targetSeeds,
                                          int lengthLeft, int consecutiveFixes, int startingIndent, int kmerNo, int allowedMismatches, int allowedGapLength, VariantMerSet variantMerSet,
                                          List<VariantMerSet> kmerSetPool, List<GrowingMerSet> growingMerPool,
                                          VariantType previousVariantType, int previousIndelMerNo, int indelBalance, int threadNo)
        {
#if TRACE
            string variants = "";
            for (int i = 0; i < variantMerSet.setSize; i++)
                variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + " " + variantMerSet.merSet[i].variantType + " ";
            trace.WriteLine(Indent(startingIndent) + "tvm @ " + kmerNo + " ll=" + lengthLeft + " fbi=" + fbi + " amm=" + allowedMismatches +
                                    " +" + ((fbi == sequence.Length) ? "end" : sequence.Bases[fbi].ToString()) +
                                    " " + previousVariantType.ToString() + " pid=" + previousIndelMerNo + " cf=" + consecutiveFixes + " lg=" + allowedGapLength + " idb=" + indelBalance + " (" + variants + ")");
#endif
            //if (kmerNo == 51 && startingIndent == 1 && allowedMismatches > 5)
            //    Debugger.Break();

            int variantGapLength = allowedGapLength;
            int noOfVariantMers = variantMerSet.setSize;

            List<Mer> variantMers = variantMerSet.merSet;
            int viableMerCount = 0;
            int OKLowestCost = int.MaxValue;
            int OKHighesttExactMatches = 0;
            int countOfBestOKVariants = 0;
            int bestIdx = -1;                       // final 'best' choice
            int bestOKIdx = -1;                     // first of the 'best OK' choices (may only be one so remembered to save searching later)

            // calculate followers and cost for each viable variant in the context of the current sequence
            for (int v = 0; v < noOfVariantMers; v++)
            {
                Mer currentVariant = variantMers[v];
                ulong vMer = currentVariant.mer;
                int indent = startingIndent;

                int exactMatches = 0;
                int mismatches = 0;

                // adjust next base to account for the change type (sub = 0 as fbi already points to next base, del = -1, ins = 1+)
                int adjustedFBI = fbi + currentVariant.followingBaseIncrement;

                int currentIndelMerNo = previousIndelMerNo;
                if (currentVariant.variantType == VariantType.del || currentVariant.variantType == VariantType.ins)
                {
                    currentIndelMerNo = kmerNo;
                    if (currentVariant.followingIndel == VariantType.none)
                        currentVariant.followingIndel = currentVariant.variantType;
                }

                // don't allow indels that are too close together. It's too easy to patch together a match using small indels and 
                // a few matching in-between bases.
                if (((previousVariantType == VariantType.del || previousVariantType == VariantType.ins) &&
                     (currentVariant.variantType == VariantType.del || currentVariant.variantType == VariantType.ins)) &&
                    (kmerNo - previousIndelMerNo) < 4) // gap <= 2
                {
#if TRACE
                    trace.WriteLine(Indent(startingIndent + 1) + "close @ " + kmerNo + " " + kMers.ExpandMer(currentVariant.mer, kmerSize) + " " + previousVariantType +
                                    " @ " + previousIndelMerNo + " (" + currentVariant.variantType + " @" + kmerNo + ")");
#endif
                    currentVariant.variantType = VariantType.invalid;
                    continue;
                }

                // don't allow indels to get out of control. Too many of the same type of indels in a row is probably a sign of going astray.
                if (currentVariant.variantType == VariantType.ins)
                    currentVariant.indelBalance = indelBalance - currentVariant.gapLength;
                if (currentVariant.variantType == VariantType.del)
                    currentVariant.indelBalance = indelBalance + currentVariant.gapLength;
                if (Math.Abs(currentVariant.indelBalance) > maxDownstream / 2)
                {
#if TRACE
                    trace.WriteLine(Indent(startingIndent + 1) + "indel @ " + kmerNo + " " + kMers.ExpandMer(currentVariant.mer, kmerSize) + " " + previousVariantType +
                                    " @ " + previousIndelMerNo + " (" + currentVariant.variantType + " @" + kmerNo + ") " + "cidb=" + currentVariant.indelBalance);
#endif
                    currentVariant.variantType = VariantType.invalid;
                    continue;
                }

                // we've got a viable k-mer variant - one that matches a k-mer in the target set, so find out how well it matches downstream
                CountFollowers(currentVariant, sequence, adjustedFBI, kmerSize, targetMers, allowedVariants, lengthLeft, consecutiveFixes + 1, currentIndelMerNo, indent + 1, kmerNo,
                               allowedMismatches, allowedGapLength, kmerSetPool, threadNo);

                exactMatches = currentVariant.exactMatchesFollowing;
                mismatches = currentVariant.mismatchesFollowing;
                int followers = exactMatches + mismatches;
                int costOfFollowing = currentVariant.costOfFollowing;

                if (followers > 0)
                {
                    // don't choose this path if we appear to be rewriting the downstream region
                    if (currentVariant.variantType != VariantType.none && mismatches > exactMatches && followers > maxDownstream / 2)
                    {
                        currentVariant.variantType = VariantType.invalid;
                        continue;
                    }

                    viableMerCount++;

                    if (followers >= lengthLeft)
                    {
                        // found a winner - all the way to the end with no mismatches 
                        if (exactMatches >= lengthLeft)
                        {
                            noOfVariantMers = v + 1;
                            bestIdx = v;
                            break;
                        }

                        // equal scores to the previous best-OK
                        if (costOfFollowing == OKLowestCost && exactMatches == OKHighesttExactMatches)
                            countOfBestOKVariants++;

                        // possibly a better best-OK
                        if (exactMatches > OKHighesttExactMatches)
                        {
                            OKLowestCost = costOfFollowing;
                            OKHighesttExactMatches = exactMatches;
                            countOfBestOKVariants = 1;
                            bestOKIdx = v;
                        }
                        if (exactMatches == OKHighesttExactMatches && costOfFollowing < OKLowestCost)
                        {
                            OKLowestCost = costOfFollowing;
                            countOfBestOKVariants = 1;
                            bestOKIdx = v;
                        }

                        // stop looking for alternatives after handling all the Subs if we get a satisfactory change (good length at low cost)
                        if (currentVariant.variantType != VariantType.sub && costOfFollowing < -5)
                        {
                            for (int i = v; i < noOfVariantMers; i++)
                                variantMers[i].variantType = VariantType.invalid;
                            noOfVariantMers = v;
                            break;
                        }

                        // if we got to the end of the match region with this variant, any other variants will have to do at least as well to be viable
                        if (allowedMismatches > mismatches)
                            allowedMismatches = mismatches + 1;
                    }
                }
            }

            // if no viable variants, continue the main scanning loop
            if (viableMerCount == 0)
            {
#if TRACE
                trace.WriteLine(Indent(startingIndent) + "tvr @ " + kmerNo + " no viable variants ");
#endif
                return noViableVariants;
            }
#if TRACE
            trace.WriteLine(Indent(startingIndent + 1) + "tvc @ " + kmerNo + " " + " ll=" + lengthLeft + " lc=" + OKLowestCost + " bem=" + OKHighesttExactMatches + " best=" + bestIdx);
            for (int v = 0; v < noOfVariantMers; v++)
                trace.WriteLine(Indent(startingIndent + 1) + "tvc @ " + kmerNo + " " +
                            kMers.ExpandMer(variantMers[v].mer, kmerSize) + " @ " + (fbi + variantMers[v].followingBaseIncrement) + " "
                            + variantMers[v].variantType + "(" + variantMers[v].gapLength + ")" +
                            " bem=" + variantMers[v].exactMatchesFollowing + " bmm=" + variantMers[v].mismatchesFollowing + " fid=" + variantMers[v].followingIndel +
                            " cost=" + variantMers[v].costOfVariant + "/" + variantMers[v].costOfFollowing + " cf=" + consecutiveFixes);
#endif
            // if we didn't find a great choice (all the way with exact matches), we'll settle for the best of the OK variants if there's only one of these
            if (bestIdx < 0 && countOfBestOKVariants == 1)
            {
                bestIdx = bestOKIdx;
            }

            // if there are multiple best OK variants (good length, lowest cost, highest exact matches), just choose the first one encountered (more likely to be Sub)
            if (bestIdx < 0 && countOfBestOKVariants > 1)
            {
                for (int v = 0; v < noOfVariantMers; v++)
                {
                    // ignore any non-viable variants
                    if (variantMers[v].variantType == VariantType.invalid)
                        continue;

                    if (variantMers[v].costOfFollowing == OKLowestCost && variantMers[v].exactMatchesFollowing == OKHighesttExactMatches &&
                        (variantMers[v].exactMatchesFollowing + variantMers[v].mismatchesFollowing >= lengthLeft))
                    {
                        bestIdx = v;
                        break;
                    }
                }
            }

            // nothing with a good length, choose the first one with the highest follower length at the lowest cost
            // This will have the side effect of preferring subs over other equivalent changes.

            if (bestIdx < 0)
            {
                int highestFollowers = 0;
                int lowestCost = int.MaxValue;

                // find the best length at the lowest cost
                for (int v = 0; v < noOfVariantMers; v++)
                {
                    // ignore any non-viable variants
                    if (variantMers[v].variantType == VariantType.invalid)
                        continue;

                    int followers = variantMers[v].exactMatchesFollowing + variantMers[v].mismatchesFollowing;
                    if (followers > highestFollowers)
                    {
                        highestFollowers = followers;
                        lowestCost = variantMers[v].costOfFollowing;
                    }
                    if (followers == highestFollowers && variantMers[v].costOfFollowing < lowestCost)
                        lowestCost = variantMers[v].costOfFollowing;
                }

                // now find the first variant that gives these values
                for (int v = 0; v < noOfVariantMers; v++)
                {
                    // ignore any non-viable variants
                    if (variantMers[v].variantType == VariantType.invalid)
                        continue;

                    int followers = variantMers[v].exactMatchesFollowing + variantMers[v].mismatchesFollowing;
                    if (followers == highestFollowers && variantMers[v].costOfFollowing == lowestCost)
                    {
                        bestIdx = v;
                        break;
                    }
                }
            }
#if TRACE
            trace.WriteLine(Indent(startingIndent) + "tvr @ " + kmerNo + " " +
                            kMers.ExpandMer(variantMers[bestIdx].mer, kmerSize) + " @ " + fbi + " "
                            + variantMers[bestIdx].variantType + "(" + variantMers[bestIdx].gapLength + ")" +
                            " bem=" + variantMers[bestIdx].exactMatchesFollowing + " bmm=" + variantMers[bestIdx].mismatchesFollowing +
                            " cost=" + variantMers[bestIdx].costOfVariant + "/" + variantMers[bestIdx].costOfFollowing + " cfix=" + consecutiveFixes);
#endif

            return bestIdx;
        }

        private static VariantMerSet GetVariantMerSet(List<VariantMerSet> kmerSetPool)
        {
            VariantMerSet variantMerSet;

            if (kmerSetPool.Count == 0)
                variantMerSet = new VariantMerSet(smallVariantSetSize);
            else
            {
                variantMerSet = kmerSetPool[0];
                kmerSetPool.RemoveAt(0);
            }
            //Console.WriteLine("get VMS# " + variantMerSet.merSetNo);
            return variantMerSet;
        }

        private static void ReturnVariantMerSet(List<VariantMerSet> kmerSetPool, VariantMerSet merSet)
        {
            merSet.Reset();
            //Console.WriteLine("return VMS# " + variantMerSet.merSetNo);
            kmerSetPool.Add(merSet);
        }

        private static GrowingMerSet GetGrowingMerSet(List<GrowingMerSet> kmerSetPool)
        {
            GrowingMerSet growingMerSet;

            if (kmerSetPool.Count == 0)
                growingMerSet = new GrowingMerSet(smallVariantSetSize);
            else
            {
                growingMerSet = kmerSetPool[0];
                kmerSetPool.RemoveAt(0);
            }
            //Console.WriteLine("get VMS# " + variantMerSet.merSetNo);
            return growingMerSet;
        }

        private static void ReturnGrowingMerSet(List<GrowingMerSet> kmerSetPool, GrowingMerSet merSet)
        {
            merSet.setSize = 0;
            //Console.WriteLine("return VMS# " + variantMerSet.merSetNo);
            kmerSetPool.Add(merSet);
        }

        // count how many exact matches there are downstream from a given k-mer in a given sequence context. 
        // The scan is length-limited to avoid going to the end of the sequence every time, and mismatch-limited so we can give up early
        private static void CountFollowers(Mer kMer, Sequence sequence, int nbi, int kmerSize, MerHashSet targetMers, VariantsAllowed allowedVariants,
                                           int lengthLeft, int consecutiveFixes, int previousIndelMerNo, int indent, int kMerNo, int allowedMismatches, int allowedGap, List<VariantMerSet> kmerSetPool, int threadNo)
        {
#if TRACE
            trace.WriteLine(Indent(indent) + "cfs @ " + kMerNo + " nbi=" + nbi + " ll=" + lengthLeft + " amm=" + allowedMismatches + " " + kMer.variantType +
                           (kMer.variantType != VariantType.sub ? "(" + kMer.gapLength + ") " : " ") + kMers.ExpandMer(kMer.mer, kmerSize) +
                           " (+" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ") " +
                           "cost=" + kMer.costOfVariant + "/" + kMer.costOfFollowing + " cf=" + consecutiveFixes);
#endif
            ulong nextMer = kMer.mer;
            int startingKMerNo = kMerNo;
            int startingLengthLeft = lengthLeft;
            int consecutiveMatches = 0;
            int exactMatches = 0;
            int costOfFollowing = kMer.costOfVariant;
            int basesToFirstChange = 0;
            int usedMismatches = 0;
            VariantType followingIndel = VariantType.none;

            //// Ins always end in an exact match which has to be counted
            //if (kMer.variantType == VariantType.ins)
            //{
            //    // count the elided matched base
            //    exactMatches = 1;
            //    // and remember we just used up a following base as well 
            //    lengthLeft--;
#if TRACE
            //if (kMer.variantType == VariantType.ins)
            //    trace.WriteLine(Indent(indent) + "cfs(adj) @ " + kMerNo + " nbi=" + nbi + " ll=" + lengthLeft + " amm=" + allowedMismatches + " umm=" + usedMismatches + " " + kMer.variantType +
            //               (kMer.variantType != VariantType.sub ? "(" + kMer.gapLength.ToString() + ") " : " ") + kMers.ExpandMer(kMer.mer, kmerSize) + " (+" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ")");
#endif
            //}

            //if (nextMer == 0x04A34AF613000000)
            //    Debugger.Break();

            while (lengthLeft > 0 && GenerateNextMer(nextMer, kmerSize, sequence, nbi, targetMers, out nextMer))
            {
                if (targetMers.Contains(nextMer))
                {
                    exactMatches++;
                    consecutiveMatches++;
                    nbi++;
                    lengthLeft--;
                    kMerNo++;
                    costOfFollowing += (int)Costs.costMatch;
                    basesToFirstChange++;
                    consecutiveFixes = 0;
#if TRACE
                    trace.WriteLine(Indent(indent + 1) + "cf= @ " + kMerNo + " nbi=" + nbi + " ll=" + lengthLeft + " " + kMers.ExpandMer(nextMer, kmerSize) + " (" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ")");
#endif
                }
                else
                {
                    lengthLeft--;
                    allowedMismatches--;
                    usedMismatches++;
                    nbi++;
                    kMerNo++;

                    if (allowedMismatches > 0 && lengthLeft > 0 && consecutiveFixes < maxConsecutiveFixes)
                    {
#if TRACE
                        //trace.WriteLine(Indent(indent + 1) + "tvm? @ " + kmerNo + " nbi=" + nbi + " ll=" + lengthLeft + " amm=" + allowedMismatches + " tf=" + targetExactMatches + " f=" + exactMatches + " ?=" + (lengthLeft >= targetExactMatches) + " " + kMers.ExpandMer(nextMer, kmerSize) + " (" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ")");
#endif
                        VariantMerSet variantMerSet = GetVariantMerSet(kmerSetPool);
                        int fbiAdjustment;
                        int lastSkipbackFBI = 0;
                        ulong overlookedMer;
                        // this is an 'extension' match, not a 'scan'
                        int noOfMerVariants = GenerateViableMerVariants(MatchType.full, allowedVariants, false, nextMer, kmerSize, false, variantMerSet, null, null, sequence, nbi, nbi, (kMerNo - previousIndelMerNo), targetMers, null, allowedGap, lastSkipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out lastSkipbackFBI);
                        int bestIdx = TryVariantMers(sequence, nbi, kmerSize, MatchType.full, allowedVariants, targetMers, null, lengthLeft, consecutiveFixes, indent + 1, kMerNo, allowedMismatches,
                                                     allowedGap, variantMerSet, kmerSetPool, null, kMer.variantType, previousIndelMerNo, kMer.indelBalance, threadNo);
                        if (bestIdx >= 0)
                        {
                            Mer bestVariant = variantMerSet.merSet[bestIdx];
                            exactMatches += bestVariant.exactMatchesFollowing;
                            usedMismatches += bestVariant.mismatchesFollowing;
                            costOfFollowing += bestVariant.costOfFollowing;
                            if (bestVariant.followingIndel == VariantType.none && (bestVariant.variantType == VariantType.ins || bestVariant.variantType == VariantType.del))
                                followingIndel = bestVariant.variantType;
                            else
                                followingIndel = bestVariant.followingIndel;
                        }
                        ReturnVariantMerSet(kmerSetPool, variantMerSet);
                    }

                    // recursion will have counted downstream from here so just exit the counting loop
                    break;
                }
            }

            // adjust cost to reward runs of exact matches
            costOfFollowing += (consecutiveMatches / 2) * (int)Costs.costMatch;
            //// compensate Ins variants for the exact-matching base subsumed by the variant
            //if (kMer.variantType == VariantType.ins && followingIndel != VariantType.del)
            //    exactMatches++;

            // copy follower results back to variant kMer
            kMer.exactMatchesFollowing = exactMatches;
            kMer.mismatchesFollowing = usedMismatches;
            kMer.costOfFollowing = costOfFollowing;
            kMer.followingIndel = followingIndel;
#if TRACE
            trace.WriteLine(Indent(indent) + "cfr @ " + startingKMerNo + " " + kMer.variantType + "(" + kMer.gapLength + ")"
                + " em=" + exactMatches + " umm=" + usedMismatches + " sll=" + startingLengthLeft +
                " ll=" + lengthLeft + (lengthLeft == 0 ? " end" : "") + " cost=" + kMer.costOfVariant + "/" + kMer.costOfFollowing + " cf=" + consecutiveFixes + " fid=" + followingIndel);
#endif
        }

        private static string Indent(int n)
        {
            return new string(' ', n * 2);
        }

        private static bool GenerateNextMer(ulong mer, int kmerSize, Sequence sequence, int nbi, MerHashSet targetMers, out ulong nextMer)
        {
            nextMer = 0;
            if (nbi >= sequence.Length)
                return false;

            char nextBase = sequence.Bases[nbi];

            long nextBasePacked = kMers.BaseCharToInt(nextBase);

            if (nextBasePacked < 0)
            {
                char chosenBase;
                bool viableBaseFound = GenerateViableNVariant(mer, kmerSize, targetMers, out nextMer, out chosenBase);
                sequence.Bases[nbi] = chosenBase;
                return viableBaseFound;
            }
            else
            {
                nextMer = mer << 2 | (ulong)nextBasePacked << (64 - kmerSize * 2);
                return true;
            }

        }

        private static long GetMersFromSequences(List<string> sequences, int kmerSize, MerHashSet distinctMers)
        {
            DateTime startOfLoad = DateTime.Now;

            distinctMers.Clear();
            long mersInSequences = 0;
            HashSet<char> validBases = new HashSet<char>() { 'A', 'C', 'G', 'T' };

            char[] degenerateBases;
            Dictionary<char, List<char>> degenerateBaseExpansions;
            kMers.InitialiseDegenerateBaseTables(out degenerateBases, out degenerateBaseExpansions);

            for (int c = 0; c < sequences.Count; c++)
            {
                // turn the gene seq into (as-read) k-mers and add them to the hash set for this organism
                List<ulong> kMersFromSeq = new List<ulong>(sequences[c].Length + 100);

                int mersFound = kMers.GenerateExpandedMersFromRead(sequences[c], kmerSize, kMersFromSeq, degenerateBases, degenerateBaseExpansions);

                mersInSequences += sequences[c].Length - kmerSize + 1;

                for (int m = 0; m < mersFound; m++)
                    distinctMers.AddIfNotPresent(kMersFromSeq[m]);

                //int mersFound = kMers.GenerateMersFromRead(sequences[c], kmerSize, ref merSet, ref merValid);
                //mersInSequences += mersFound;
                //for (int m = 0; m < mersFound; m++)
                //{
                //    if (merValid[m])
                //    {
                //        ulong mer = merSet[m];
                //        //if (mer == (ulong)0xa7d19d1624000000)
                //        //    Debugger.Break();
                //        distinctMers.AddIfNotPresent(mer);
                //    }
                //}
            }

            distinctMers.Optimise();
            //Console.WriteLine("Load took " + (DateTime.Now - startOfLoad).TotalSeconds.ToString("#.0") + "s");
            return mersInSequences;
        }

        private static void GetSequences(string genomeFN, int kmerSize, List<string> headers, List<string> sequences, int minLength)
        {
            StreamReader genome = new StreamReader(genomeFN);
            bool EOF = false;

            //if (genomeFN.EndsWith("Cloacimonetes_GCA_002402105.1_ASM240210v1_genomic.fasta"))
            //    Debugger.Break();

            while (!EOF)
            {
                string header;
                string sequence = SeqFiles.ReadRead(genome, SeqFiles.formatFNA, out header);
                if (sequence == null)
                    break;
                if (minLength > 0 && sequence.Length < minLength)
                    continue;
                int startIdx = 0;

                //int targetIdx = sequence.IndexOf("YWRYWMWGCTRGNMRKAATNANWTR");
                //if (targetIdx != -1)
                //    Debugger.Break();

                // remember this sequence so it can be tiled in the matching phase (and break scaffolds into contigs as needed)
                while (true)
                {
                    int nextDegenerateIdx = FindNextDegenerateBase(sequence, startIdx);
                    if (nextDegenerateIdx >= 0)
                    {
                        // check how bad the damage is - we'll pass on just a few degenerate (and expand them)
                        // too many and we'll break the scaffold into contigs
                        int degenerates = int.MaxValue;
                        // only try counting if there's a kMer left in the contig
                        if (nextDegenerateIdx + kmerSize < sequence.Length)
                        {
                            string startOfPossibleBreak = sequence.Substring(nextDegenerateIdx, kmerSize);
                            //if (startOfPossibleBreak == "YWRYWMWGCTRGNMRKAATNANWTR")
                            //    Debugger.Break();
                            degenerates = CountDegenerateBases(startOfPossibleBreak);
                        }
                        if (degenerates > maxDegenerate)
                        {
                            // write what we have so far
                            string contig = sequence.Substring(0, nextDegenerateIdx);
                            if (contig.Length > kmerSize)
                            {
                                sequences.Add(contig);
                                if (headers != null)
                                    headers.Add(header);
                            }

                            // skip until we get a kMer without degenerate bases
                            startIdx = nextDegenerateIdx;
                            while (true)
                            {
                                // make sure there's space for a kMer
                                startIdx++;
                                if (sequence.Length - startIdx < kmerSize)
                                    break;

                                // check degenerate bases in the next kMer
                                string nextMer = sequence.Substring(startIdx, kmerSize);
                                int degeneratesInMer = CountDegenerateBases(nextMer);
                                if (degeneratesInMer == 0)
                                    break;
                            }

                            // if the degenerate region didn't extend to the end to the contig/scaffold
                            // break it here and start scanning/copying again
                            if (startIdx < sequence.Length - kmerSize)
                            {
                                sequence = sequence.Substring(startIdx);
                                startIdx = 0;
                            }
                            else
                                break;
                        }
                        else
                            startIdx = nextDegenerateIdx + 1;
                    }
                    else
                    {
                        // single clean contig left - just add it
                        sequences.Add(sequence);
                        if (headers != null)
                            headers.Add(header);
                        break;
                    }
                }
            }

            Console.WriteLine("loaded contigs from " + genomeFN);
        }

        private static int FindNextDegenerateBase(string sequence, int startIdx)
        {
            int foundIdx = -1;
            HashSet<char> ACGT = new HashSet<char> { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' };
            for (int i = startIdx; i < sequence.Length; i++)
                if (!ACGT.Contains(sequence[i]))
                {
                    foundIdx = i;
                    break;
                }
            return foundIdx;
        }

        private static int CountDegenerateBases(string sequence)
        {
            int counter = 0;
            HashSet<char> ACGT = new HashSet<char> { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' };
            for (int i = 0; i < sequence.Length; i++)
                if (!ACGT.Contains(sequence[i]))
                    counter++;
            return counter;
        }

        private static int GenerateViableMerVariants(MatchType matchType, VariantsAllowed allowedVariants, bool doubleSub, ulong kMer, int kmerSize, bool tryCurrentMer, VariantMerSet variantMerSet, List<VariantMerSet> kmerSetPool, List<GrowingMerSet> seedGrowingPool,
                                                     Sequence sequence, int fbi, int lastMatchingFBI, int distanceToLastIndel, MerHashSet targetMers, MerHashSet targetSeeds, int allowedGapLength, int previousScanbackFBI, int threadNo,
                                                     out int fbiAdjustment, out ulong overlookedMer, out int lastScanbackFBI)
        {
            int subVariantsAdded = 0;
            int insVariantsAdded = 0;
            int delVariantsAdded = 0;
            fbiAdjustment = 0;
            overlookedMer = kMer;
            lastScanbackFBI = previousScanbackFBI;

            if (tryCurrentMer)
                variantMerSet.AddViable(VariantType.none, kMer, 0, 0, (int)Costs.costMatch);

            switch (matchType)
            {
                case MatchType.full:
                    // generate variants in this order - will be order of preference later - and don't allow consecutive indels
                    subVariantsAdded = GenerateViableMerSubVariants(kMer, kmerSize, variantMerSet, VariantsWanted.lastVariantOnly, targetMers);
                    if (allowedVariants == VariantsAllowed.subInsDel && distanceToLastIndel > 1)
                    {
                        delVariantsAdded = GenerateViableMerDelVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength);
                        insVariantsAdded = GenerateViableMerInsVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength);
                    }
                    return variantMerSet.setSize;

                case MatchType.scanning:
                    // if the genomes/genes are small enough, use seed matching/extension during scanning. Seeds have to be both fairly small and sparse
                    if (useSeeds)
                    {
#if PERF
                        perfCounts[threadNo][(int)PC.pcSeedTest]++;
#endif
                        bool seedViable = targetSeeds.Contains(kMer & seedMask);

                        if (seedViable)
                        {
#if PERF
                            perfCounts[threadNo][(int)PC.pcSeedViable]++;
                            perfTimers[threadNo][(int)PT.ptGrowSeed].Start();
#endif
                            subVariantsAdded = GrowSeed(kMer, kmerSize, targetMers, targetSeeds, variantMerSet, seedGrowingPool);
#if PERF
                            if (subVariantsAdded > 0)
                                perfCounts[threadNo][(int)PC.pcSeedGrown]++;
                            perfTimers[threadNo][(int)PT.ptGrowSeed].Stop();
#endif
                        }
                    }
                    // not using seeds so try for a multiple sub (this will only be called every scanCheckInterval bases
                    else
                    {
#if PERF
                        perfCounts[threadNo][(int)PC.pcScanMTests]++;
#endif
                        if (doubleSub)
                            subVariantsAdded = GenerateViableMerTwoSubVariants(kMer, kmerSize, 10, variantMerSet, targetMers);
                        else
                            subVariantsAdded = GenerateViableMerSubVariants(kMer, kmerSize, variantMerSet, VariantsWanted.allSubvariants, targetMers);
#if PERF
                        if (subVariantsAdded > 0)
                            perfCounts[threadNo][(int)PC.pcScanMMatches]++;
                        if (subVariantsAdded > 1)
                            perfCounts[threadNo][(int)PC.pcScanMMatchedM]++;
#endif
                        // if we got a match, we may have skipped over preceding matching kMers, so backtrack and look for any of these
                        if (subVariantsAdded > 0)
                        {
                            ulong skippedMer = 0;
                            int startOfSkippedInterval = fbi - kmerSize - scanCheckInterval;
                            int endOfSkippedInterval = fbi - kmerSize;
#if PERF
                            if (fbi <= previousScanbackFBI)
                                perfCounts[threadNo][(int)PC.pcSkipbackProbs]++;
#endif
                            if (startOfSkippedInterval >= 0 && fbi > previousScanbackFBI)
                            {
                                // find the first non-N kMer in the interval
                                bool kMerOK = false;
                                lastScanbackFBI = fbi;
                                while (startOfSkippedInterval < endOfSkippedInterval)
                                {
                                    kMerOK = Sequence.CondenseMer(sequence, startOfSkippedInterval, kmerSize, out skippedMer);
                                    if (kMerOK)
                                        break;
                                    startOfSkippedInterval++;
                                }

                                if (kMerOK)
                                {
                                    // find the first viable kMer/variant in the skipped interval
                                    bool foundMatch = false;
                                    VariantMerSet btVMS = GetVariantMerSet(kmerSetPool);
                                    int btSubVariants = 0;
                                    while (!foundMatch && startOfSkippedInterval < endOfSkippedInterval)
                                    {
                                        btSubVariants = GenerateViableMerTwoSubVariants(skippedMer, kmerSize, 0, btVMS, targetMers);
                                        if (btSubVariants > 0)
                                            foundMatch = true;
                                        else
                                        {
                                            startOfSkippedInterval++;
                                            kMerOK = kMers.CondenseMerIncremental(kmerSize, skippedMer, sequence.Bases[startOfSkippedInterval + kmerSize - 1], out skippedMer);
                                            if (!kMerOK)
                                                startOfSkippedInterval = endOfSkippedInterval;
                                        }
                                    }

                                    fbiAdjustment = endOfSkippedInterval - startOfSkippedInterval;
                                    overlookedMer = skippedMer;

                                    // if backtrack succeeded in finding matching kMers, return this set
                                    if (fbiAdjustment > 0)
                                    {
                                        variantMerSet.CopyfromMerSet(btVMS);
#if PERF
                                        perfCounts[threadNo][(int)PC.pcScanBackTracks]++;
#endif
#if TRACE
                                        string variants = "";
                                        for (int i = 0; i < variantMerSet.setSize; i++)
                                            variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + " ";
                                        trace.WriteLine("skipped back " + fbiAdjustment + " bases on scan match. " + kMers.ExpandMer(skippedMer, kmerSize) + " --> " + variants);
#endif
                                    }
                                    ReturnVariantMerSet(kmerSetPool, btVMS);
                                }
                            }
                        }
                    }
                    return variantMerSet.setSize;

                case MatchType.spanGap:
                    // try seeing if we've gone astray and need to go back to the source sequence
                    ulong merFromSequence;
                    int seqMerVariantsAdded = 0;
                    if (Sequence.CondenseMer(sequence, fbi - kmerSize, kmerSize, out merFromSequence))
                    {
                        if (targetMers.Contains(merFromSequence))
                        {
                            variantMerSet.AddViable(VariantType.sub, merFromSequence, 0, 0, (int)Costs.costSub);
                            seqMerVariantsAdded = 1;
                        }
                        seqMerVariantsAdded = GenerateViableMerTwoSubVariants(merFromSequence, kmerSize, 0, variantMerSet, targetMers);
                    }

                    subVariantsAdded = GenerateViableMerSubVariants(kMer, kmerSize, variantMerSet, VariantsWanted.allSubvariants, targetMers);
                    if (allowedVariants == VariantsAllowed.subInsDel && distanceToLastIndel > 1)
                    {
                        delVariantsAdded = GenerateViableMerDelVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength);
                        insVariantsAdded = GenerateViableMerInsVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength);
                    }
                    return variantMerSet.setSize;
            }
            return 0;
        }

        private static int GrowSeed(ulong mer, int kmerSize, MerHashSet targetMers, MerHashSet targetSeeds, VariantMerSet variantMerSet, List<GrowingMerSet> growingMerPool)
        {
            // starts with a kMer with its start (the seed) matching a short starting kMer in the target set.
            // The seed kMer is extended base-by-base using the target seed set. Mismatches against the starting kMer
            // are tracked and branches are abandoned when the mismatch count gets too high.

            GrowingMerSet growingMers = GetGrowingMerSet(growingMerPool);

            growingMers.Add(seedLength - 1, mer, 0, false);

            // grow/extend the seed kMer until it is a full kMer. previous generation is g-1, generation being added is g
            for (int g = seedLength; g < kmerSize; g++)
            {
                ulong currentGenBaseMask = ~(0xc000000000000000 >> g * 2);
                int currentGenSeedShift = (g - seedLength + 1) * 2;
                int currentGenIdx = growingMers.setSize;

                // iterate over all growing kMers in the 'current' generation (g-1)
                for (int v = 0; v < currentGenIdx; v++)
                {
                    GrowingMer gmv = growingMers.merSet[v];

                    //Console.WriteLine(gmv.gen + "-" + v + ": " + kMers.ExpandMer(gmv.mer, kmerSize));

                    // try the next seed in the current kMer - no need to worry about too many mismatches in this case
                    ulong seedFromMer = (gmv.mer << currentGenSeedShift) & seedMask;
                    if (targetSeeds.Contains(seedFromMer))
                    {
                        growingMers.Add(g, gmv.mer, gmv.mismatches, false);
                        //Console.WriteLine(g + " A: " + kMers.ExpandMer(gmv.mer, kmerSize));
                        // and don't look for alternative bases if the next one in the kMer is OK. Continuing to look for variants in this case can lead to tree explosions.
                        //continue;
                    }

                    // don't try variants if we've reached max mismatches for this path
                    if (gmv.mismatches == maxSeedMismatches)
                        continue;

                    // try substituting all 4 possible bases (and skip the non-variant we just checked)
                    for (ulong b = 0; b < 4; b++)
                    {
                        ulong variant = (gmv.mer & currentGenBaseMask) | (b << (32 - g - 1) * 2);
                        ulong seedFromVariant = (variant << currentGenSeedShift) & seedMask;

                        //Console.WriteLine(kMers.ExpandMer(variant, kmerSize));
                        //Console.WriteLine(kMers.ExpandMer(seedFromVariant, seedLength));

                        // skip the non-variant (already handled above)
                        if (seedFromVariant == seedFromMer)
                            continue;

                        // if we found the seed from this extended kMer, add it to the list for the next generation
                        if (targetSeeds.Contains(seedFromVariant))
                        {
                            growingMers.Add(g, variant, gmv.mismatches + 1, false);
                            //Console.WriteLine(g + " V: " + kMers.ExpandMer(variant, kmerSize));
                        }
                    }
                }

                // compress the growingMers list (keep next generation [g] and discard/overwrite the previous generation [g-1]
                int condensedIdx = 0;
                bool foundMatchingMer = false;
                for (int v = 0; v < growingMers.setSize; v++)
                {
                    GrowingMer ngmv = growingMers.merSet[v];
                    // current generation, so overwrite the earliest remaining previous generation
                    if (ngmv.gen == g)
                    {
                        growingMers.CopyfromMer(v, condensedIdx);
                        bool merMatched = targetMers.Contains(ngmv.mer);
                        growingMers.merSet[condensedIdx].matchedFullMer = merMatched;
                        foundMatchingMer = foundMatchingMer | merMatched;
                        condensedIdx++;
                    }
                }
                growingMers.setSize = condensedIdx;

                //Console.WriteLine(currentGenIdx + " viable");

                // break if there are no still-growing kMers
                if (currentGenIdx == 0 || foundMatchingMer)
                    break;
            }

            int variantsAdded = 0;
            for (int v = 0; v < growingMers.setSize; v++)
            {
                if (growingMers.merSet[v].matchedFullMer)
                {
                    variantMerSet.AddViable(VariantType.sub, growingMers.merSet[v].mer, 0, 0, (int)Costs.costSub);
                    variantsAdded++;
                }
            }

            ReturnGrowingMerSet(growingMerPool, growingMers);

            return variantsAdded;
        }

        private static void DumpGMV(List<GrowingMer> growingMers, int kmerSize)
        {
            foreach (GrowingMer gm in growingMers) Console.WriteLine(kMers.ExpandMer(gm.mer, kmerSize));
        }

        private static void DumpSeeds(MerHashSet seeds, int kmerSize)
        {
            StreamWriter seedsDump = new StreamWriter("seeds.txt");
            foreach (ulong mer in seeds)
                seedsDump.WriteLine(kMers.ExpandMer(mer, kmerSize));
            seedsDump.Close();
        }

        private static string DumpVMS(VariantMerSet vms, int kmerSize)
        {
            string result = "";
            for (int i = 0; i < vms.setSize; i++)
                result += i + "\t" + kMers.ExpandMer(vms.merSet[i].mer, kmerSize) + "\n";
            return result;
        }

        // generate all single-sub variants of the kMer (excluding the starting kMer)
        private static int GenerateViableMerSubVariants(ulong kMer, int kmerSize, VariantMerSet variantMerSet, VariantsWanted variantsWanted, MerHashSet targetMers)
        {
            int start = 0;
            int last = kmerSize;
            if (variantsWanted == VariantsWanted.lastVariantOnly)
                start = kmerSize - 1;
            if (variantsWanted == VariantsWanted.firstVariantOnly)
                last = 1;
            int variantsAdded = 0;

            ulong baseMask = 0xc000000000000000;

            for (int m = start; m < last; m++)
            {
                ulong kMerWithHole = kMer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong kMerVariant = kMerWithHole | newBase;

                    if (kMerVariant == kMer)
                        continue;

                    if (targetMers.Contains(kMerVariant))
                    {
                        variantMerSet.AddViable(VariantType.sub, kMerVariant, 0, 0, (int)Costs.costSub);
                        variantsAdded++;
                    }
                }
            }
            return variantsAdded;
        }

        private static int GenerateViableMerTwoSubVariants(ulong mer, int kmerSize, int from, VariantMerSet variantMerSet, MerHashSet targetMers)
        {
            int variantsAdded = 0;
            ulong baseMask = 0xc000000000000000;
            variantMerSet.distinctMers.Clear();

            // vary first of the two bases
            for (int m1 = from; m1 < kmerSize; m1++)
            {
                ulong firstBaseMask = baseMask >> (m1 * 2);
                ulong startingMerWithHole = mer & ~firstBaseMask;
                ulong firstStartingBase = mer & firstBaseMask;
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong firstVariantBase = b << (64 - (m1 + 1) * 2);
                    ulong firstMerVariant = startingMerWithHole | firstVariantBase;

                    // vary the second of the bases
                    for (int m2 = m1 + 1; m2 < kmerSize; m2++)
                    {
                        ulong secondBaseMask = baseMask >> (m2 * 2);
                        ulong variantMerWithHole = firstMerVariant & ~secondBaseMask;
                        ulong secondStartingBase = mer & secondBaseMask;
                        for (ulong b2 = 0; b2 <= 3; b2++)
                        {
                            ulong secondVariantBase = b2 << (64 - (m2 + 1) * 2);
                            ulong secondMerVariant = variantMerWithHole | secondVariantBase;

                            if (!variantMerSet.distinctMers.Contains(secondMerVariant))
                            {
                                variantMerSet.distinctMers.Add(secondMerVariant);
                                if (secondMerVariant != mer && targetMers.Contains(secondMerVariant))
                                {
                                    variantMerSet.AddViable(VariantType.sub, secondMerVariant, 0, 0, (int)Costs.costSub);
                                    variantsAdded++;
                                }
                            }

                        }
                    } // second varying base
                } // first varying base
            }

            return variantsAdded;
        }

        private static int GenerateViableMerThreeSubVariants(ulong mer, int kmerSize, VariantMerSet variantMerSet, MerHashSet targetMers)
        {
            int variantsAdded = 0;
            ulong baseMask = 0xc000000000000000;
            variantMerSet.distinctMers.Clear();

            // vary first of the three bases
            for (int m1 = 0; m1 < kmerSize; m1++)
            {
                ulong firstBaseMask = baseMask >> (m1 * 2);
                ulong firstMerWithHole = mer & ~firstBaseMask;
                ulong firstStartingBase = mer & firstBaseMask;
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong firstVariantBase = b << (64 - (m1 + 1) * 2);
                    ulong firstMerVariant = firstMerWithHole | firstVariantBase;

                    // vary the second of the bases
                    for (int m2 = m1 + 1; m2 < kmerSize; m2++)
                    {
                        ulong secondBaseMask = baseMask >> (m2 * 2);
                        ulong secondMerWithHole = firstMerVariant & ~secondBaseMask;
                        ulong secondStartingBase = mer & secondBaseMask;
                        for (ulong b2 = 0; b2 <= 3; b2++)
                        {
                            ulong secondVariantBase = b2 << (64 - (m2 + 1) * 2);
                            ulong secondMerVariant = secondMerWithHole | secondVariantBase;

                            // vary the third of the bases
                            for (int m3 = m2 + 1; m3 < kmerSize; m3++)
                            {
                                ulong thirdBaseMask = baseMask >> (m3 * 2);
                                ulong thirdMerWithHole = firstMerVariant & ~secondBaseMask;
                                ulong thirdStartingBase = mer & secondBaseMask;
                                for (ulong b3 = 0; b3 <= 3; b3++)
                                {
                                    ulong thirdVariantBase = b3 << (64 - (m3 + 1) * 2);
                                    ulong thirdMerVariant = thirdMerWithHole | thirdVariantBase;

                                    if (!variantMerSet.distinctMers.Contains(thirdMerVariant))
                                    {
                                        variantMerSet.distinctMers.Add(thirdMerVariant);
                                        if (thirdMerVariant != mer && targetMers.Contains(thirdMerVariant))
                                        {
                                            variantMerSet.AddViable(VariantType.sub, thirdMerVariant, 0, 0, (int)Costs.costSub);
                                            variantsAdded++;
                                        }
                                    }
                                }
                            }
                        } // third varying base
                    } // second varying base
                } // first varying base
            }

            return variantsAdded;
        }

        // Inserts bases to produce viable variants. Multiple del variants can be generated, corresponding to filling a gap. 
        // All del variants returned will be followed by an exact match.
        private static int GenerateViableMerDelVariants(ulong mer, int kmerSize, Sequence sequence, int fbi, VariantMerSet variantMerSet, MerHashSet targetMers, int allowedGapLength)
        {
            int variantsAdded = 0;
            ulong upperMerMask = 0xffffffffffffffff << (64 - (kmerSize - 1) * 2);
            int firstDelVariantAdded = variantMerSet.setSize;
            int startingDelVariant = firstDelVariantAdded;
            ulong nextBaseShifted = ulong.MaxValue;
            if (fbi < sequence.Length)
            {
                char nextBase = sequence.Bases[fbi - 1];
                nextBaseShifted = (ulong)kMers.BaseCharToInt(nextBase) << (64 - kmerSize * 2);
            }

            // add viable single-insertion variants
            AddViableDelVariants(mer & upperMerMask, 1, kmerSize, targetMers, variantMerSet);

            // add extended variants of the del variants already added
            for (int g = 1; g < allowedGapLength; g++)
            {
                // remember last one added
                int lastDelVariantAdded = variantMerSet.setSize;

                for (int v = firstDelVariantAdded; v < lastDelVariantAdded; v++)
                    AddViableDelVariants(variantMerSet.merSet[v].mer << 2, (g + 1), kmerSize, targetMers, variantMerSet);

                firstDelVariantAdded = lastDelVariantAdded;
            }

            // now condense the Del variants by deleting any without a following exact match
            int nextViableDelIdx = startingDelVariant;
            for (int v = startingDelVariant; v < variantMerSet.setSize; v++)
            {
                ulong followingMer = variantMerSet.merSet[v].mer << 2 | nextBaseShifted;
                if (targetMers.Contains(followingMer))
                {
                    variantMerSet.merSet[nextViableDelIdx].CopyFrom(variantMerSet.merSet[v]);
                    nextViableDelIdx++;
                    variantsAdded++;
                }
            }
            variantMerSet.setSize = nextViableDelIdx;
            //for (int i = 0; i < variantMerSet.setSize; i++)
            //    trace.WriteLine(variantMerSet.merSet[i].variantType.ToString() + " " + kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize));

            return variantsAdded;
        }

        private static int AddViableDelVariants(ulong maskedMer, int gap, int kmerSize, MerHashSet targetMers, VariantMerSet variantMerSet)
        {
            int variantsAdded = 0;

            for (ulong b = 0; b <= 3; b++)
            {
                ulong newBase = b << (64 - (kmerSize * 2));
                ulong merVariant = maskedMer | newBase;

                // see if the generated kMer (with the inserted base) matches the target
                if (targetMers.Contains(merVariant))
                {
                    variantMerSet.AddViable(VariantType.del, merVariant, -1, gap, gap > 1 ? (int)Costs.costGapExtend : (int)Costs.costGapOpen);
                    variantsAdded++;
                }
            }

            return variantsAdded;
        }

        // deletes bases to produce viable variants.
        private static int GenerateViableMerInsVariants(ulong mer, int kmerSize, Sequence sequence, int fbi, VariantMerSet variantMerSet, MerHashSet targetMers, int allowedGapLength)
        {
            int variantsAdded = 0;

            int consecutiveIns = 0;
            for (int g = 0; g < allowedGapLength; g++)
            {
                char followingBase;
                if (fbi < sequence.Length)
                {
                    followingBase = sequence.Bases[fbi];
                    ulong newBase = (ulong)kMers.BaseCharToInt(followingBase) << (64 - kmerSize * 2);
                    ulong merWithoutLastBaseMask = (ulong)0xfffffffffffffff << ((64 - (kmerSize) * 2) + 2);
                    ulong merWithoutLastBase = mer & merWithoutLastBaseMask;
                    ulong merVariant = merWithoutLastBase | newBase;
                    consecutiveIns++;
                    if (targetMers.Contains(merVariant))
                    {
                        variantMerSet.AddViable(VariantType.ins, merVariant, consecutiveIns, g + 1, g > 0 ? (int)Costs.costGapExtend : (int)Costs.costGapOpen);
                        variantsAdded++;
                    }
                }
                fbi++;
            }

            return variantsAdded;

        }

        // Takes a previous (matched) k-mer that is followed by N in the sequence and produces a viable variant if possible (for GetNextMer)
        private static bool GenerateViableNVariant(ulong previousMer, int kmerSize, MerHashSet targetMers, out ulong viableMer, out char chosenBase)
        {
            ulong partialMer = previousMer << 2;                // shift across 1 base to make room for the variants
            viableMer = partialMer;                             // default 'viable' mer has rhs base set to A
            chosenBase = 'A';
            bool foundViableMer = false;

            // try all possible variants, stopping when we find a viable one. There could be multiple viable kmers but we'll take the first one (and all upstream paths will suffer equally)
            for (ulong b = 0; b < 4; b++)
            {
                ulong sb = b << (64 - kmerSize * 2);
                ulong merVariant = partialMer | sb;
                if (targetMers.Contains(merVariant))
                {
                    viableMer = merVariant;
                    foundViableMer = true;
                    chosenBase = kMers.baseToChar[b];
                    break;
                }
            }

            return foundViableMer;
        }

        // Takes a previous k-mer that is followed by N in the sequence and produces a *non-viable* variant if possible (for subsequent correction).
        // Return N-->A if all are viable (unlikely).
        private static ulong GenerateNonViableNVariant(ulong previousMer, int kmerSize, MerHashSet targetMers, out char chosenBase)
        {
            ulong partialMer = previousMer << 2;                    // shift across 1 base to make room for the variants
            ulong nonViableMer = partialMer;                        // default mer has rhs base set to A
            chosenBase = 'A';

            // try all possible variants, stopping when we find a non-viable one
            for (ulong b = 0; b < 4; b++)
            {
                ulong sb = b << (64 - kmerSize * 2);
                ulong merVariant = partialMer | sb;
                if (!targetMers.Contains(merVariant))
                {
                    nonViableMer = merVariant;
                    chosenBase = kMers.baseToChar[b];
                    break;
                }
            }

            return nonViableMer;
        }

        private static int CountConsecutiveNs(Sequence sequence, int lastBaseIdx)
        {
            int countNs = 0;
            while (lastBaseIdx < sequence.Length && kMers.BaseCharToInt(sequence.Bases[lastBaseIdx]) == -1)
            {
                lastBaseIdx++;
                countNs++;
            }

            return countNs;
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
    }

    public class Mer
    {
        public ulong mer;
        public VariantType variantType;
        public int exactMatchesFollowing;
        public int mismatchesFollowing;
        public int followingBaseIncrement;
        public int gapLength;
        public int indelBalance;
        public int costOfVariant;
        public int costOfFollowing;
        public VariantType followingIndel;

        internal void CopyFrom(Mer mer)
        {
            this.mer = mer.mer;
            this.variantType = mer.variantType;
            this.exactMatchesFollowing = mer.exactMatchesFollowing;
            this.mismatchesFollowing = mer.mismatchesFollowing;
            this.followingBaseIncrement = mer.followingBaseIncrement;
            this.gapLength = mer.gapLength;
            this.costOfVariant = mer.costOfVariant;
            this.costOfFollowing = mer.costOfFollowing;
            this.followingIndel = VariantType.none;
        }
    }

    public class VariantMerSet
    {
        //public int merSetNo;
        public int setSize;
        public List<Mer> merSet;
        public HashSet<ulong> distinctMers;

        public VariantMerSet(int capacity)
        {
            setSize = 0;
            merSet = new List<Mer>(capacity);
            distinctMers = new HashSet<ulong>();
        }

        public void Reset()
        {
            setSize = 0;
            //distinctMers.Clear();
        }

        public void AddViable(VariantType variantType, ulong mer, int followingBaseIncrement, int gapLength, int cost)
        {
            if (setSize == merSet.Count)
            {
                Mer newMer = new Mer();
                merSet.Add(newMer);
            }

            merSet[setSize].variantType = variantType;
            merSet[setSize].mer = mer;
            merSet[setSize].followingBaseIncrement = followingBaseIncrement;
            merSet[setSize].exactMatchesFollowing = 0;
            merSet[setSize].mismatchesFollowing = 0;
            merSet[setSize].gapLength = gapLength;
            merSet[setSize].costOfVariant = cost;
            merSet[setSize].costOfFollowing = 0;
            setSize++;
        }

        public void CopyfromMer(int fromIdx, int i)
        {
            merSet[i].variantType = merSet[fromIdx].variantType;
            merSet[i].mer = merSet[fromIdx].mer;
            merSet[i].followingBaseIncrement = merSet[fromIdx].followingBaseIncrement;
            merSet[i].exactMatchesFollowing = merSet[fromIdx].exactMatchesFollowing;
            merSet[i].mismatchesFollowing = merSet[fromIdx].mismatchesFollowing;
        }

        public void CopyfromMerSet(VariantMerSet fromSet)
        {
            if (merSet.Count < fromSet.merSet.Count)
                for (int i = merSet.Count; i < fromSet.merSet.Count; i++)
                    merSet.Add(new Mer());

            for (int i = 0; i < fromSet.setSize; i++)
            {
                merSet[i].variantType = fromSet.merSet[i].variantType;
                merSet[i].mer = fromSet.merSet[i].mer;
                merSet[i].followingBaseIncrement = fromSet.merSet[i].followingBaseIncrement;
                merSet[i].exactMatchesFollowing = fromSet.merSet[i].exactMatchesFollowing;
                merSet[i].mismatchesFollowing = fromSet.merSet[i].mismatchesFollowing;
            }
            setSize = fromSet.setSize;
        }
    }

    public class GrowingMerSet
    {
        //public int merSetNo;
        public int setSize;
        public List<GrowingMer> merSet;

        public GrowingMerSet(int capacity)
        {
            setSize = 0;
            merSet = new List<GrowingMer>(capacity);
            //merSetNo = Interlocked.Increment(ref nextVariantSetNo);
        }

        public void Add(int gen, ulong mer, int mismatches, bool matchedFullMer)
        {
            if (setSize == merSet.Count)
            {
                GrowingMer newMer = new GrowingMer();
                merSet.Add(newMer);
            }

            merSet[setSize].gen = gen;
            merSet[setSize].mer = mer;
            merSet[setSize].mismatches = mismatches;
            merSet[setSize].matchedFullMer = matchedFullMer;

            setSize++;
        }

        public void CopyfromMer(int fromIdx, int i)
        {
            merSet[i].gen = merSet[fromIdx].gen;
            merSet[i].mer = merSet[fromIdx].mer;
            merSet[i].mismatches = merSet[fromIdx].mismatches;
            merSet[i].matchedFullMer = merSet[fromIdx].matchedFullMer;
        }
    }

    public class GrowingMer
    {
        public int gen;
        public ulong mer;
        public int mismatches;
        public bool matchedFullMer;
    }

    public class MatchingRequest
    {
        public int source;
        public int target;

        public MatchingRequest(int s, int t)
        {
            source = s;
            target = t;
        }
    }

    public class MatchingThreadParams
    {
        public int threadNo;
        public bool exactMatch;
        public VariantsAllowed allowedVariants;
        public int kmerSize;
        public Queue<MatchingRequest> queuedRequests;
        public List<string> sourceSeqs;
        public List<string> targetFNs;
        public List<string> targetSeqs;
        public string[] targetNames;
        public long[] targetSizes;
        public int minLength;
        public long[,] matches;
        public long[,] timesForMatches;
    }

    public class TiledGenome
    {
        public int refCount;
        public MerHashSet mersFromGenomePlus;
        public MerHashSet mersFromGenomeMinus;
        public MerHashSet mersFromGenomeEither;
        public MerHashSet seedsFromGenome;

        public TiledGenome(int genomesCount)
        {
            this.refCount = genomesCount;
            this.mersFromGenomePlus = null;
            this.mersFromGenomeMinus = null;
            this.mersFromGenomeEither = null;
            this.seedsFromGenome = null; 
        }
    }

}


