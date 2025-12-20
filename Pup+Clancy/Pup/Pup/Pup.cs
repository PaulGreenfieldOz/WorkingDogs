using System.Diagnostics;
using WorkingDogsCore;
using PupClancy;
using static PupClancy.SharedCode;

namespace Pup
{
    // Cross-compares (at a k-mer level) the set of sequences from a set of genes/genomes/draft genomes (from an assembler or NCBI/GenBank).
    // All files are assumed to be in FASTA format.
    //
    // The input (target) sequences can come in 3 confusing forms:
    // -dirs        a set of directories, each containing a set of files containing one or more sequences
    //              a directory for each genome, and a file for each contig/sequence (such as a .fna file)
    // -files       a set of files, each containing a number of sequences (each file being a genome/draft genome)
    //              one file per genome, each containing a number of sequences
    // -single      a single file containing a number of sequences (used for gene-level comparison)
    //              sequences within the file are compared
    //
    // The comparisons can be exact or single-base variants (inexact/fuzzy), and either just subs or subs/dels/ins.
    //
    // The default behaviour is an NxN comparison from all the (target) genomes onto the same target genomes.
    // The -source option allows for the source and destinations to be different, useful if matching a set of genomes against a set of reference genomes.
    //
    // The result is a matrix of each organism compared to each other organism in 3 formats: counts, symmetric similarity PIM, sparse similarity PIM
    //
    // -fuzzy -k 25 -t 24 -f GxG_test_1.txt GxG_test.fa
    // -fuzzy -k 20 -t 16 -f Hancockii_GxGindel_20.txt *.fa
    // -fuzzy -k 20 -t 16 -s alkB_fuzzy_20.txt alkBFromGenBankDB.fa
    // -fuzzy -k 20 -t 16 -d cholera_fuzzy_20.txt *.fa
    // -fuzzy -k 20 -t 2 -mm  5 -f c23o_GxG_v13txt  c23oTests.fa
    // -fuzzy -k 20 -t 2 -mm  5  -f 846-858_GxG_f20_v13.txt 846-858.fa
    // -k 20 -f -fuzzy -t 4 -out c23o_GxG_fz_V24_GrowSeed-L10M5-3-K20-SpS.txt c23oTests.fa
    // -k 20 -f -fuzzy -t 4 -out Acido_GxG_fz_V24_GrowSeed-L10M5-ES-K20.txt Acido*.fasta -v
    // -k 20 -d -subs -t 5 -v -out Solanums_GxG_exact_k20.txt *.fsa_nt
    // -k 21 -t 16 -inexact -sf -out c23o_v131_k21_inexact_GxG.txt c23oTests.fa -mm 10
    // -k 21 -t 16 -inexact -files -out Mc_v131_k21_inexact_mm10_GxG.txt *.fna -mm 10 -quiet
    // -k 21 -t 16 -inexact -files -out Mc_bourg_MAB_inexact_mm10_GxG.txt Methanoculleus_bourgensis_GCF_000304355.2_Mb_MS2_genomic.fna Methanoculleus_sp._MAB1_GCF_900036045.1_Methanoculleus_sp_MAB1_genomic.fna -mm 10 -sl 10 -sm 5 -quiet

    class Program
    {
        public const string release = "v1.4.2";

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

            string command = "Pup " + release;
            for (int p = 0; p < args.Length; p++)
                command += " " + args[p];

            string sourceFlag = "-f";
            SearchOption so = SearchOption.TopDirectoryOnly;
            bool exactMatch = false;
            VariantsAllowed allowedVariants = VariantsAllowed.subInsDel;
            int kmerSize = 21;
            int noThreads = Environment.ProcessorCount / 2;
            List<string> targetFNParams = new List<string>();
            string sourceFNArg = null;
            string countsFN = null;
            string countsPrefix = null;
            string countsSuffix = ".txt";
            bool verbose = false;
            bool quiet = false;
            float minPIM = 0.0f;
            int minLength = 0;
            int slParam = 0;
            int smParam = 0;

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
                            kmerSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        if (kmerSize > 31 || kmerSize <= 0)
                        {
                            Console.WriteLine("kMer length must be 0>=k<= 31");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-d" || args[p] == "-dirs")
                    {
                        sourceFlag = "-d";
                        continue;
                    }

                    if (args[p] == "-f" || args[p] == "-files")
                    {
                        sourceFlag = "-f";
                        continue;
                    }

                    if (args[p] == "-sf" || args[p] == "-singlefile" || args[p] == "-genes")
                    {
                        sourceFlag = "-s";
                        continue;
                    }

                    if (args[p] == "-s")
                    {
                        so = SearchOption.AllDirectories;
                        continue;
                    }

                    if (args[p] == "-source" || args[p] == "-sources")
                    {
                        sourceFNArg = args[p + 1];
                        p++;
                        continue;
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
                            slParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -sl|-seedLength parameter: " + args[p + 1]);
                            return;
                        }

                        if (slParam == 0)
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
                            smParam = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the -sm|-seedMismatch parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-min" || args[p] == "-minPIM")
                    {
                        if (!CheckForParamValue(p, args.Length, "Number expected after -min|-minPIM"))
                            return;
                        try
                        {
                            minPIM = Convert.ToSingle(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("Expected a number for the --min|-minPIM parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-out" || args[p] == "-output" || args[p] == "-o")
                    {
                        if (!CheckForParamValue(p, args.Length, "file name expected after -out"))
                            return;

                        countsFN = args[p + 1];
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

                    if (args[p] == "-inexact" || args[p] == "-fuzzy")
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

                    if (args[p] == "-v" || args[p] == "-verbose")
                    {
                        verbose = true;
                        continue;
                    }

                    if (args[p] == "-q" || args[p] == "-quiet")
                    {
                        quiet = true;
                        continue;
                    }

                    Console.WriteLine("Unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                targetFNParams.Add(args[p]);
            }

            smallVariantSetSize = 3 + 1 + 4;
            largeVariantSetSize = 3 * kmerSize + kmerSize + 4 * kmerSize;

            List<string> sourceFNs = new List<string>();
            if (sourceFNArg != null)
            {
                string path;
                string FNP;
                GetPathFN(sourceFNArg, out path, out FNP);
                string[] foundFNs = Directory.GetFiles(path, FNP, SearchOption.TopDirectoryOnly);
                if (foundFNs.Length == 0)
                    Console.WriteLine("No files found matching " + sourceFNArg);

                foreach (string foundFN in foundFNs)
                    sourceFNs.Add(foundFN);
            }

            List<string> targetGenomeFNs = new List<string>();
            foreach (string FN in targetFNParams)
            {
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
                    targetGenomeFNs.Add(foundFN);
            }

            int noOfTargetGenomes = targetGenomeFNs.Count;
            if (noOfTargetGenomes == 0)
            {
                Console.WriteLine("no seq files found to compare");
                return;
            }

            targetGenomeFNs.Sort();

#if (LLTRACE || HLTRACE)
            SharedCode.trace = new StreamWriter("trace.txt");
#endif
            DateTime comparisonStart = DateTime.Now;

            // source sequences - will be same as target sequences for usual GxG comparisons
            // a nice form of the source sequence file names
            string[] sourceGenomeNames;
            // and how many valid k-mers there were in each genome
            long[] sourceGenomeSizes;
            // and the actual sequences/contigs from each genome
            List<string>[] sourceGenomeSequences;

            // and the same for the target genomes
            // a nice form of the sequence file names
            string[] targetGenomeNames;
            // and how many valid k-mers there were in each genome
            long[] targetGenomeSizes;
            // and the actual sequences/contigs from each genome
            List<string>[] targetGenomeSequences;

            char[] FNSeparator = new char[] { Path.DirectorySeparatorChar };
            Dictionary<string, int> directoryNames = new Dictionary<string, int>();
            int nextFileIdx = 0;

            bool betweenFiles = true;
            if (sourceFlag == "-f" && targetGenomeFNs.Count == 1)
            {
                Console.WriteLine("Only one file found. Comparing within this file");
                betweenFiles = false;
            }
            if (sourceFlag == "-s")
                betweenFiles = false;

            if (betweenFiles)
            {
                if (sourceFlag == "-d")
                {
                    int organismIdx = 0;
                    foreach (string genomeFN in targetGenomeFNs)
                    {
                        string[] parsedFN = genomeFN.Split(FNSeparator);
                        string dirName = parsedFN[parsedFN.Length - 2];
                        if (!directoryNames.ContainsKey(dirName))
                        {
                            directoryNames.Add(dirName, organismIdx);
                            organismIdx++;
                        }
                    }
                    noOfTargetGenomes = directoryNames.Count;
                }

                if (sourceFlag == "-f")
                {
                    noOfTargetGenomes = targetGenomeFNs.Count;
                }

                targetGenomeSequences = new List<string>[noOfTargetGenomes];
                targetGenomeNames = new string[noOfTargetGenomes];
                targetGenomeSizes = new long[noOfTargetGenomes];

                // read all the sequences from each genome 
                foreach (string genomeFN in targetGenomeFNs)
                {
                    if (sourceFlag == "-d")
                    {
                        string[] parsedFN = genomeFN.Split(FNSeparator);
                        string dirName = parsedFN[parsedFN.Length - 2];
                        int organismIdx = directoryNames[dirName];
                        targetGenomeNames[organismIdx] = dirName;
                        if (targetGenomeSequences[organismIdx] == null)
                            targetGenomeSequences[organismIdx] = new List<string>();
                        targetGenomeSizes[organismIdx] = GetSequences(genomeFN, kmerSize, null, targetGenomeSequences[organismIdx], minLength, quiet, verbose);
                    }

                    if (sourceFlag == "-f")
                    {
                        targetGenomeSequences[nextFileIdx] = new List<string>();
                        targetGenomeSizes[nextFileIdx] = GetSequences(genomeFN, kmerSize, null, targetGenomeSequences[nextFileIdx], minLength, quiet, verbose);
                        string lastName = genomeFN.Substring(genomeFN.LastIndexOf('\\') + 1);
                        targetGenomeNames[nextFileIdx] = lastName.Substring(0, lastName.LastIndexOf('.'));
                        nextFileIdx++;
                    }
                }
            }
            else
            // within files
            {
                List<string> geneNames = new List<string>();
                List<string> genes = new List<string>();
                for (int i = 0; i < targetGenomeFNs.Count; i++)
                    GetSequences(targetGenomeFNs[i], kmerSize, geneNames, genes, minLength, quiet, verbose);

                targetGenomeSequences = new List<string>[genes.Count];
                targetGenomeNames = new string[genes.Count];
                targetGenomeSizes = new long[genes.Count];
                noOfTargetGenomes = genes.Count;

                for (int i = 0; i < genes.Count; i++)
                {
                    targetGenomeSequences[i] = new List<string>();
                    targetGenomeSequences[i].Add(genes[i]);
                    targetGenomeNames[i] = geneNames[i].Replace('\t', ' ');
                    targetGenomeSizes[i] = genes[i].Length;
                }
            }

            // parse the sources if they are separate from the targets - otherwise just make the source data structures point to the targets (GxG)
            int noOfSourceGenomes = 0;
            if (sourceFNs.Count != 0)
            {
                noOfSourceGenomes = sourceFNs.Count;
                sourceGenomeSequences = new List<string>[noOfSourceGenomes];
                sourceGenomeNames = new string[noOfSourceGenomes];
                sourceGenomeSizes = new long[noOfSourceGenomes];
                nextFileIdx = 0;
                long maxSourceGenomeSize = 0;

                foreach (string sourceFN in sourceFNs)
                {
                    sourceGenomeSequences[nextFileIdx] = new List<string>();
                    long sourceSize = GetSequences(sourceFN, kmerSize, null, sourceGenomeSequences[nextFileIdx], 0, quiet, verbose);
                    if (sourceSize > maxSourceGenomeSize)
                        maxSourceGenomeSize = sourceSize;
                    string lastName = sourceFN.Substring(sourceFN.LastIndexOf('\\') + 1);
                    sourceGenomeNames[nextFileIdx] = lastName.Substring(0, lastName.LastIndexOf('.'));
                    sourceGenomeSizes[nextFileIdx] = sourceSize;
                    nextFileIdx++;
                }
            }
            else
            {
                sourceGenomeSequences = targetGenomeSequences;
                sourceGenomeNames = targetGenomeNames;
                sourceGenomeSizes = targetGenomeSizes;
                noOfSourceGenomes = noOfTargetGenomes;
            }

            //Console.WriteLine("SL=" + seedLength + "; SM=" + maxSeedMismatches + "; K=" + kmerSize);

            if (betweenFiles)
                Console.WriteLine("Read " + noOfTargetGenomes + " genome files");
            else
                Console.WriteLine("Read " + noOfTargetGenomes + " genes");

            // counts matrix
            long[,] matches = new long[noOfTargetGenomes, noOfTargetGenomes];
            // timings matrix
            long[,] timesForMatches = new long[noOfTargetGenomes, noOfTargetGenomes];

#if PERF
            SharedCode.perfCounts = new long[noThreads][];
            perfTimers = new Stopwatch[noThreads][];
            for (int i = 0; i < noThreads; i++)
            {
                perfCounts[i] = new long[noOfCounts];
                perfTimers[i] = new Stopwatch[noOfPerfTimers];
                for (int t = 0; t < noOfPerfTimers; t++)
                    perfTimers[i][t] = new Stopwatch();
            }
#endif
            if (useSeeds)
            {
                if (slParam != 0)
                {
                    seedLength = slParam;
                    if (smParam != 0) 
                        maxSeedMismatches = smParam;
                    else
                        maxSeedMismatches = (kmerSize - seedLength) / 2;
                }
                seedMask = 0xffffffffffffffff << (64 - seedLength * 2);         
            }

            // test if output file name is safe to use 
            if (countsFN == null)
            {
                Console.WriteLine("no results file specified (-out)");
                return;
            }

            int countsLastDotIdx = countsFN.LastIndexOf('.');
            if (countsLastDotIdx < 0)
            {
                countsPrefix = countsFN;
                countsSuffix = "";
            }
            else
            {
                countsPrefix = countsFN.Substring(0, countsLastDotIdx);
                countsSuffix = countsFN.Substring(countsLastDotIdx);
            }

            StreamWriter counts;
            try
            {
                counts = new StreamWriter(countsPrefix + countsSuffix);
            }
            catch
            {
                Console.WriteLine("Failed to open counts file: " + countsPrefix + countsSuffix);
                return;
            }
            counts.Close();

#if PERF
            stats = new StreamWriter(countsPrefix + "_stats.txt");
            stats.Write("pup (" + release + ") ");
            foreach (string arg in args)
                stats.Write(" " + arg);
            stats.WriteLine();
            stats.WriteLine("\t" + perfHeaders + "\t" + perfTimersHeaders);
#endif
            command += " [k=" + kmerSize + " mm=" + maxMismatches;
            if (exactMatch)
                command += " exact]";
            else
                command += " inexact" + (allowedVariants == VariantsAllowed.subOnly ? " subs" : " all") + " sl=" + seedLength + " sm=" + maxSeedMismatches + "]";
            Console.WriteLine(command);

                // queue matching requests for all of our source genomes/genes against all of the target genome kMer sets
                //   in target order - S1T1, S2T1, S1T2, S2T2
                Queue < MatchingRequest > queuedRequests = new Queue<MatchingRequest>(noOfSourceGenomes * noOfTargetGenomes);
            for (int t = 0; t < noOfTargetGenomes; t++)
            {
                for (int s = 0; s < noOfSourceGenomes; s++)
                {
                    MatchingRequest mr = new MatchingRequest(s, t);
                    queuedRequests.Enqueue(mr);
                }
                MatchingRequest mre = new MatchingRequest(-1, t);
                queuedRequests.Enqueue(mre);
            }

            // cache of previously tiled target genomes (empty until needed, removed after all genomes have been matched against it)
            TiledGenome[] targets = new TiledGenome[noOfTargetGenomes];
            for (int t = 0; t < noOfTargetGenomes; t++)
                targets[t] = new TiledGenome(noOfTargetGenomes);

#if (DEBUG || LLTRACE || HLTRACE)
            noThreads = 1;
#endif
#if (LLTRACE || HLTRACE)                        
            trace.AutoFlush = true;
            trace.WriteLine(command);
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
                matchingParams[b].kmerSize = kmerSize;
                matchingParams[b].queuedRequests = queuedRequests;
                matchingParams[b].sourceGenomeNames = sourceGenomeNames;
                matchingParams[b].sourceGenomeSizes = sourceGenomeSizes;
                matchingParams[b].sourceGenomeSequences = sourceGenomeSequences;
                matchingParams[b].targetGenomeNames = targetGenomeNames;
                matchingParams[b].targetGenomeSizes = targetGenomeSizes;
                matchingParams[b].targetGenomeSequences = targetGenomeSequences;
                matchingParams[b].targets = targets;
                matchingParams[b].matches = matches;
                matchingParams[b].verbose = verbose;
                matchingParams[b].quiet = quiet;
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

            WriteResults(countsPrefix, countsSuffix, sourceGenomeNames, sourceGenomeSizes, targetGenomeNames, targetGenomeSizes, matches, minPIM);
#if PERF
            WriteRectangle(countsPrefix + "_timings_" + countsSuffix, noOfSourceGenomes, sourceGenomeNames, sourceGenomeSizes, noOfTargetGenomes, targetGenomeNames, targetGenomeSizes, timesForMatches, null, 0.0f);
            stats.Close();
#endif

#if (LLTRACE || HLTRACE)
            trace.Close();
#endif

            Console.WriteLine("Finished cross-comparing " + noOfTargetGenomes + " genes/genomes in " + (DateTime.Now - comparisonStart).TotalSeconds.ToString("#.0") + "s");
#if DEBUG   
            if (!useSeeds)
                Console.WriteLine((skipping ? "skip " + scanCheckInterval : "no skip") + ", " + (doubleSubs ? "double" : "single") + " subs");
#endif
        }

        private static void WriteResults(string resultsPrefix, string resultsSuffix, string[] sourceGenomeNames, long[] sourceGenomeSizes, string[] targetGenomeNames, long[] targetGenomeSizes, long[,] matches, float minPIM)
        {
            int noOfSourceGenomes = sourceGenomeNames.Length;
            int noOfTargetGenomes = targetGenomeNames.Length;

            // write out the raw counts matrix (asymmetric)
            WriteRectangle(resultsPrefix + resultsSuffix, noOfSourceGenomes, sourceGenomeNames, sourceGenomeSizes, noOfTargetGenomes, targetGenomeNames, targetGenomeSizes, matches, null, 0.0f);

            // convert the counts to a symmetric (max) percent simularity matrix
            // start by generating an asymmetric PIM
            float[,] PIM = new float[matches.GetLength(0), matches.GetLength(1)];
            for (int i = 0; i < noOfSourceGenomes; i++)
            {
                float sourceSize = (float)sourceGenomeSizes[i];
                for (int j = 0; j < noOfTargetGenomes; j++)
                    PIM[i, j] = (float)matches[i, j] / sourceSize;
            }
            // and now make it symmetrical by taking just the highest of each [i,j]/[j,i] pair
            for (int s = 0; s < noOfSourceGenomes; s++)
                for (int t = 0; t < noOfTargetGenomes; t++)
                {
                    float ijPIM = PIM[s, t];
                    float jiPIM = PIM[t, s];
                    float highestPIM = Math.Max(ijPIM, jiPIM);
                    PIM[s, t] = highestPIM;
                    PIM[t, s] = highestPIM;
                }

            WriteRectangle(resultsPrefix + "_PIM" + resultsSuffix, noOfSourceGenomes, sourceGenomeNames, sourceGenomeSizes, noOfTargetGenomes, targetGenomeNames, targetGenomeSizes, null, PIM, minPIM);

            // and lastly write out the PIM as a sparse matrix (only possible if symmetric GxG)
            if (noOfSourceGenomes == noOfTargetGenomes)
                WriteSparsePIM(resultsPrefix + "_sparsePIM" + resultsSuffix, noOfTargetGenomes, targetGenomeNames, PIM, minPIM);

        }

        // write out the PIM in sparse, inverted form (for usearch cluster_aggd). 100% identity is 0.0, no identity is 1.0. No identity cells are skipped (and assumed to be 1.0)
        private static void WriteSparsePIM(string resultsFN, int noOfGenomes, string[] genomeNames, float[,] PIM, float min)
        {
            StreamWriter spm = new StreamWriter(resultsFN);

            // first write out the diagonal
            for (int i = 0; i < noOfGenomes; i++)
                spm.WriteLine(genomeNames[i] + "\t" + genomeNames[i] + "\t" + (1.0 - PIM[i, i]).ToString("F3"));

            // then the non-zero entries
            for (int i = 0; i < noOfGenomes; i++)
                for (int j = i; j < noOfGenomes; j++)
                {
                    if (i != j && PIM[i, j] >= min && PIM[i, j] != 0.0)
                        spm.WriteLine(genomeNames[i] + "\t" + genomeNames[j] + "\t" + (1.0 - PIM[i, j]).ToString("F3"));
                }

            spm.Close();
        }

        private static void WriteRectangle(string resultsFN, int noOfSourceGenomes, string[] sourceGenomeNames, long[] sourceGenomeSizes, int noOfTargetGenomes, string[] targetGenomeNames, long[] targetGenomeSizes, long[,] matches, float[,] percentages, float minPIM)
        {
            StreamWriter results = new StreamWriter(resultsFN);

            results.Write("\t\t");
            for (int i = 0; i < noOfTargetGenomes; i++)
            {
                results.Write(targetGenomeNames[i]);
                results.Write("\t");
            }
            results.WriteLine();

            results.Write("\t\t");
            for (int i = 0; i < noOfTargetGenomes; i++)
            {
                results.Write(targetGenomeSizes[i]);
                results.Write("\t");
            }
            results.WriteLine();

            for (int s = 0; s < noOfSourceGenomes; s++)
            {
                results.Write(sourceGenomeNames[s]);
                results.Write("\t");
                results.Write(sourceGenomeSizes[s]);
                results.Write("\t");

                for (int t = 0; t < noOfTargetGenomes; t++)
                {
                    if (matches != null)
                        results.Write(matches[s, t]);
                    if (percentages != null)
                    {
                        float percentage = percentages[s, t];
                        if (percentage < minPIM)
                            percentage = 0.0f;
                        results.Write(percentage.ToString("F3"));
                    }
                    results.Write("\t");
                }
                results.WriteLine();
            }

            results.Close();
        }

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Pup [-k kMerSize] [-dirs|-files|-singleFile] [-exact|-inexact|-subs] [-t #threads] [-source sourceFNP] -out countsFN genomeFNs/Dirs " + release);
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
            string[] sourceGenomeNames = theseParams.sourceGenomeNames;
            long[] sourceGenomeSizes = theseParams.sourceGenomeSizes;
            List<string>[] sourceGenomeSequences = theseParams.sourceGenomeSequences;
            string[] targetGenomeNames = theseParams.targetGenomeNames;
            long[] targetGenomeSizes = theseParams.targetGenomeSizes;
            TiledGenome[] targets = theseParams.targets;
            List<string>[] targetGenomeSequences = theseParams.targetGenomeSequences;
            long[,] matches = theseParams.matches;
            int minLength = theseParams.minLength;
            long[,] timesForMatches = theseParams.timesForMatches;
            bool verbose = theseParams.verbose;
            bool quiet = theseParams.quiet;

            bool finishThread = false;
            MatchingRequest matchingRequest = null;

            List<VariantMerSet> kmerSetPool = new List<VariantMerSet>();
            VariantMerSet outerMerSet = new VariantMerSet(largeVariantSetSize);

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

                    // breakpoint
                    //if (s == 0 && t == 1)
                    //    Debugger.Break();

                    //Console.WriteLine("starting match of " + genomeNames[s] + " to " + genomeNames[t] + " on #" + threadNo);
                    //trace.WriteLine("starting match of " + genomeNames[s] + " to " + genomeNames[t] + " on #" + threadNo);

                    long targetSize = PrepTargetMerSet(targets, t, null, targetGenomeNames[t], kmerSize, targetGenomeSequences[t], minLength, verbose, quiet);
                    if (targetSize > 0)
                        targetGenomeSizes[t] = targetSize;
                    List<string> sourceGenome = sourceGenomeSequences[s];

                    DateTime matchingStart = DateTime.Now;
                    int timeTaken;

                    matches[s, t] = MatchSourceToTarget(exactMatch, allowedVariants, kmerSize, sourceGenome, targets[t], s, t, targetGenomeNames[s], targetGenomeNames[t], kmerSetPool, outerMerSet, threadNo, out timeTaken, quiet, verbose);
                    timesForMatches[s, t] = timeTaken;

                    if (verbose)
                        Console.WriteLine("Matched " + matches[s, t] + "/" + sourceGenomeSizes[s] + " " + kmerSize + "-mers (" + sourceGenomeSequences[s].Count + " contigs) from " + sourceGenomeNames[s] + " to " + targetGenomeNames[t] + " in " + (DateTime.Now - matchingStart).TotalSeconds.ToString("#.0") + "s");

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
                    ReleaseTargetMers(t, targetGenomeNames[t], targets);
                }
            }

            //Console.WriteLine("finished matching thread " + threadNo);
        }

        private static long MatchSourceToTarget(bool exactMatch, VariantsAllowed allowedVariants, int kmerSize, List<string> sourceGenome, TiledGenome target,
                                               int sourceIdx, int targetIdx, string sourceName, string targetName,
                                               List<VariantMerSet> kmerSetPool, VariantMerSet outerMerSet, int threadNo, out int timeTaken, bool quiet, bool verbose)
        {
            long matchingMers = 0;
            Stopwatch sw = new Stopwatch();
            sw.Start();

            //Console.WriteLine("matching " + sourceName + " to " + targetName);

            // tile and match each sequence in the source genome/gene against the current target (re-tiling every time to support fix-on-fuzzy-match)
            for (int c = 0; c < sourceGenome.Count; c++)
            {
                Sequence sequence = new Sequence(sourceGenome[c]);
                if (sequence.Length < kmerSize)
                    continue;
#if (LLTRACE || HLTRACE )
                trace.WriteLine("matching " + sourceName + "[" + c + "] to " + targetName);
#endif
                int mersMatchedForSeq = MatchSequenceToTarget(exactMatch, allowedVariants, ScanType.forwardScan, kmerSize, sequence, 0, target, Strand.either, sourceIdx, targetIdx, c, kmerSetPool, outerMerSet, threadNo);
                matchingMers += mersMatchedForSeq;
            } // for each source sequence

            sw.Stop();
            timeTaken = (int)sw.ElapsedMilliseconds;
            if (!quiet)
                Console.WriteLine("matched " + sourceName + " to " + targetName);
            return matchingMers;
        }
    }
    public class MatchingThreadParams
    {
        public int threadNo;
        public bool exactMatch;
        public VariantsAllowed allowedVariants;
        public int kmerSize;
        public Queue<MatchingRequest> queuedRequests;
        public string[] sourceGenomeNames;
        public long[] sourceGenomeSizes;
        public List<string>[] sourceGenomeSequences;
        public string[] targetGenomeNames;
        public long[] targetGenomeSizes;
        public List<string>[] targetGenomeSequences;
        public TiledGenome[] targets;
        public int minLength;
        public long[,] matches;
        public bool verbose;
        public bool quiet;
        public long[,] timesForMatches;
    }
}


