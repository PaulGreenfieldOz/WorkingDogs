using System;
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Diagnostics;
using WorkingDogsCore;
using PupClancy;
using static PupClancy.SharedCode;


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
    // The comparisons can be exact or single-base variants (inexact/fuzzy), and either just subs or subs/dels/ins.
    //
    // The result is a table with each row being a called gene/contig followed by how many matches it got from each of the genomes
    //
    // usage: Clancy [-k kmerSize] [-exact|-inexact|-subs] [-t #threads] [-s] [-genes|-genomes] -out countsFN sourceFN targetFNs
    //
    // NM-KH31_C3.ffn Priestia_aryabhattai_genomes\*.fna Priestia_megaterium_genomes\*.fna -k 25 -s -subs -t 32 -out Clancy_counts_PaPm_subs.txt
    // NM-KH31_C3.ffn Close_relatives\*.fna -k 25 -s -subs -t 32 -out Clancy_counts_CR_subs.txt
    // NM-KH31_C3.ffn "Priestia plasmids.fna" -genes -k 25 -s -subs -t 32 -out Clancy_counts_plasmid_subs.txt
    // clancy -inexact -t 16 Methanoculleus_bourgensis_GCF_000304355.2_Mb_MS2.ffn -targets genomes *.fna -o Mb_MS2_onto_allMc_gxG_PC142.txt
    // clancy -t 24 Aspergillus_flavus_annotations\ncbi_dataset\data\GCF_009017415.1\cds_from_genomic.ffn -targets genomes -inexact -s Aspergillus_flavus_genomes/*.fna Aspergillus_oryzae_genomes/*.fna -out A.flavus_NRRL_3357_to_A.flavus_oryzae_gxG.txt
    // clancy -t 24 "B:\Pup\Aspergillus\Aspergillus_flavus_annotations\ncbi_dataset\data\GCF_009017415.1\cds_from_genomic.ffn" -targets genomes -inexact -s *.fna -out A.flavus_NRRL_3357_to_A.sample_gxG.txt
    class Program
    {

        public const string release = "v1.4.2";

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                WriteUsage();
                return;
            }

            string command = "Clancy " + release;
            for (int p = 0; p < args.Length; p++)
                command += " " + args[p];
            SearchOption so = SearchOption.TopDirectoryOnly;
            bool exactMatch = false;
            VariantsAllowed allowedVariants = VariantsAllowed.subInsDel;
            int kmerSize = 21;
            int noThreads = Environment.ProcessorCount / 2;
            List<string> FNParams = new List<string>();
            string? countsFN = null;
            string? countsPrefix = null;
            string countsSuffix = ".txt";
            // each target file is treated as a single target, or the single target file is treated as a file of target seqs (genes/contigs)
            bool mapToGenome = true;
            bool verbose = false;
            bool quiet = false;
            // min length for target - useful when matching to metagenome
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
                            Console.WriteLine("kMer length must be 0>k<=31");
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
                        string targetsArg = args[p + 1].ToLower();
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

                    if (args[p] == "-out" || args[p] == "-output" || args[p] == "-o")
                    {
                        if (!CheckForParamValue(p, args.Length, "file name expected after -out"))
                            return;

                        countsFN = args[p + 1];
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

                FNParams.Add(args[p]);
            }

            smallVariantSetSize = 3 + 1 + 4;
            largeVariantSetSize = 3 * kmerSize + kmerSize + 4 * kmerSize;

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

#if (LLTRACE || HLTRACE)
            StreamWriter trace = new StreamWriter("trace.txt");
            SharedCode.trace = trace;
#endif

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
                    GetSequences(FN, kmerSize, th, ts, minLength, quiet, verbose);
                    foreach (string h in th)
                        targetHeaders.Add(h);
                    foreach (string s in ts)
                        targetSeqs.Add(s);
                }
                noOfTargets = targetHeaders.Count;
            }

            DateTime comparisonStart = DateTime.Now;

            char[] FNSeparator = new char[] { Path.DirectorySeparatorChar };

            //Console.WriteLine("SL=" + seedLength + "; SM=" + maxSeedMismatches + "; K=" + kmerSize);

            // get all the source (gene/contig) seqs
            List<string> geneHeaders = new List<string>(50000);
            List<string> geneSeqs = new List<string>(50000);
            GetSequences(sourceFN, kmerSize, geneHeaders, geneSeqs, 0, quiet, verbose);

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

            // cache of previously tiled target genomes (empty until needed, removed after all genomes have been matched against it)
            TiledGenome[] targets = new TiledGenome[noOfTargets];
            for (int t = 0; t < noOfTargets; t++)
                targets[t] = new TiledGenome(geneSeqs.Count);

            Queue<List<ulong>> seedListPool = new Queue<List<ulong>>(100000);

#if (DEBUG || LLTRACE || HLTRACE)
            noThreads = 1;
#endif
#if (LLTRACE || HLTRACE)                        
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
                matchingParams[b].kmerSize = kmerSize;
                matchingParams[b].queuedRequests = queuedRequests;
                matchingParams[b].sourceNames = geneHeaders;
                matchingParams[b].sourceSeqs = geneSeqs;
                matchingParams[b].targetFNs = targetFNs;
                matchingParams[b].targetSeqs = targetSeqs;
                matchingParams[b].targetNames = targetNames;
                matchingParams[b].targetSizes = targetSizes;
                matchingParams[b].targets = targets;
                matchingParams[b].minLength = minLength;
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

            WriteResults(countsFN, geneHeaders, geneSeqs, kmerSize, targetNames, targetSizes, matches);

#if (LLTRACE || HLTRACE)
            trace.Close();
#endif

            Console.WriteLine("Finished comparing " + geneSeqs.Count + " genes to " + noOfTargets + " genes/genomes in " + (DateTime.Now - comparisonStart).TotalSeconds.ToString("#.0") + "s");
            //if (!useSeeds)
            //    Console.WriteLine((skipping ? "skip " + scanCheckInterval : "no skip") + ", " + (doubleSubs ? "double" : "single") + " subs");
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
                results.Write(seqs[s].Length - kMerSize + 1);
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
            Console.WriteLine("usage: Clancy [-k kmerSize] [-exact|-inexact|-subs] [-t #threads] [-targets genes|genomes] -out countsFN [-s] sourceFN targetFNs " + release);
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
            List<string> sourceNames = theseParams.sourceNames;
            List<string> sourceSeqs = theseParams.sourceSeqs;
            List<string> targetFNs = theseParams.targetFNs;
            List<string> targetSeqs = theseParams.targetSeqs;
            string[] targetNames = theseParams.targetNames;
            long[] targetSizes = theseParams.targetSizes;
            TiledGenome[] targets = theseParams.targets;
            int minLength = theseParams.minLength;
            long[,] matches = theseParams.matches;
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

                    //Console.WriteLine("starting match of " + targetNames[s] + " to " + targetNames[t] + " on #" + threadNo);
                    //trace.WriteLine("starting match of " + targetNames[s] + " to " + targetNames[t] + " on #" + threadNo);

                    List<string> genomeSequences = new List<string>(10000);

                    // -genes: target sequences already known so just copy it over
                    if (targetSeqs.Count > 0)
                        genomeSequences.Add(targetSeqs[t]);

                    long targetSize = PrepTargetMerSet(targets, t, targetFNs[t], targetNames[t], kmerSize, genomeSequences, minLength, verbose, quiet);
                    if (targetSize > 0)
                        targetSizes[t] = targetSize;
                    string sourceSeq = sourceSeqs[s];

                    DateTime matchingStart = DateTime.Now;
                    int timeTaken;

#if HLTRACE
                    trace.WriteLine("starting match [" + s + "," +  t + "]: " + sourceNames[s] + " to " + targetNames[t]);
#endif
                    long StoTMatches = MatchSourceToTarget(exactMatch, allowedVariants, kmerSize, sourceSeqs[s], targets[t], s, t, sourceNames[s], targetNames[t], kmerSetPool, outerMerSet, threadNo, out timeTaken);

                    //if (matches[s, t] != 0 && matches[s, t] != StoTMatches)
                    //    Console.WriteLine("[" + s + "," + t + "]: " + matches[s, t] + "-->" + StoTMatches);
                    matches[s, t] = StoTMatches;
                    timesForMatches[s, t] = timeTaken;
                    //Console.WriteLine("matched " + s + ", " + t);

                    ReleaseTargetMers(t, targetNames[t], targets);
                }
            }

            //Console.WriteLine("finished matching thread " + threadNo);
        }
    }
    public class MatchingThreadParams
    {
        public int threadNo;
        public bool exactMatch;
        public VariantsAllowed allowedVariants;
        public int kmerSize;
        public Queue<MatchingRequest> queuedRequests;
        public List<string> sourceNames;
        public List<string> sourceSeqs;
        public List<string> targetFNs;
        public List<string> targetSeqs;
        public string[] targetNames;
        public long[] targetSizes;
        public TiledGenome[] targets;
        public int minLength;
        public long[,] matches;
        public bool verbose;
        public bool quiet;
        public long[,] timesForMatches;
    }
}


