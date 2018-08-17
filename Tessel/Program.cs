using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Runtime;
using WorkingDogsCore;
using MerCollections;

namespace Tessel
{
    public class Program
    {
        // Tiles a set of read files into canonical k-mers and writes out a file containing these k-mers and how many times each of them appears.
        //
        // usage: Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-tmp tempDir] [-min minCount] [-tq] [-s] [pirs] [-nopairs] [-canonical|asread] [-text textFN] [-textFormat pairs|sum|faPairs|faSum] cbtName readsFN... or readsPattern
        //
        // Tessel  -k 25 -g 6000000 -t 4 CsporParamTest s_1_sequence_merged.fastq
        // -k 25 -g 60000000 -t 16 -min 1  -tmp g:\temp\ TesselTest-16 s_1_?_sequence.fastq

        const string version = "2.1.0";                         // version - update on every significant change

        static Stopwatch timeSpentFilling = new Stopwatch();
        static TimeSpan timeSpentSorting = TimeSpan.Zero;
        static TimeSpan timeSpentWriting = TimeSpan.Zero;

        // progress monitors
        const int reportInterval = 60000;           // report back at this interval (ms)
        static bool stopMonitor = false;            // tell monitor thread to finish

        // current statistics
        enum phases
        {
            tilingMers,
            writingMers,
            tilingPairs,
            writingPairs
        }
        static phases programPhase = phases.tilingMers;     // for rate reporter
        static long progressReadsRead = 0;                  // no. of reads
        static long progressMersToWrite = 0;
        static long progressMersWritten = 0;
        static long progressPairsTiled = 0;
        static long progressPairsWritten = 0;
        static int kLength = 0;                             // k-mer length - only for reporting

        const int maxReadSize = 1000;
        const int batchSize = 1000;
        const int targetBufferSize = 150000;
        const int defaultReadLength = 300;
        const int defaultHeaderLength = 100;
        const int readBufferSize = 1000000;

        // text output formats (driven by options - default is binary only)
        enum TextFormat
        {
            tabPairs,
            tabSummed,
            faPairs,
            faSummed
        }; 

        public enum TypeOfCount
        {
            ulongType,          // for .cbt files. Values are a pair of counts packed into a ulong
            uintType            // for .prs file. Values are a single count (uint)
        }

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                // Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-pairs|-nopairs] [-trim nn] [-trimQual nn] [-tmp tempDir] [-min minCount] [-s] cbtName readsFN... or readsPattern
                WriteUsage();
                return;
            }

            int readsFormat = SeqFiles.formatNone; // format e.g. fasta, fasta, sfa
            string cbtName = null;               // output binary file
            bool textMers = false;               // write out k-mers to a text file
            string textFileName = null;          // and the name of the text file
            TextFormat wantedtextFormat = TextFormat.tabPairs;  // and the format wanted

            int merSize = 0;                     // tiling length
            long genomeSize = 0;                 // estimated genome size
            int trimLength = 0;                  // how much to trim from end of reads before tiling
            int trimQual = 0;                    // min qual value for end of read qual trimming
            int qualOffset = 0;                  // convert from qual character to score by subtracting this value
            int noThreads = 1;                   // no. of parallel threads used in counting and writing
            int minCount = 1;                    // only write k-mers with counts >= to .cbt file
            string tempDir = "";                 // temporary directory for saving flushed singletons
            bool flushSingletons = true;         // (normal) flush full singletons tables to files
            bool inexactCount = false;           // generate rough count for HiFreq kMers (ala Jellyfish) by discarding singleton tables
            bool canonical = true;               // save canonical Kmers - or as they appear in reads/genome
            System.IO.SearchOption readsSearchOption = SearchOption.TopDirectoryOnly;  // default is top-directory only; -s option forces recursive search
            bool generatePairs = true;           // generate kMer pairs file (.prs) after kMers

            string myProcessNameAndArgs = null;  // extract the command arguments from the program
            List<string> FNParams = new List<string>();    // the set of file names or patterns to be (expanded and) tiled. The first will be the .cbt name.
            string fnSeparator = Path.DirectorySeparatorChar.ToString();  // \ for Windows; / for Unix/Linux

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-k")
                    {
                        if (!CheckForParamValue(p, args.Length, "k-mer length number expected after -k"))
                            return;
                        try
                        {
                            merSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        if (merSize > 32)
                        {
                            Console.WriteLine("k-mer length must be <= 32");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-trim")
                    {
                        if (!CheckForParamValue(p, args.Length, "no. of bases to trim expected after -trim"))
                            return;
                        try
                        {
                            trimLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number after the -trim parameter: " + args[p + 1]);
                            return;
                        }
                        if (trimLength < 0)
                        {
                            Console.WriteLine("trim length must be >= 0");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-trimqual" || args[p] == "-tq")
                    {
                        if (!CheckForParamValue(p, args.Length, "min qual expected after -trimqual"))
                            return;
                        try
                        {
                            trimQual = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number after the -trimqual parameter: " + args[p + 1]);
                            return;
                        }
                        if (trimQual < 0 || trimQual > 63)
                        {
                            Console.WriteLine("min qual must be between 0 and 63");
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-g" || args[p] == "-genome")
                    {
                        if (!CheckForParamValue(p, args.Length, "genome length number expected after -g|-genome"))
                            return;
                        try
                        {
                            genomeSize = Convert.ToInt64(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -g|-genome parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
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

                    if (args[p] == "-t" || args[p] == "-threads")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -t|-threads"))
                            return;
                        try
                        {
                            noThreads = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -t|-threads parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-m" || args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -m|-min"))
                            return;
                        try
                        {
                            minCount = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -m|-min parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-tmp")
                    {
                        if (!CheckForParamValue(p, args.Length, "temp directory name expected after -tmp"))
                            return;
                        tempDir = args[p + 1];
                        p++;
                        continue;
                    }

                    if (args[p] == "-s")
                    {
                        readsSearchOption = SearchOption.AllDirectories;
                        continue;
                    }

                    if (args[p] == "-pairs")
                    {
                        generatePairs = true;
                        continue;
                    }

                    if (args[p] == "-nopairs")
                    {
                        generatePairs = false;
                        continue;
                    }

                    if (args[p] == "-noflush")
                    {
                        flushSingletons = false;
                        continue;
                    }

                    if (args[p] == "-inexact")
                    {
                        inexactCount = true;
                        continue;
                    }

                    if (args[p] == "-canonical")
                    {
                        canonical = true;
                        continue;
                    }

                    if (args[p] == "-asread")
                    {
                        canonical = false;
                        continue;
                    }

                    if (args[p] == "-text")
                    {
                        if (!CheckForParamValue(p, args.Length, "file name expected after -text"))
                            return;
                        textMers = true;
                        textFileName = args[p + 1];
                        wantedtextFormat = TextFormat.tabPairs;
                        p++;
                        continue;
                    }

                    if (args[p] == "-textformat" || args[p] == "-tf")
                    {
                        if (!textMers)
                        {
                            Console.Write("-text option must be set before specifying text format");
                            return;
                        }
                        if (!CheckForParamValue(p, args.Length, "file format expected after -textformat"))
                            return;
                        string textFormat = args[p + 1].ToLower();
                        switch (textFormat)
                        {
                            case "pairs":
                                wantedtextFormat = TextFormat.tabPairs;
                                break;
                            case "sum":
                                wantedtextFormat = TextFormat.tabSummed;
                                break;
                            case "fapairs":
                                wantedtextFormat = TextFormat.faPairs;
                                break;
                            case "fasum":
                                wantedtextFormat = TextFormat.faSummed;
                                break;
                            default:
                                Console.WriteLine("unknown text format: " + textFormat);
                                break;
                        }
                        p++;
                        continue;
                    }

                    Console.WriteLine("unrecognised option: " + args[p]);
                    WriteUsage();
                    return;
                }

                FNParams.Add(args[p]);
            }


            // track what program and args produced the result files
            Process myProcess = Process.GetCurrentProcess();
            myProcessNameAndArgs = myProcess.ProcessName;
            foreach (string a in args)
                myProcessNameAndArgs = myProcessNameAndArgs + " " + a;

            // validate the temp directory
            if (tempDir != "")
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

            if (FNParams.Count < 2)
            {
                Console.WriteLine("Expected a tile-file name and at least one reads file name or pattern");
                return;
            }

            kLength = merSize;

            string seqFilesDir = Directory.GetCurrentDirectory();
            string singletonsDir = tempDir == "" ? (seqFilesDir + fnSeparator) : tempDir;

            // the .cbt name is expected to be the first of the non-option-like parameters
            string cbtSuffix = ".cbt";
            if (!canonical)
                cbtSuffix = ".abt";
            cbtName = FNParams[0];
            FNParams.RemoveAt(0);
            // and remove the .cbt suffix if it was present
            if (cbtName.ToLower().EndsWith(cbtSuffix))
                cbtName = cbtName.Substring(0, cbtName.Length - cbtSuffix.Length);

            // find out if we're to write the .cbt to a directory - and make sure the directory is present
            if (cbtName.Contains(fnSeparator))
            {
                string cbtDir = cbtName.Substring(0, cbtName.LastIndexOf(fnSeparator));
                if (!Directory.Exists(cbtDir))
                    Directory.CreateDirectory(cbtDir);
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
                string[] matchedReadsFNs = Directory.GetFiles(readsFilePaths[f], readsFileNames[f], readsSearchOption);
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

            if (readsFNs.Length == 0)
            {
                Console.WriteLine("No matching read files found");
                return;
            }

            // if we weren't told what format the reads are in, use the format from the first file 
            if (readsFormat == SeqFiles.formatNone)
            {
                StreamReader formatTester = new StreamReader(readsFNs[0]);
                string firstLine = formatTester.ReadLine();
                if (firstLine[0] == '>')
                    readsFormat = SeqFiles.formatFNA;
                if (firstLine[0] == '@')
                    readsFormat = SeqFiles.formatFASTQ;
                formatTester.Close();
            }

            if (trimQual > 0 && readsFormat != SeqFiles.formatFASTQ)
            {
                Console.WriteLine("qual trimming only supported with fastq files");
                return;
            }

            bool fullQualHeader;
            if (trimQual > 0)
                qualOffset = SeqFiles.ResolveFastqQualAmbiguity(readsFNs[0], out fullQualHeader);

            // start the monitor/synchronising thread
            programPhase = phases.tilingMers;
            Thread monitorProgress = new Thread(RateReporter);
            monitorProgress.Priority = ThreadPriority.AboveNormal;
            monitorProgress.Start();

            Stopwatch countingTimer = new Stopwatch();
            Stopwatch mergingTimer = new Stopwatch();
            Stopwatch sortingTimer = new Stopwatch();
            Stopwatch pairsTimer = new Stopwatch();
            Stopwatch totalTimer = new Stopwatch();

            totalTimer.Start();

            // allocate the overall stats (updated here when threads finish)
            ThreadStats tesselStats = new ThreadStats();

            // allocate the kMer tables
            MerTables merTables = new MerTables(genomeSize, merSize, singletonsDir, noThreads, flushSingletons, inexactCount, canonical);

            //Console.WriteLine("memory at start of counting: " + GC.GetTotalMemory(false));

            // Tiling and Counting phase   
            // =========================
            countingTimer.Start();

            for (int f = 0; f < readsFNs.Length; f++)
            {
                TileAndCount(readsFNs, f, readsFormat, merSize, trimLength, trimQual, qualOffset, noThreads, merTables, tesselStats);
            }

            countingTimer.Stop();

            //Console.WriteLine("memory at end of counting: " + GC.GetTotalMemory(false));

            //Console.WriteLine(countingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts counting");
            //for (int i = 0; i < merTables.repeatedMers.Length; i++)
            //{
            //    Dictionary<int, int> repeatedMersStats = merTables.repeatedMers[i].CheckHashTable();
            //    int[] lengths = new int[repeatedMersStats.Count];
            //    int[] counts = new int[repeatedMersStats.Count];
            //    repeatedMersStats.Keys.CopyTo(lengths, 0);
            //    repeatedMersStats.Values.CopyTo(counts, 0);
            //    Array.Sort(lengths, counts);
            //    for (int j = 0; j < repeatedMersStats.Count; j++)
            //        Console.WriteLine(i + " " + lengths[j] + " " + counts[j]);
            //}

            // sorting in-memory tables prior to merge (flushed tables already sorted)
            // =======================================================================
            sortingTimer.Start();

            SortMerTables(noThreads, merTables);

            sortingTimer.Stop();           
            //Console.WriteLine("memory at end of sorting: " + GC.GetTotalMemory(false));

            // wait for any pending singleton writes to complete (overlapped with sorting to give the writes more time)
            merTables.Stabilise();

            // merging tables and writing distinct kMers and counts (in .cbt format)
            // =====================================================================
            mergingTimer.Start();
            programPhase = phases.writingMers;

            string cbtFN = cbtName + "_" + merSize + cbtSuffix;
            Dictionary<int, long> sumReps = new Dictionary<int, long>(50000);        // for generating histogram of sum counts
            progressMersToWrite = tesselStats.noOfTiledMers;

            //Console.WriteLine(GC.GetTotalMemory(false) + " memory at end of counting");
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //    Console.WriteLine("merTable[" + p + "].length = " + merTables.repeatedMers[p].lengthEntries + "\tcount= " + merTables.repeatedMers[p].Count);

            //histo.WriteLine(GC.GetTotalMemory(false) + " memory at end of counting");
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //    histo.WriteLine("merTable[" + p + "].length = " + merTables.repeatedMers[p].length + "\t" + merTables.repeatedMers[p].Count);

            MergeAndWrite(noThreads, cbtFN, TypeOfCount.ulongType, minCount, merSize, merTables.repeatedMers, merTables.overflowMers, merTables.singletonFilters, merTables.singletonsWaitingFlush,
                            merTables.flushedSingletonFNs, merTables.firstFlushedSingletonMer, merTables.lastFlushedSingletonMer, merTables.flushedSingletonsCount,
                            merTables.flushedLowRepsFNs, merTables.firstFlushedLowRepsMer, merTables.lastFlushedLowRepsMer, merTables.flushedLowRepsCount, tesselStats, sumReps);

            mergingTimer.Stop();

            //Console.WriteLine("memory at end of merging/writing: " + GC.GetTotalMemory(false));

            // Dump expanded k-mers and counts to a text file
            // ==============================================
            if (textMers)
            {
                Console.WriteLine("Writing k-mers to text file...");
                StreamWriter results = new StreamWriter(File.Open(textFileName, FileMode.Create, FileAccess.Write));
                BinaryReader cbtReader = new BinaryReader(new FileStream(cbtFN, FileMode.Open, FileAccess.Read, FileShare.Read, 100000, FileOptions.SequentialScan));
                WriteExpandedMers(cbtReader, results, wantedtextFormat);
                results.Close();
                cbtReader.Close();
            }

            // generate kMers pairs if this option is set
            if (generatePairs)
            {
                programPhase = phases.tilingPairs;

                // clean up memory used by counting/merging/writing
                merTables = null;
                GC.Collect();

                // tables used for pairs generation
                kMerTable kMersTable;                   // the collection of kMers read back from the .cbt file
                MerCollections.PairTable distinctPairs; // shared pairs dictionary (ulong, long) - thread safe as there are no deletions or resizes

                // pairs constants 
                const int endGuard = 16;                // don't get too close to the error-prone end
                const int minGap = 16;                  // don't let the gap get too small
                const int maxGap = 64;                  // and don't let the gap be too big
                int pairGap = 0;                        // (constant) gap between mers
                int minReadLength = 0;                  // minimum length for a read to be tiled for pairs

                string pairsFN = cbtFN.Replace(".cbt", ".prs");

                // look at the first file to determine a likely read length (for gap calculations)
                StreamReader testReader = new StreamReader(readsFNs[0]);
                int readLength = 0;
                for (int i = 0; i < 20; i++)
                {
                    string nextRead = SeqFiles.ReadRead(testReader, readsFormat);
                    if (nextRead == null)
                        break;
                    int nextLength = nextRead.Length;
                    if (nextLength > readLength)
                        readLength = nextLength;
                }
                testReader.Close();

                // load the .cbt file into a kMerTable
                kMersTable = new kMerTable(cbtFN, minCount, 1000);
                int averageDepth = kMersTable.averageDepthLoaded;
                long loadedDistinctMers = kMersTable.distinctMersLoaded;

                if (merSize < kMerPairs.pairFragmentSize)
                {
                    Console.WriteLine("mers in .cbt file are shorter than pair fragment size: " + merSize + " < " + kMerPairs.pairFragmentSize);
                    return;
                }

                // have to able to fit at least two full mers into the read (no overlaps)
                minReadLength = 2 * merSize;
                if (readLength < minReadLength)
                {
                    Console.WriteLine("reads too short to generate pairs: " + readLength);
                    return;
                }

                distinctPairs = new MerCollections.PairTable(loadedDistinctMers, noThreads, 2 * kMerPairs.pairFragmentSize);

                // calculate a gap size based on the first read
                pairGap = (readLength - endGuard) / 2 - (kMerPairs.pairFragmentSize * 2);
                if (pairGap < minGap)
                    pairGap = minGap;
                if (pairGap > maxGap)
                    pairGap = maxGap;

                // start read counter again
                progressReadsRead = 0;

                foreach (string readsFN in readsFNs)
                    CountPairsInReadsFile(readsFN, readsFormat, noThreads, merSize, minReadLength, pairGap, minCount, kMersTable, distinctPairs, tesselStats);

                progressPairsTiled = tesselStats.noOfTiledPairs;

                SortPairsTables(distinctPairs, noThreads, tesselStats);

                programPhase = phases.writingPairs;

                MergeAndWrite(noThreads, pairsFN, TypeOfCount.uintType, minCount, pairGap, distinctPairs.repeatedMers, distinctPairs.overflowMers,
                              null, null, null, null, null, null,null, null, null, null, tesselStats, null);

            }

            // write kMer histogram and Tessel stats
            // =====================================
            string histoFN = cbtName + "_" + merSize + "_histo.txt";

            WriteHistogram(histoFN, sumReps, myProcessNameAndArgs, countingTimer, mergingTimer, sortingTimer, tesselStats);

            StopMonitorThread(monitorProgress);

            totalTimer.Stop();
            Console.WriteLine("kMer/pairs tiling took " + totalTimer.Elapsed.TotalSeconds.ToString("#.0") + "s");
        }

        private static void SortPairsTables(MerCollections.PairTable pairsTable, int noThreads, ThreadStats tesselStats)
        {
            // convert all the in-memory structures from hash tables to sorted arrays 
            int noOfSortingTasks = pairsTable.repeatedMers.Length + pairsTable.overflowMers.Length;   // a good enough approximation
            Queue<SorterTask> sortingTasks = new Queue<SorterTask>(noOfSortingTasks);

            // prepare the sorting tasks
            for (int i = 0; i < pairsTable.repeatedMers.Length; i++)
            {
                SorterTask st = new SorterTask();
                st.merDictionary = pairsTable.repeatedMers[i];
                sortingTasks.Enqueue(st);
            }
            for (int i = 0; i < pairsTable.overflowMers.Length; i++)
                if (pairsTable.overflowMers[i] != null)
                {
                    SorterTask st = new SorterTask();
                    st.merDictionary = pairsTable.overflowMers[i];
                    sortingTasks.Enqueue(st);
                }

            int noSortingThreads = 0;
            if (sortingTasks.Count > noThreads)
                noSortingThreads = noThreads;
            else
                noSortingThreads = sortingTasks.Count;

            Thread[] sortingThreads = new Thread[noSortingThreads];
            SorterParams sorterParams = new SorterParams();

            for (int t = 0; t < noSortingThreads; t++)
            {
                sorterParams.waitingTasks = sortingTasks;
                sortingThreads[t] = new Thread(new ParameterizedThreadStart(Program.PairSortingThread));
                sortingThreads[t].Priority = ThreadPriority.BelowNormal;
                sortingThreads[t].Start(sorterParams);
            }

            for (int t = 0; t < noSortingThreads; t++)
            {
                sortingThreads[t].Join();
                sortingThreads[t] = null;
                //Console.WriteLine("finished sorting thread " + (t + 1));
            }
        }

        private static void CountPairsInReadsFile(string readsFN, int readsFormat, int noThreads, int merSize, int minReadLength, int pairGap, 
                                                  int minCount, kMerTable kMersTable, MerCollections.PairTable pairsTable, ThreadStats tesselStats)
        {
            Console.WriteLine("Generating pairs from " + readsFN);
            StreamReader reads = new StreamReader(readsFN, Encoding.ASCII, false, 1000000);
            BufferedReader bufferedReads = new BufferedReader(readsFormat, reads, null);

            PairCountingThreadParams[] countingParams = new PairCountingThreadParams[noThreads];
            Thread[] countingThreads = new Thread[noThreads];
            for (int t = 0; t < noThreads; t++)
            {
                countingParams[t] = new PairCountingThreadParams();
                countingParams[t].threadNumber = t;
                countingParams[t].readsFile = bufferedReads;
                countingParams[t].merSize = merSize;
                countingParams[t].minReadLength = minReadLength;
                countingParams[t].pairGap = pairGap;
                countingParams[t].minCount = minCount;
                countingParams[t].kMersTable = kMersTable;
                countingParams[t].pairsTable = pairsTable;
                countingParams[t].threadStats = new ThreadStats();
                countingThreads[t] = new Thread(new ParameterizedThreadStart(PairCountingThread));
                countingThreads[t].Priority = ThreadPriority.BelowNormal;
                countingThreads[t].Start(countingParams[t]);
            }

            //  and wait for them all to finish
            for (int t = 0; t < noThreads; t++)
            {
                countingThreads[t].Join();

                tesselStats.noOfTiledPairs += countingParams[t].threadStats.noOfTiledPairs;

                countingThreads[t] = null;
                //Console.WriteLine("finished pair counting thread " + (t + 1));
            }

            bufferedReads.Close();
        }

        private static void SortMerTables(int noThreads, MerTables merTables)
        {
            // convert all the in-memory structures from hash tables to sorted arrays 
            int noOfSortingTasks = merTables.repeatedMers.Length + merTables.singletonFilters.Length + merTables.repeatedMers.Length;   // a good enough approximation
            Queue<SorterTask> sortingTasks = new Queue<SorterTask>(noOfSortingTasks);

            // prepare the sorting tasks
            for (int i = 0; i < merTables.singletonFilters.Length; i++)
            {
                SorterTask ufs = new SorterTask();
                ufs.merCollection = merTables.singletonFilters[i];
                ufs.merDictionary = null;
                sortingTasks.Enqueue(ufs);
                //Console.WriteLine("queued singletons[" + i + "] (" + merTables.singletonFilters[i].Count + ") for sorting");
            }
            for (int i = 0; i < merTables.singletonsWaitingFlush.Length; i++)
            {
                if (merTables.singletonsWaitingFlush[i] != null)
                    foreach (SingletonCollection waitingSingletons in merTables.singletonsWaitingFlush[i])
                    {
                        SorterTask ufs = new SorterTask();
                        ufs.merCollection = waitingSingletons;
                        ufs.merDictionary = null;
                        sortingTasks.Enqueue(ufs);
                        //Console.WriteLine("queued unflushed waiting singletons[" + i + "] (" + waitingSingletons.Count + ") for sorting");
                    }
            }
            for (int i = 0; i < merTables.repeatedMers.Length; i++)
            {
                SorterTask rmd = new SorterTask();
                rmd.merCollection = null;
                rmd.merDictionary = merTables.repeatedMers[i];
                sortingTasks.Enqueue(rmd);
            }
            for (int i = 0; i < merTables.overflowMers.Length; i++)
                if (merTables.overflowMers[i] != null)
                {
                    SorterTask rmo = new SorterTask();
                    rmo.merCollection = null;
                    rmo.merDictionary = merTables.overflowMers[i];
                    sortingTasks.Enqueue(rmo);
                }

            int noSortingThreads = 0;
            if (sortingTasks.Count > noThreads)
                noSortingThreads = noThreads;
            else
                noSortingThreads = sortingTasks.Count;

            Thread[] sortingThreads = new Thread[noSortingThreads];
            SorterParams sorterParams = new SorterParams();

            for (int t = 0; t < noSortingThreads; t++)
            {
                sorterParams.waitingTasks = sortingTasks;
                sortingThreads[t] = new Thread(new ParameterizedThreadStart(TableSortingThread));
                sortingThreads[t].Priority = ThreadPriority.BelowNormal;
                sortingThreads[t].Start(sorterParams);
            }

            for (int t = 0; t < noSortingThreads; t++)
            {
                sortingThreads[t].Join();
                sortingThreads[t] = null;
                //Console.WriteLine("finished sorting thread " + (t + 1));
            }

            // clean up memory from just-sorted data structures 
            GC.Collect();
        }

        private static void PairSortingThread(object threadParams)
        {
            SorterParams theseParams = (SorterParams)threadParams;
            Queue<SorterTask> sortingTasks = theseParams.waitingTasks;

            bool allDone = false;

            while (!allDone)
            {
                SorterTask sortTask = null;

                lock (sortingTasks)
                {
                    if (sortingTasks.Count == 0)
                        allDone = true;
                    else
                    {
                        sortTask = sortingTasks.Dequeue();
                    }
                }

                if (!allDone)
                {
                    sortTask.merDictionary.Sort();
                }
            }
        }

        private static void TileAndCount(string[] readsFNs, int fileNo, int readsFormat, int merSize, int trimLength, int trimQual, int qualOffset, int noThreads, MerTables merTables, ThreadStats tesselStats)
        {
            string readsFN = readsFNs[fileNo];
            Console.WriteLine("Reading and tiling " + readsFN);
            StreamReader reads = new StreamReader(new FileStream(readsFN, FileMode.Open, FileAccess.Read, FileShare.Read, readBufferSize, FileOptions.SequentialScan), Encoding.ASCII);
            //StreamReader reads = new StreamReader(readsFN, Encoding.ASCII, false, 1000000);
            BufferedReader bufferedReads = new BufferedReader(readsFormat, reads, null);
            fileNo++;

            CountingThreadParams[] countingParams = new CountingThreadParams[noThreads];
            Thread[] countingThreads = new Thread[noThreads];

            for (int t = 0; t < noThreads; t++)
            {
                countingParams[t] = new CountingThreadParams();
                countingParams[t].threadNumber = t;
                countingParams[t].readsFile = bufferedReads;
                countingParams[t].merSize = merSize;
                countingParams[t].merTables = merTables;
                countingParams[t].trimLength = trimLength;
                countingParams[t].trimQual = trimQual;
                countingParams[t].qualOffset = qualOffset;
                countingParams[t].threadStats = new ThreadStats();
                countingThreads[t] = new Thread(new ParameterizedThreadStart(CountingThread));
                countingThreads[t].Priority = ThreadPriority.BelowNormal;
                countingThreads[t].Start(countingParams[t]);
            }

            for (int t = 0; t < noThreads; t++)
            {
                countingThreads[t].Join();

                tesselStats.noOfReadsRead += countingParams[t].threadStats.noOfReadsRead;
                tesselStats.noOfTiledMers += countingParams[t].threadStats.noOfTiledMers;
                tesselStats.noOfNMers += countingParams[t].threadStats.noOfNMers;
                tesselStats.noOfPoorQuals += countingParams[t].threadStats.noOfPoorQuals;

                countingThreads[t] = null;
                //Console.WriteLine("finished counting thread " + (t + 1));
            }

            bufferedReads.Close();

            //long repeatedUniqueMers = 0;
            //long repeatedMerCount = 0;
            //long totalRepeatedUniqueMers = 0;
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //{
            //    repeatedUniqueMers += merTables.repeatedMers[p].Count;
            //    totalRepeatedUniqueMers += merTables.repeatedMers[p].entries.Length;
            //    foreach (KeyValuePair<ulong, long> kvp in merTables.repeatedMers[p])
            //    {
            //        long count = kvp.Value;
            //        long rcCount = (long)count & 0xffffffff;
            //        long plusCount = (long)count >> 32;
            //        repeatedMerCount += rcCount + plusCount;
            //    }
            //}
            //Console.WriteLine("repeats: " + repeatedUniqueMers + "/" + totalRepeatedUniqueMers + " unique mers with " + repeatedMerCount + " repeats");
            //for (int p = 0; p < merTables.noOfPartitions; p++)
            //    Console.WriteLine("p" + p + ": " + merTables.repeatedMers[p].Count + "/" + merTables.repeatedMers[p].entries.Length + " repeats");
            //for (int s = 0; s < merTables.noSingletonPartitions; s++)
            //    Console.WriteLine("s" + s + ": " + merTables.singletonFilters[s].Count + "/" + merTables.singletonFilters[s].entries.Length + " singles, " + merTables.flushSingletonNumber[s] + " flushes");
            //for (int t = 0; t < noThreads; t++)
            //{
            //    if (merTables.overflowMers[t] != null)
            //        Console.WriteLine("o" + t + ": " + "overflow=" + merTables.overflowMers[t].Count);
            //}

            // don't clean up for either the first or last files - seems to take longer but not reduce memory usage by much at all
            if (fileNo < readsFNs.Length && fileNo != 1)
            {
                //Console.WriteLine("memory before condensing tables: " + GC.GetTotalMemory(false));
                merTables.FlushLowRepMers(merTables, fileNo);
                // and tidy up the heap after all this table clean-up
                GC.Collect(2, GCCollectionMode.Optimized);
                //Console.WriteLine("memory after condensing tables: " + GC.GetTotalMemory(false));
            }

            // wait for flushing to finish between files (don't let the queue grow indefinitely)
            if (fileNo < readsFNs.Length)
                merTables.WaitForSingletonBufferFlush();
        }

        private static void WriteHistogram(string histoFN, Dictionary<int, long> sumReps, string myProcessNameAndArgs, Stopwatch countingTimer, Stopwatch mergingTimer, Stopwatch sortingTimer, ThreadStats tesselStats)
        {
            // Extract, sort and write histogram for original and sum counts to a text file
            Console.WriteLine("Generating histogram and stats...");
            int[] sums = new int[sumReps.Count];
            long[] repsReps = new long[sumReps.Count];
            sumReps.Keys.CopyTo(sums, 0);
            sumReps.Values.CopyTo(repsReps, 0);
            Array.Sort(sums, repsReps);

            StreamWriter histo = new StreamWriter(File.Open(histoFN, FileMode.Create, FileAccess.Write));

            histo.WriteLine(">" + myProcessNameAndArgs + " (Tessel " + version + ")");
            long totalMers = tesselStats.noOfMersWritten + tesselStats.noOfMersDropped;
            histo.WriteLine(">copies\tcounts\t" + totalMers);
            for (int hi = 0; hi < sums.Length; hi++)
            {
                histo.Write(sums[hi]);
                histo.Write('\t');
                histo.Write(repsReps[hi]);
                histo.Write('\t');
                long mersInBucket = sums[hi] * repsReps[hi];
                histo.Write(mersInBucket);
                histo.Write('\t');
                histo.Write((((float)mersInBucket / (float)totalMers) * 100.0).ToString("F2") + "%");
                histo.WriteLine();
            }

            histo.WriteLine();
            histo.WriteLine(tesselStats.noOfReadsRead + "\treads");
            histo.WriteLine(tesselStats.noOfTiledMers + "\tmers tiled from reads");
            histo.WriteLine(tesselStats.noOfNMers + "\tN-mers");
            histo.WriteLine(tesselStats.noOfValidMers + "\tvalid mers");
            histo.WriteLine(tesselStats.noOfDistinctMers + "\tdistinct mers written to cbt file");
            histo.WriteLine(tesselStats.noOfMersWritten + "\ttotal mers written to cbt file");
            histo.WriteLine(tesselStats.noOfMersDropped+ "\tmers dropped (too few reps)");
            histo.WriteLine(tesselStats.noOfPoorQuals + "\tpoor quality bases trimmed from ends of reads");
            histo.WriteLine(tesselStats.noOfDistinctPairs + "\tdistinct pairs written to prs file");
            histo.WriteLine(tesselStats.noOfPairsWritten + "\ttotal pairs written to prs file");
            histo.WriteLine(countingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts counting");
            histo.WriteLine(mergingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts merging");
            histo.WriteLine(sortingTimer.Elapsed.TotalSeconds.ToString("#.0") + "\ts sorting");
            histo.WriteLine("Filling took " + timeSpentFilling.Elapsed.TotalSeconds.ToString("#.0") + "s");
            histo.WriteLine("Sorting took " + timeSpentSorting.TotalSeconds.ToString("#.0") + "s");
            histo.WriteLine("Writing took " + timeSpentWriting.TotalSeconds.ToString("#.0") + "s");

            //histo.WriteLine(merTables.singles + "\tsingles\t" + merTables.promoted + "\tpromoted\t" + merTables.repeatsFull + "\trepeatsFull\t" + merTables.repeatsQuick + "\trepeatsQuick\t" + merTables.flushes + "\tflushes");
            //for (int t = 0; t < noThreads; t++)
            //{
            //    histo.WriteLine(threadLockTimers[t].Elapsed.TotalSeconds.ToString("#.000") + "\t secs waiting by thread " + t);
            //}

            //histo.WriteLine("Partition stats");
            //histo.WriteLine("p#\tsingles\tdictionary\trepeats\tpromoted\tflushes\tresizes\tAs\tdict0\twaiting");
            //for (int p = 0; p < merSets.noOfPartitions; p++)
            //{
            //    histo.Write(p + "\t");
            //    histo.Write(merSets.pSingles[p] + "\t");
            //    histo.Write(merSets.pDictionary[p] + "\t");
            //    histo.Write(merSets.repeatedMers[p].Count + "\t");
            //    histo.Write(merSets.pPromoted[p] + "\t");
            //    histo.Write(merSets.pFlushes[p] + "\t");
            //    histo.Write(merSets.repeatedMers[p].resizeCount + "\t");
            //    histo.Write(merSets.As[p] + "\t");
            //    if (merSets.repeatedMers[p].ContainsKey(0))
            //        histo.Write((merSets.repeatedMers[p][0]>>32) + "\t" + (merSets.repeatedMers[p][0] & 0xffffffff) + "\t");
            //    else
            //        histo.Write("\t\t");
            //    histo.WriteLine(merSets.pWaitingTimer[p].Elapsed.TotalSeconds.ToString("#.000"));
            //}

            histo.Close();
        }

        private static void TableSortingThread(object threadParams)
        {
            SorterParams theseParams = (SorterParams)threadParams;
            Queue<SorterTask> sortingTasks = theseParams.waitingTasks;

            bool allDone = false;

            while (!allDone)
            {
                SorterTask sortTask = null;

                lock (sortingTasks)
                {
                    if (sortingTasks.Count == 0)
                        allDone = true;
                    else
                    {
                       sortTask = sortingTasks.Dequeue();
                    }
                }

                if (!allDone)
                {
                    if (sortTask.merCollection != null)
                    {
                        //Console.WriteLine("sorting singletons (" + sortTask.merCollection.Count + ")");
                        sortTask.merCollection.Condense();              // remove inactive entries
                        sortTask.merCollection.Sort();                  // sort the active entries
                        sortTask.merCollection.ExtendAndSaveRC();       // extend the singletons to full length and save the RC bits
                        //Console.WriteLine("sorted singletons (" + sortTask.merCollection.Count + ")");
                    }
                    if (sortTask.merDictionary != null)
                        sortTask.merDictionary.Sort();
                }
            }
        }

        private static void StopMonitorThread(Thread monitorProgress)
        {
            stopMonitor = true;
            monitorProgress.Abort();
            monitorProgress.Join();
        }
        
        private static void CountingThread(object threadParams)
        {
            CountingThreadParams theseParams = (CountingThreadParams)threadParams;
            int threadNo = theseParams.threadNumber;
            BufferedReader readsFile = theseParams.readsFile;
            int merSize = theseParams.merSize;
            MerTables merTables = theseParams.merTables;
            int trimLength = theseParams.trimLength;
            int trimQual = theseParams.trimQual;
            int qualOffset = theseParams.qualOffset;
            ThreadStats threadStats = theseParams.threadStats;

            ulong[] merSet = new ulong[maxReadSize];
            bool[] merValid = new bool[maxReadSize];

            Sequence[] readHeaderBatch = new Sequence[batchSize];
            Sequence[] readBatch = new Sequence[batchSize];
            Sequence[] qualHeaderBatch = new Sequence[batchSize];
            Sequence[] qualBatch = new Sequence[batchSize];
            for (int i = 0; i < batchSize; i++)
            {
                readHeaderBatch[i] = new Sequence(defaultHeaderLength);
                readBatch[i] = new Sequence(defaultReadLength);
                qualHeaderBatch[i] = new Sequence(defaultHeaderLength);
                qualBatch[i] = new Sequence(defaultReadLength);             
            }

            //ulong testMer;
            //MerStrings.CondenseMer("CCTAAAAAAAAAAAAAAAAAAAAAA", out testMer);
            //ulong testMerRC = MerStrings.ReverseComplement(testMer);

            bool EOF = false;
            int readsInBatch = 0;

            while (!EOF)
            {
                lock (readsFile)
                {
                    readsInBatch = readsFile.ReadReads(batchSize, readHeaderBatch, readBatch, qualHeaderBatch, qualBatch);

                    if (readsInBatch != batchSize)
                        EOF = true;

                    progressReadsRead += readsInBatch;
                }

                for (int r = 0; r < readsInBatch; r++)
                {
                    threadStats.noOfReadsRead++;

                    int readLength = readBatch[r].Length;
                    if (readLength < merSize)
                        continue;
                    if (trimLength != 0)
                        readBatch[r].Length = readLength - trimLength;

                    if (trimQual > 0)
                        threadStats.noOfPoorQuals += SeqFiles.TrimTrailingPoorQuals(readBatch[r], qualBatch[r], trimQual, qualOffset);

                    int mersInRead = Sequence.GenerateMersFromRead(readBatch[r], merSize, ref merSet, ref merValid);

                    for (int i = 0; i < mersInRead; i++)
                    {
                        threadStats.noOfTiledMers++;
                       
                        if (merValid[i])
                        {
                            merTables.AddOrIncrement(merSet[i], merSize, threadNo);
                            threadStats.noOfValidMers++;
                        }
                        else
                            threadStats.noOfNMers++;
                    } // for each mer in the read

                } // for each read in the batch

            } // while there are still reads in the file

        }

        static void PairCountingThread(object param)
        {
            PairCountingThreadParams threadParams = (PairCountingThreadParams)param;
            int threadNo = threadParams.threadNumber;
            BufferedReader readsFile = threadParams.readsFile;
            int merSize = threadParams.merSize;
            int minReadLength = threadParams.minReadLength;
            MerCollections.PairTable pairTable = threadParams.pairsTable;
            kMerTable kMerTable = threadParams.kMersTable;
            int pairGap = threadParams.pairGap;
            int minCount = threadParams.minCount;
            ThreadStats threadStats = threadParams.threadStats;

            bool EOF = false;

            Sequence[] readHeaderBatch = new Sequence[batchSize];
            Sequence[] readBatch = new Sequence[batchSize];
            Sequence[] qualHeaderBatch = new Sequence[batchSize];
            Sequence[] qualsBatch = new Sequence[batchSize];
            for (int i = 0; i < batchSize; i++)
            {
                readHeaderBatch[i] = new Sequence(defaultHeaderLength);
                readBatch[i] = new Sequence(defaultReadLength);
                qualHeaderBatch[i] = new Sequence(defaultHeaderLength);
                qualsBatch[i] = new Sequence(defaultReadLength);
            }
            int readsInBatch = 0;
            long readsRead = 0;
            long pairsTiled = 0;

            ulong[] mersFromRead = new ulong[1000];
            bool[] mersValid = new bool[1000];
            ulong[] pairsFromRead = new ulong[1000];
            bool[] pairsValid = new bool[1000];

            while (!EOF)
            {
                lock (readsFile)
                {
                    readsInBatch = readsFile.ReadReads(batchSize, readHeaderBatch, readBatch, qualHeaderBatch, qualsBatch);
                }

                if (readsInBatch != batchSize)
                    EOF = true;

                readsRead += readsInBatch;
                progressReadsRead += readsInBatch;

                for (int r = 0; r < readsInBatch; r++)
                {
                    Sequence read = readBatch[r];
                    int readLength = read.Length;

                    if (readLength < minReadLength)
                    {
                        continue;
                    }

                    //string target = "CTGGAAGTTGCCAGCTGGCTGGAAGATGTAGATGGCAAACAGGAAACGCGTTATGCCTTTATTGATGAGGCCGATAATAAAACAGAGGATTCTCTGAAGGCTGCGAAGGAGAAAATTTTCGCCGCGTTCCCGGGGCTGAAAGAGTGTACT";
                    //if (read.ToString() == target)
                    //    Debugger.Break();

                    // get kMers and depths from read so we can check presence in the tiled kMers table
                    int merCount = Sequence.GenerateMersFromRead(read, merSize, ref mersFromRead, ref mersValid);
                    int pairsInRead = kMerPairs.GeneratePairsFromRead(read, pairGap, ref pairsFromRead, ref pairsValid);

                    for (int i = 0; i < pairsInRead; i++)
                    {
                        ulong pair = pairsFromRead[i];
                        bool pairValid = pairsValid[i];

                        if (pairValid)
                        {
                            // pair table contains canonical kMers
                            ulong rcPair = kMers.ReverseComplement(pair, 32);
                            if (rcPair < pair)
                                pair = rcPair;

                            // increment the pair entry if it exists
                            if (pairTable.IncrementIfPresent(pair, threadNo))
                            {
                                pairsTiled++;
                                continue;
                            }

                            // check that the pair is based on sound kMers, and add it to the table 
                            int firstMerIdx = i;
                            int secondMerIdx = i + kMerPairs.pairFragmentSize + pairGap - (merSize - kMerPairs.pairFragmentSize);

                            int firstMerDepth = kMerTable.GetDepthSum(mersFromRead[firstMerIdx]);
                            int secondMerDepth = kMerTable.GetDepthSum(mersFromRead[secondMerIdx]);

                            // slip this pair if either of the kMer depths is below minCount
                            if (firstMerDepth < minCount || secondMerDepth < minCount)
                                continue;

                            //if (pair == 0x000015100F010015)
                            //    Debugger.Break();

                            // count the number of times each pair is seen
                            pairTable.AddOrIncrement(pair, threadNo);
                            pairsTiled++;
                        } // valid pair
                    } // for each pair
                } // for each read in batch
            } // EOF loop

            threadStats.noOfTiledPairs = pairsTiled;
        }

        private static void GetDepthsForRead(kMerTable kMerTable, int merCount, ulong[] mersFromRead, bool[] merValid, int[] depths)
        {
            for (int m = 0; m < merCount; m++)
            {
                int depth;

                if (merValid[m])
                {
                    ulong mer = mersFromRead[m];
                    depth = kMerTable.GetDepthSum(mer);
                }
                else
                {
                    depth = 0;
                }

                depths[m] = depth;
            }
        }


        private static void MergeAndWrite(int noThreads, string bwFN, TypeOfCount typeOfCount, int minCount, int flagInt, MerDictionary[] repeatedMers, MerDictionary[] overflowMers, SingletonCollection[] singletons, List<SingletonCollection>[] deferredSingletons, 
                                            List<string>[] singletonFNs, List<ulong>[] lowestFlushedMer, List<ulong>[] highestFlushedMer, List<int>[] flushedSingletonsCount,
                                            List<string> flushedLowRepsFNs, List<ulong> firstFlushedLowRepsMer, List<ulong> lastFlushedLowRepsMer, List<int> flushedLowRepsCount, ThreadStats tesselStats, Dictionary<int, long> sumReps)
        {
            // open the .cbt/prs file 
            BinaryWriter bwFile = new BinaryWriter(new FileStream(bwFN, FileMode.Create, FileAccess.Write, FileShare.None, 100000));

            // first 4 bytes of the .cbt file is the k-mer size (could version the file format using upper bits of this length word)
            // first 4 bytes of a .prs file is the pairGap
            bwFile.Write(flagInt);

            int noOfSingletons = 0;
            if (singletons != null)
                noOfSingletons = singletons.Length;
            int noOfFlushFiles = 0;
            if (singletonFNs != null)
                for (int p = 0; p < singletonFNs.Length; p++)
                {
                    int flushFilesInPartition = singletonFNs[p].Count;
                    noOfFlushFiles += flushFilesInPartition;
                }
            int noOfOverflows = 0;
            if (overflowMers != null)
                for (int p = 0; p < overflowMers.Length; p++)
                    if (overflowMers[p] != null)
                        noOfOverflows++;
            int noOfDeferred = 0;
            if (deferredSingletons != null)
            for (int p = 0; p < deferredSingletons.Length; p++)
                if (deferredSingletons[p] != null)
                    noOfDeferred += deferredSingletons[p].Count;
            int noOfFlushedLowRepsFNs = 0;
            if (flushedLowRepsFNs != null)
                noOfFlushedLowRepsFNs = flushedLowRepsFNs.Count;

            //                 shared mers           overflow        unflushed        flushed          flushed lowrep mers     deferred singletons
            int noMerSources = repeatedMers.Length + noOfOverflows + noOfSingletons + noOfFlushFiles + noOfFlushedLowRepsFNs + noOfDeferred;
            MerSource[] merSources = new MerSource[noMerSources];
            List<MerSource> openMerSources = new List<MerSource>(noMerSources);

            int nextSource = 0;
            // shared repeated mers partitions
            if (repeatedMers != null)
                for (int i = 0; i < repeatedMers.Length; i++)
                {
                    //Console.WriteLine(nextSource + ": repeatedMers[" + i + "] holds " + repeatedMers[i].Count + "/" + repeatedMers[i].Capacity + " k-mers");
                    merSources[nextSource] = new MerDictionarySource(repeatedMers[i]);
                    merSources[nextSource].sourceNo = nextSource;         
                    merSources[nextSource].Open();
                    openMerSources.Add(merSources[nextSource]);
                    nextSource++;
                }
            //Console.WriteLine("created " + repeatedMers.Length + " MerDictionarySources (primary)");

            // all the overflow mer tables
            if (overflowMers != null)
                for (int i = 0; i < overflowMers.Length; i++)
                {
                    if (overflowMers[i] != null)
                    {
                        //Console.WriteLine(nextSource + ": overflowMers[" + i + "] holds " + overflowMers[i].Count + "/" + overflowMers[i].Capacity + " k-mers");
                        merSources[nextSource] = new MerDictionarySource(overflowMers[i]);
                        merSources[nextSource].sourceNo = nextSource;           
                        merSources[nextSource].Open();
                        openMerSources.Add(merSources[nextSource]);
                        nextSource++;
                    }
                }
            //Console.WriteLine("created " + noOfOverflows + " MerDictionarySources (overflows)");

            // all the unflushed singleton tables (but only open the first one)
            if (singletons != null)
                for (int p = 0; p < singletons.Length; p++)
                {
                    if (singletons[p].Count > 0)
                    { 
                        //Console.WriteLine(nextSource + ": singletons[" + p + "] holds " + singletons[p].Count + "/" + singletons[p].Capacity + " k-mers");
                        merSources[nextSource] = new MerSingletonsSource(singletons[p], p);
                        merSources[nextSource].sourceNo = nextSource;     
                        // only open the first partition
                        if (p == 0)
                        {
                            merSources[nextSource].Open();

                            openMerSources.Add(merSources[nextSource]);
                        }
                        nextSource++;
                    }
                }
            //Console.WriteLine("created " + singletons.Length + " unflushed singleton tables");

            // any deferred flushed singleton tables (but only open those from the first partition)
            int deferredSingletonsCount = 0;
            if (singletons != null)
                for (int p = 0; p < singletons.Length; p++)
                {
                    if (deferredSingletons[p] != null)
                    {
                        for (int i = 0; i < deferredSingletons[p].Count; i++)
                        {
                            if (deferredSingletons[p][i].Count > 0)
                            {
                                //Console.WriteLine(nextSource + ": deferred singletons[" + p + "," + i + "] holds " + deferredSingletons[p][i].Count + "/" + deferredSingletons[p][i].Capacity + " k-mers");
                                merSources[nextSource] = new MerSingletonsSource(deferredSingletons[p][i], p);
                                merSources[nextSource].sourceNo = nextSource;
                                deferredSingletonsCount++;
                                // only open those from first partition
                                if (p == 0)
                                {
                                    merSources[nextSource].Open();
                                    openMerSources.Add(merSources[nextSource]);
                                }
                                nextSource++;
                            }
                        }
                    }
                }
            //Console.WriteLine("created " + deferredSingletonsCount + " deferred singleton tables");

            // all the flushed low-rep mers
            if (flushedLowRepsFNs != null)
                for (int i = 0; i < flushedLowRepsFNs.Count; i++)
                {
                    //Console.WriteLine(nextSource + ": flushedLowReps[" + i + "] holds " + flushedLowRepsCount[i] + " k-mers");
                    merSources[nextSource] = new MerFlushedLowRepSource(flushedLowRepsFNs[i], firstFlushedLowRepsMer[i], lastFlushedLowRepsMer[i]);
                    merSources[nextSource].sourceNo = nextSource;
                    merSources[nextSource].Open();
                    openMerSources.Add(merSources[nextSource]);
                    nextSource++;
                }
            //Console.WriteLine("created " + flushedLowRepsFNs.Count + " MerFlushedLowRepSources");

            // and the flushed singleton files (first partition starts open, all others start being closed, and opened only when needed)
            int totalSingletonFiles = 0;
            if (singletonFNs != null)
                for (int p = 0; p < singletonFNs.Length; p++)
                {
                    for (int f = 0; f < singletonFNs[p].Count; f++)
                    {
                        //Console.WriteLine(nextSource + ": flushedSingletons[" + p + "," + f + "] holds " + flushedSingletonsCount[p][f] + " k-mers");
                        string singletonFN = singletonFNs[p][f];
                        merSources[nextSource] = new MerFlushedSingletonSource(singletonFN, p, lowestFlushedMer[p][f], highestFlushedMer[p][f], flushedSingletonsCount[p][f]);        
                        merSources[nextSource].sourceNo = nextSource;
                        // only open the merSource if this is the first partition
                        if (p == 0)
                        {
                            merSources[nextSource].Open();
                            openMerSources.Add(merSources[nextSource]);
                        }
                        nextSource++;
                        totalSingletonFiles++;
                    }
                }
            //Console.WriteLine("created " + totalSingletonFiles + " MerFlushedSingletonSources");

            //for (int i = 0; i < nextSource; i++)
            //    Console.WriteLine(i + " " + merSources[i].sourceTypeString + merSources[i].opened + " " + merSources[i].lowestMer.ToString("x16") + " " + merSources[i].firstMer.ToString("x16") + " " + merSources[i].lastMer.ToString("x16"));

            //Console.WriteLine("memory at start of merging: " + GC.GetTotalMemory(false));

            int noSortingThreads = 2;
            int noMerBuffers = noSortingThreads + 4;

            Queue<MerBuffer> buffersToBeSorted = new Queue<MerBuffer>(noMerBuffers);
            Queue<MerBuffer> buffersToBeWritten = new Queue<MerBuffer>(noMerBuffers);
            MerBuffer[] merBuffers = new MerBuffer[noMerBuffers];
            for (int t = 0; t < noMerBuffers; t++)
                merBuffers[t] = new MerBuffer(t, targetBufferSize);

            //Console.WriteLine("memory after merBuffers allocated: " + GC.GetTotalMemory(false));

            // event for signalling sorters when buffer placed on queue. AutoReset to wake only one Sorter at a time
            EventWaitHandle sortingEWH = new EventWaitHandle(false, EventResetMode.AutoReset);
            // event for signalling writer that a buffer has been sorted and placed on the writing queue
            EventWaitHandle writingEWH = new EventWaitHandle(false, EventResetMode.AutoReset);
            // event for signalling that there are free buffers available (set by Writer, waited on by Fill loop)
            EventWaitHandle freedEWH = new EventWaitHandle(true, EventResetMode.AutoReset);

            // start the sorting threads
            SortingThreadParams[] sortingParams = new SortingThreadParams[noSortingThreads];
            Thread[] sortingThreads = new Thread[noSortingThreads];

            for (int t = 0; t < noSortingThreads; t++)
            {
                sortingParams[t] = new SortingThreadParams();
                sortingParams[t].threadNumber = t;
                sortingParams[t].sortingEWH = sortingEWH;
                sortingParams[t].writingEWH = writingEWH;
                sortingParams[t].buffersToBeSorted = buffersToBeSorted;
                sortingParams[t].buffersToBeWritten = buffersToBeWritten;
                sortingThreads[t] = new Thread(new ParameterizedThreadStart(BufferSortingThread));
                sortingThreads[t].Priority = ThreadPriority.BelowNormal;
                sortingThreads[t].Start(sortingParams[t]);
            }

            // and start the writing thread
            WritingThreadParams writingParams = new WritingThreadParams();
            writingParams.threadNumber = 0;
            writingParams.buffersToBeWritten = buffersToBeWritten;
            writingParams.writingEWH = writingEWH;
            writingParams.freedEWH = freedEWH;
            writingParams.bwFile = bwFile;
            writingParams.typeOfCount = typeOfCount;
            writingParams.minCount = minCount;
            writingParams.tesselStats = tesselStats;
            writingParams.sumReps = sumReps;
            Thread writingThread = new Thread(new ParameterizedThreadStart(WritingThread));
            writingThread.Priority = ThreadPriority.Normal;
            writingThread.Start(writingParams);

            // now just merge and write until all mers have been written
            bool mersLeft = true;
            ulong highestMerInBuffer = 0;
            int bufferNo = 0;
            int finalBufferNo = 0;
            int nextSingletonPartitionNo = 0;

            // loop around filling buffers and sending them down the sort & write pipeline
            while (mersLeft)
            {
                // get an empty merBuffer (and wait if none are available)
                int bi = GetFreeMerBuffer(freedEWH, merBuffers);
                merBuffers[bi].bufferNo = bufferNo;
                finalBufferNo = bufferNo;
                bufferNo++;
                //ulong startingMerInBuffer = highestMerInBuffer;  // trace only

                // and get the next load of mers to be written from all sources
                mersLeft = FillBuffer(openMerSources, merSources, ref nextSingletonPartitionNo, ref merBuffers[bi], highestMerInBuffer, out highestMerInBuffer);

                // reached the end - the final merBuffer was filled last time around this loop
                if (!mersLeft)
                    break;
                
                lock (buffersToBeSorted)
                {
                    //Console.WriteLine("Fill: queueing buffer #" + merBuffers[bi].bufferNo + ". " + startingMerInBuffer.ToString("X16") + " to " + highestMerInBuffer.ToString("X16"));
                    buffersToBeSorted.Enqueue(merBuffers[bi]);
                }

                // tell the sorters that a buffer is available
                sortingEWH.Set();
            }

            //Console.WriteLine("Fill: end reached with buffer #" + bufferNo);
            // tell all the sorters to finish up by queueing up termination buffers 
            MerBuffer finishBuffer = new MerBuffer(0, 1);
            finishBuffer.finishSorter = true;
            finishBuffer.bufferNo = finalBufferNo;
            for (int t = 0; t < noSortingThreads; t++)
            {
                lock (buffersToBeSorted)
                {
                    //Console.WriteLine("Fill: enqueuing sorting stop buffer #" + finishBuffer.bufferNo);
                    buffersToBeSorted.Enqueue(finishBuffer);
                }
                sortingEWH.Set();
            }
            // and wait for them all to finish
            for (int t = 0; t < noSortingThreads; t++)
            {
                sortingThreads[t].Join();
                sortingThreads[t] = null;
                //Console.WriteLine("finished sorting thread " + (t + 1));
            }

            // tell the Writer to finish 
            finishBuffer.bufferNo = finalBufferNo;
            finishBuffer.finishWriter = true;
            lock (buffersToBeWritten)
            {
                buffersToBeWritten.Enqueue(finishBuffer);
            }
            writingEWH.Set();

            // wait for the writing thread to finish
            writingThread.Join();
            writingThread = null;

            bwFile.Close();
        }

        private static int GetFreeMerBuffer(EventWaitHandle freedEWH, MerBuffer[] merBuffers)
        {
            bool foundFreeBuffer = false;
            int freeBufferIdx = -1;

            while (!foundFreeBuffer)
            {
                for (int b = 0; b < merBuffers.Length; b++)
                    if (merBuffers[b].free) 
                    {
                        freeBufferIdx = b;
                        foundFreeBuffer = true;
                        break;
                    }
                if (!foundFreeBuffer)
                {
                    //Console.WriteLine("waiting for free buffer");
                    freedEWH.WaitOne();
                }
            }

            merBuffers[freeBufferIdx].free = false;
            //Array.Clear(merBuffers[freeBufferIdx].mers, 0, merBuffers[freeBufferIdx].mers.Length);
            //Console.WriteLine("allocated buffer " + freeBufferIdx + " [" + merBuffers[freeBufferIdx].mers.Length + "]");
            return freeBufferIdx;
        }

        private static bool FillBuffer(List<MerSource> openMerSources, MerSource[] allMerSources, ref int nextSingletonPartitionNo, ref MerBuffer merBuffer, ulong startingMer, out ulong endingMer)
        {
            timeSpentFilling.Start();

            merBuffer.activeMers = 0;
            int bufferIdx = 0;

            bool foundSomeMers = false;
            ulong lastLoadedInitialMer = ulong.MaxValue;

            int firstOpenSingletonSourceIdx = -1;
            int firstOpenSourceIdx = -1;
            int openSingletonSourcesCount = 0;

            // find the first open, valid singleton source (if there is one)
            for (int s = 0; s < openMerSources.Count; s++)
            {
                if (openMerSources[s].opened && (openMerSources[s] is MerFlushedSingletonSource || openMerSources[s] is MerSingletonsSource))
                {
                    if (firstOpenSingletonSourceIdx < 0)
                    {
                        firstOpenSingletonSourceIdx = s;
                        firstOpenSourceIdx = s;
                    }
                    openSingletonSourcesCount++;
                }
            }

            // no open singletons - could have run into the end of a partition so open the next partition if this has happened
            if (openSingletonSourcesCount == 0)
            {
                nextSingletonPartitionNo++;

                // remove any now-closed sources from the MerSources set (and allow their memory to be GCed)
                for (int ams = 0; ams < allMerSources.Length; ams++)
                    if (allMerSources[ams] != null && allMerSources[ams].closed)
                        allMerSources[ams] = null;

                // open the next set of singletons partitions
                for (int ams = 0; ams < allMerSources.Length; ams++)
                    if (allMerSources[ams] != null && allMerSources[ams].partition == nextSingletonPartitionNo)
                            allMerSources[ams].Open();

                // and rebuild the open sources list
                openMerSources.Clear();
                for (int ams = 0; ams < allMerSources.Length; ams++)
                    if (allMerSources[ams] != null && allMerSources[ams].opened)
                        openMerSources.Add(allMerSources[ams]);

                // and try again to find the first open singleton source (or any other source in this case)
                for (int s = 0; s < openMerSources.Count; s++)
                {
                    if (firstOpenSourceIdx < 0)
                        firstOpenSourceIdx = s;
                    if (openMerSources[s] is MerFlushedSingletonSource || openMerSources[s] is MerSingletonsSource)
                    {
                        if (firstOpenSingletonSourceIdx < 0)
                            firstOpenSingletonSourceIdx = s;
                        openSingletonSourcesCount++;
                    }
                }

                // prefer singleton source
                if (firstOpenSingletonSourceIdx > 0)
                    firstOpenSourceIdx = firstOpenSingletonSourceIdx;
            }

            // no open sources so return
            if (firstOpenSourceIdx < 0)
            {
                endingMer = ulong.MaxValue;
                return false;
            }

            // initial buffer loading from the first active source
            MerSource firstSource = openMerSources[firstOpenSourceIdx];
            foundSomeMers = true;

            //Console.WriteLine("Filling #" + merBuffer.bufferNo + " from source " + initialLoadingSourceNo + " " + firstSource.sourceTypeString + " " + firstSource.lowestMer.ToString("X16") + " (" + firstSource.firstMer.ToString("X16") + " - " + firstSource.lastMer.ToString("X16") + ")");

            // resize the buffers if necessary
            int initialBufferSize = targetBufferSize / (openSingletonSourcesCount + 4);

            // get a buffer load of mers from this initial source
            for (int i = 0; i < initialBufferSize; i++)
            {
                merBuffer.mers[bufferIdx] = firstSource.lowestMer;
                merBuffer.counts[bufferIdx] = firstSource.countPair;
                lastLoadedInitialMer = firstSource.lowestMer;
                bufferIdx++;

                // reached the end of this source (will set 'opened' to false)
                if (!firstSource.MoveToNextMer())
                    break;
            }

            // and load any duplicates of the highest value as well - only possible if there was an insert race at just this mer
            while (firstSource.opened && firstSource.lowestMer == lastLoadedInitialMer)
            {
                merBuffer.mers[bufferIdx] = firstSource.lowestMer;
                merBuffer.counts[bufferIdx] = firstSource.countPair;
                bufferIdx++;
                if (!firstSource.MoveToNextMer())
                    break;
            }

            endingMer = lastLoadedInitialMer;
            //Console.WriteLine("buffer: " + startingMer.ToString("X16") + " " + endingMer.ToString("X16"));

            // now add the matching mers/counts from all the other open sources 
            for (int s = 0; s < openMerSources.Count; s++)
            {
                MerSource ms = openMerSources[s];        
                
                // don't load the initial source twice
                if (s == firstOpenSourceIdx)
                    continue;

                // and skip any that are actually closed (but are in the open list until it is rebuilt when the next partition is added)
                if (!ms.opened)
                    continue;

                while (lastLoadedInitialMer >= ms.lowestMer)
                {
                    if (bufferIdx == merBuffer.mers.Length)
                    {
                        Array.Resize<ulong>(ref merBuffer.mers, (merBuffer.mers.Length + merBuffer.mers.Length / 4));
                        Array.Resize<ulong>(ref merBuffer.counts, merBuffer.mers.Length);
                        //Console.WriteLine("resizing buffer " + merBuffer.bufferNo + " during fill to " + merBuffer.mers.Length);
                    }

                    merBuffer.mers[bufferIdx] = ms.lowestMer;
                    merBuffer.counts[bufferIdx] = ms.countPair;
                    bufferIdx++;

                    if (!ms.MoveToNextMer())
                    {
                        //Console.WriteLine("end reached for source " + s + " " + ms.sourceTypeString + " (" + ms.firstMer.ToString("X16") + " - " + ms.lastMer.ToString("X16"));
                        break;
                    }
                }
            }

            merBuffer.activeMers = bufferIdx;

            timeSpentFilling.Stop();

            return foundSomeMers;
        }

        private static void BufferSortingThread(object threadParams)
        {
            SortingThreadParams theseParams = (SortingThreadParams)threadParams;
            int threadNo = theseParams.threadNumber;
            EventWaitHandle sortingEWH = theseParams.sortingEWH;
            EventWaitHandle writingEWH = theseParams.writingEWH;
            Queue<MerBuffer> buffersToBeSorted = theseParams.buffersToBeSorted;
            Queue<MerBuffer> buffersToBeWritten = theseParams.buffersToBeWritten;
            Stopwatch localTimeSpentSorting = new Stopwatch();

            bool stopSorter = false;

            //Console.WriteLine("starting sorter #" + threadNo);

            // loop until we're told to stop (and the sorting queue is empty)
            while (!stopSorter)
            {
                // wait for something to be enqueued on the sorting queue 
                sortingEWH.WaitOne();
                //Console.WriteLine("sorter[" + threadNo + "]: woken");

                MerBuffer bufferToSort = null;      
          
                // get the buffer to sort - there will always be a buffer in the queue 
                lock (buffersToBeSorted)
                {
                    if (buffersToBeSorted.Count == 0)
                    {
                        //Console.WriteLine("sort[" + threadNo + "]: queue empty");
                        continue;
                    }
                    bufferToSort = buffersToBeSorted.Dequeue();
                    //Console.WriteLine("sort[" + threadNo + "]: dequeued buffer #" + bufferToSort.bufferNo);
                    // signal another Sorter if there are other buffers still queued (as only this thread will have been woken)
                    if (buffersToBeSorted.Count > 0)
                    {
                        //Console.WriteLine("sort[" + threadNo + "]: signalling next sorter");
                        sortingEWH.Set();
                    }
                }

                if (bufferToSort.finishSorter)
                {
                    //Console.WriteLine("sort[" + threadNo + "]: found stop buffer #" + bufferToSort.bufferNo);
                    stopSorter = true;
                }

                // sort the buffer entries in mer order (unless this is just the termination buffer)
                if (!stopSorter)
                {
                    localTimeSpentSorting.Start();
                    //Array.Sort<ulong, ulong>(bufferToSort.mers, bufferToSort.counts, 0, bufferToSort.activeMers);
                    bufferToSort.Sort();
                    localTimeSpentSorting.Stop();

                    // and put the sorted buffer on the writing queue and wake up any waiting writer
                    lock (buffersToBeWritten)
                    {
                        //Console.WriteLine("sort[" + threadNo + "]: enqueuing buffer #" + bufferToSort.bufferNo + " for writing. "
                        //                    + bufferToSort.mers[0].ToString("X16") + " - " + bufferToSort.mers[bufferToSort.activeMers - 1].ToString("X16"));
                        buffersToBeWritten.Enqueue(bufferToSort);
                    }
                    writingEWH.Set();
                }
            }

            lock (buffersToBeSorted)
            {
                timeSpentSorting += localTimeSpentSorting.Elapsed;
            }
            //Console.WriteLine("exiting sorter #" + threadNo);
        }

        // in-lined for performance
        //private static void GetLowestMer(ulong[] bufferMers, ulong[] bufferCountPairs, ref int lowestIdx, int bufferCount, out ulong lowestMer, out ulong lowestCountPair)
        //{
        //    lowestMer = bufferMers[lowestIdx];
        //    lowestCountPair = bufferCountPairs[lowestIdx];
        //
        //    lowestIdx++;
        //    while (lowestIdx < bufferCount && bufferMers[lowestIdx] == lowestMer)
        //    {
        //        lowestCountPair += bufferCountPairs[lowestIdx];
        //        lowestIdx++;
        //    }
        //}

        private static void WritingThread(object threadParams)
        {
            WritingThreadParams theseParams = (WritingThreadParams)threadParams;
            int threadNo = theseParams.threadNumber;
            EventWaitHandle writingEWH = theseParams.writingEWH;
            Queue<MerBuffer> buffersToBeWritten = theseParams.buffersToBeWritten;
            int minCount = theseParams.minCount;
            ThreadStats tesselStats = theseParams.tesselStats;
            Dictionary<int, long> sumReps = theseParams.sumReps;
            BinaryWriter cbtFile = theseParams.bwFile;
            TypeOfCount typeOfCount = theseParams.typeOfCount;
            EventWaitHandle freedEWH = theseParams.freedEWH;
            //List<string> bufferTrace = new List<string>();

            //Console.WriteLine("starting writer thread");

            // buffers have to be written in order
            int nextBufferNoToWrite = 0;
            // so out-of-order buffers are saved until their turn comes
            SortedList<int, MerBuffer> pendingBuffers = new SortedList<int, MerBuffer>(10);

            // last mer written from previous call - used to check ordering constraint
            ulong lastMerWritten = 0;

            Stopwatch localTimeSpentWriting = new Stopwatch();

            bool stopWriter = false;

            // loop until we've written the final buffer
            while (!stopWriter)
            {
                // wait for something to be enqueued on the writing queue 
                writingEWH.WaitOne();
                //Console.WriteLine("writer: woken");

                // copy buffers from the to-be-written queue (unordered) to the (ordered) list of pending buffers
                lock (buffersToBeWritten)
                {
                    bool bufferQueueEmpty = buffersToBeWritten.Count == 0;
                    while (!bufferQueueEmpty)
                    {
                        // get the buffer to write - there will always be at least one buffer in the queue if the event has been set
                        MerBuffer bufferToWrite = buffersToBeWritten.Dequeue();
                        // add it to the list of buffers to be written
                        pendingBuffers.Add(bufferToWrite.bufferNo, bufferToWrite);
                        bufferQueueEmpty = buffersToBeWritten.Count == 0;
                        //Console.WriteLine("writer: added buffer #" + bufferToWrite.bufferNo + " to pending list (" + buffersToBeWritten.Count + " still queued)");
                    }
                }

                // and write as many of the pending buffers as possible (we could have pending buffers waiting for the new arrivals)
                bool writtenAll = false;
                while (!writtenAll)
                {
                    // is the next buffer to be written in the pending list?
                    if (pendingBuffers.ContainsKey(nextBufferNoToWrite))
                    {
                        MerBuffer nextBufferToWrite = pendingBuffers[nextBufferNoToWrite];
                        // only write the buffer if it is a real one
                        if (nextBufferToWrite.finishWriter)
                        {
                            //Console.WriteLine("writer: received stop buffer");
                            pendingBuffers.Remove(nextBufferNoToWrite);
                            stopWriter = true;
                        }
                        else
                        {
                            //string writeTrace = "buffer #" + nextBufferToWrite.bufferNo + ". " +
                            //                         nextBufferToWrite.mers[0].ToString("X16") + " - " + nextBufferToWrite.mers[nextBufferToWrite.activeMers - 1].ToString("X16");
                            //bufferTrace.Add(writeTrace);
                            //if (bufferTrace.Count == 10)
                            //    bufferTrace.RemoveAt(0);
                            if (nextBufferToWrite.mers[0] < lastMerWritten)
                            {
                                Console.WriteLine("first mer of next buffer less than last mer of previous buffer: " + nextBufferToWrite.mers[0].ToString("X16") + " < " + lastMerWritten.ToString("X16"));
                                //foreach (string traceLine in bufferTrace)
                                //    Console.WriteLine(traceLine);
                                //Console.WriteLine("source " + nextBufferToWrite.counts[0]);
                            }
                            localTimeSpentWriting.Start();
                            if (typeOfCount == TypeOfCount.ulongType)
                                WriteULongBuffer(cbtFile, minCount, nextBufferToWrite, tesselStats, sumReps);
                            if (typeOfCount == TypeOfCount.uintType)
                                WriteUIntBuffer(cbtFile, minCount, nextBufferToWrite, tesselStats);
                            localTimeSpentWriting.Stop();
                            lastMerWritten = nextBufferToWrite.mers[nextBufferToWrite.activeMers - 1];

                            // finished with this buffer
                            pendingBuffers.Remove(nextBufferNoToWrite);
                            // now ready to write the next buffer
                            nextBufferNoToWrite++;
                            // and tell the Filler in case it's waiting for this buffer
                            nextBufferToWrite.free = true;
                            freedEWH.Set();
                        }
                    }
                    else
                        // none of the pending buffers are ready to be written yet...  
                        writtenAll = true;
                    if (pendingBuffers.Count == 0)
                        writtenAll = true;
                }
            }

            timeSpentWriting = localTimeSpentWriting.Elapsed;

            //Console.WriteLine("exiting writer #" + threadNo);
        }

        private static void WriteULongBuffer(BinaryWriter cbtFile, int minCount, MerBuffer merBuffer, ThreadStats tesselStats, Dictionary<int, long> sumReps)
        {
            int nextIdx = 0;
            ulong nextMer;
            ulong nextCountPair;
            byte[] byteBuffer = merBuffer.byteBuffer;
            int byteIdx = 0;

            while (nextIdx < merBuffer.activeMers)
            {
                nextMer = merBuffer.mers[nextIdx];
                nextCountPair = merBuffer.counts[nextIdx];

                nextIdx++;
                // merge all counts for this kMer
                while (nextIdx < merBuffer.activeMers && merBuffer.mers[nextIdx] == nextMer)
                {
                    nextCountPair += merBuffer.counts[nextIdx];
                    nextIdx++;
                }

                //WriteMerToFile(lowestMer, lowestCountPair, minCount, cbtFile);
                int countAsRead = (int)(nextCountPair >> 32);
                int countRC = (int)(nextCountPair & 0x00000000FFFFFFFF);
                int sumCount = countAsRead + countRC;

                if (sumCount >= minCount)
                {
                    //cbtFile.Write(lowestMer);
                    //cbtFile.Write(countAsRead);
                    //cbtFile.Write(countRC);
                    //Console.WriteLine(kMers.ExpandMer(lowestMer, 32) + " " + countAsRead + " " + countRC);
                    // effectively the same code but tests showed this was faster
                    byteBuffer[byteIdx++] = (byte)nextMer;
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 8);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 16);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 24);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 32);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 40);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 48);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 56);
                    byteBuffer[byteIdx++] = (byte)countAsRead;
                    byteBuffer[byteIdx++] = (byte)(countAsRead >> 8);
                    byteBuffer[byteIdx++] = (byte)(countAsRead >> 16);
                    byteBuffer[byteIdx++] = (byte)(countAsRead >> 24);
                    byteBuffer[byteIdx++] = (byte)countRC;
                    byteBuffer[byteIdx++] = (byte)(countRC >> 8);
                    byteBuffer[byteIdx++] = (byte)(countRC >> 16);
                    byteBuffer[byteIdx++] = (byte)(countRC >> 24);

                    if (byteIdx >= MerBuffer.byteBufferlen)
                    {
                        cbtFile.Write(byteBuffer, 0, byteIdx);
                        byteIdx = 0;
                    }

                    tesselStats.noOfDistinctMers++;
                    tesselStats.noOfMersWritten += sumCount;
                }
                else
                {
                    tesselStats.noOfDistinctMersDropped++;
                    tesselStats.noOfMersDropped += sumCount;
                }

                if (sumReps != null)
                {
                    if (sumReps.ContainsKey(sumCount))
                        sumReps[sumCount]++;
                    else
                        sumReps.Add(sumCount, 1);
                }
            }

            progressMersWritten = tesselStats.noOfMersWritten;

            if (byteIdx > 0)
                cbtFile.Write(byteBuffer, 0, byteIdx);
        }

        // in-lined/expanded for performance
        //private static void WriteMerToFile(ulong mer, ulong count, int minCount, BinaryWriter cbtFile)
        //{
        //    int countAsRead = (int)(count >> 32);
        //    int countRC = (int)(count & 0x00000000FFFFFFFF);
        //    int sumCount = countAsRead + countRC;

        //    if (sumCount >= minCount)
        //    {
        //        cbtFile.Write(mer);
        //        cbtFile.Write(countAsRead);
        //        cbtFile.Write(countRC);
        //        progressMersWritten += sumCount;
        //        totalUniqueMers++;
        //        totalMersWritten += sumCount;
        //    }
        //    else
        //    {
        //        progressMersDropped += sumCount;
        //        totalUniqueMersDropped++;
        //        totalMersDropped += sumCount;
        //    }

        //    if (sumReps.ContainsKey(sumCount))
        //        sumReps[sumCount]++;
        //    else
        //        sumReps.Add(sumCount, 1);
        //}

        private static void WriteUIntBuffer(BinaryWriter cbtFile, int minCount, MerBuffer merBuffer, ThreadStats tesselStats)
        {
            int nextIdx = 0;
            ulong nextMer;
            ulong nextCountPair;
            byte[] byteBuffer = merBuffer.byteBuffer;
            int byteIdx = 0;

            while (nextIdx < merBuffer.activeMers)
            {
                nextMer = merBuffer.mers[nextIdx];
                nextCountPair = merBuffer.counts[nextIdx];

                nextIdx++;
                // merge all counts for this kMer
                while (nextIdx < merBuffer.activeMers && merBuffer.mers[nextIdx] == nextMer)
                {
                    nextCountPair += merBuffer.counts[nextIdx];
                    nextIdx++;
                }

                int count = (int)(nextCountPair & 0x00000000FFFFFFFF);

                if (count >= minCount)
                {
                    //cbtFile.Write(lowestMer);
                    //cbtFile.Write(countAsRead);
                    //cbtFile.Write(countRC);
                    //Console.WriteLine(kMers.ExpandMer(lowestMer, 32) + " " + countAsRead + " " + countRC);
                    // effectively the same code but tests showed this was faster
                    byteBuffer[byteIdx++] = (byte)nextMer;
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 8);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 16);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 24);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 32);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 40);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 48);
                    byteBuffer[byteIdx++] = (byte)(nextMer >> 56);
                    byteBuffer[byteIdx++] = (byte)count;
                    byteBuffer[byteIdx++] = (byte)(count >> 8);
                    byteBuffer[byteIdx++] = (byte)(count >> 16);
                    byteBuffer[byteIdx++] = (byte)(count >> 24);

                    if (byteIdx >= MerBuffer.byteBufferlen)
                    {
                        cbtFile.Write(byteBuffer, 0, byteIdx);
                        byteIdx = 0;
                    }

                    tesselStats.noOfDistinctPairs++;
                    tesselStats.noOfPairsWritten += count;
                }
            }

            progressPairsWritten = tesselStats.noOfPairsWritten;

            if (byteIdx > 0)
                cbtFile.Write(byteBuffer, 0, byteIdx);
        }

        public abstract class MerSource
        {
            public int sourceNo;                // source no. 
            public ulong lowestMer;             // current mer available from this source (next to be consumed)
            public ulong countPair;             // and corresponding count pair
            public bool opened;                 // has this source been opened yet?
            public bool closed;                 // has this been exhausted (and closed)?
            public ulong firstMer;              // lowest mer in this source
            public ulong lastMer;               // highest mer in this source
            public string sourceTypeString;     // used for tracing
            public int partition;               // singletons partition no.

            protected const ulong asReadSingleCount = 0x0000000100000000;
            protected const ulong rcSingleCount =     0x0000000000000001;

            public abstract bool MoveToNextMer();
            public abstract void Open();
        }

        public class MerDictionarySource : MerSource
        {
            //MerDictionary.Entry[] entries;
            ulong[] mers;
            long[] counts;
            int idx = 0;
            int maxIdx;
            MerDictionary source;

            public MerDictionarySource(MerDictionary source)
            {
                this.source = source;
                opened = false;
                closed = false;
                sourceTypeString = "MerDictionary";
            }

            public override void Open()
            {
                //entries = source.entries;
                mers = source.keys;
                counts = source.values;
                idx = 0;
                maxIdx = source.Count - 1;

                if (maxIdx > 0)
                {
                    //lowestMer = entries[0].key;
                    //countPair = (ulong)entries[0].value;
                    //firstMer = entries[0].key;
                    //lastMer = entries[maxIdx].key;
                    lowestMer = mers[0];
                    countPair = (ulong)counts[0];
                    firstMer = mers[0];
                    lastMer = mers[maxIdx];

                    opened = true;
                    closed = false;
                }
                else
                {
                    opened = false;
                    closed = true;
                }


                //Console.WriteLine("opened mer dictionary. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override bool MoveToNextMer()
            {
                if (closed)
                    return false;

                if (idx == maxIdx)
                {
                    lowestMer = 0;
                    countPair = 0;
                    opened = false;
                    closed = true;
                    //Console.WriteLine("closed MerDictionary");
                    return false;
                }

                idx++;
                //lowestMer = entries[idx].key;
                //countPair = (ulong)entries[idx].value;
                lowestMer = mers[idx];
                countPair = (ulong)counts[idx];
                //count--;
                return true;
            }
        }

        // unflushed singletons
        public class MerSingletonsSource : MerSource
        {
            //MerCollection.Entry[] entries;
            long[] mers;
            bool[] merIsRC;
            int idx = 0;
            int maxIdx;
            SingletonCollection source;

            public MerSingletonsSource(SingletonCollection source, int partition)
            {
                this.source = source;
                this.partition = partition;
                opened = false;
                closed = false;
                sourceTypeString = "Singletons";
            }

            public override void Open()
            {
                //count = source.Count;
                //entries = source.entries;
                mers = source.keys;
                merIsRC = source.rcFlags;
                idx = 0;
                maxIdx = source.Count - 1;

                if (maxIdx >= 0)
                {
                    firstMer = (ulong)mers[0];
                    lastMer = (ulong)mers[maxIdx - 1];
                    lowestMer = firstMer;
                    if (merIsRC[0])
                        countPair = rcSingleCount;
                    else
                        countPair = asReadSingleCount;

                    opened = true;
                    closed = false;
                }
                else
                {
                    // empty unflushed singletons partition
                    opened = false;
                    closed = true;
                }
                //Console.WriteLine("opened unflushed singletons. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override bool MoveToNextMer()
            {
                if (closed)
                    return false;

                if (idx == maxIdx)
                {
                    lowestMer = 0;
                    countPair = 0;
                    opened = false;
                    closed = true;
                    //Console.WriteLine("closed singletons (" + firstMer.ToString("x16") + " - " + lastMer.ToString("x16") + ")");  
                    return false;
                }

                idx++;
                //lowestMer = (ulong)entries[idx].key;
                lowestMer = (ulong)mers[idx];
                if (merIsRC[idx])
                    countPair = rcSingleCount;
                else
                    countPair = asReadSingleCount;

                return true;
            }
        }

        // flushed singletons from files
        public class MerFlushedSingletonSource : MerSource
        {
            BinaryReader flushedSingletons;
            int length;
            int mersRead = 0;
            bool lowestMerIsRC = false;
            string flushedFN;

            public MerFlushedSingletonSource(string flushedFN, int partition, ulong firstMer, ulong lastMer, int flushedSingletonsCount)
            {
                this.flushedFN = flushedFN;
                this.firstMer = firstMer;
                this.lastMer = lastMer;
                this.partition = partition;
                opened = false;
                closed = false;
                sourceTypeString = "FlushedSingletons";
            }

            public override void Open()
            {
                if (opened)
                    return;

                flushedSingletons = new BinaryReader(new FileStream(flushedFN, FileMode.Open, FileAccess.Read, FileShare.None, 100000, FileOptions.SequentialScan));
                length = flushedSingletons.ReadInt32();

                lowestMer = flushedSingletons.ReadUInt64();
                lowestMerIsRC = flushedSingletons.ReadBoolean();
                if (lowestMerIsRC)
                    countPair = rcSingleCount;
                else
                    countPair = asReadSingleCount;
                mersRead = 1;
                opened = true;
                closed = false;
                //Console.WriteLine("opened " + flushedFN + " with " + length + " singletons. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override bool MoveToNextMer()
            {
                if (closed)
                    return false;

                if (mersRead == length)
                {
                    lowestMer = 0;
                    countPair = 0;
                    opened = false;
                    closed = true;
                    flushedSingletons.Close();
                    File.Delete(flushedFN);
                    //Console.WriteLine("closed flushed singletons " + flushedFN + " (" + firstMer.ToString("x16") + " - " + lastMer.ToString("x16") + ")");
                    return false;
                }
                else
                {
                    mersRead++;
                    lowestMer = flushedSingletons.ReadUInt64();
                    lowestMerIsRC = flushedSingletons.ReadBoolean();
                    if (lowestMerIsRC)
                        countPair = rcSingleCount;
                    else
                        countPair = asReadSingleCount;

                    return true;
                }
            }
        }

        public class MerFlushedLowRepSource : MerSource
        {
            BinaryReader flushedLowReps;
            int length;
            int mersRead;
            string flushedFN;

            public MerFlushedLowRepSource(string flushedFN, ulong firstMer, ulong lastMer)
            {
                this.flushedFN = flushedFN;
                this.firstMer = firstMer;
                this.lastMer = lastMer;
                opened = false;
                closed = false;
                //count = flushedLowRepCount;
                sourceTypeString = "FlushedLowReps";
                partition = -1;
            }

            public override void Open()
            {
                if (opened)
                    return;

                flushedLowReps = new BinaryReader(new FileStream(flushedFN, FileMode.Open, FileAccess.Read, FileShare.None, 100000, FileOptions.SequentialScan));
                length = flushedLowReps.ReadInt32();
                lowestMer = flushedLowReps.ReadUInt64();
                countPair = flushedLowReps.ReadUInt64();

                opened = true;
                closed = false;
                mersRead = 1;
                //Console.WriteLine("opened " + flushedFN + " with " + length + " lowreps. lowest=" + lowestMer.ToString("x16") + " " + firstMer.ToString("x16") + " " + lastMer.ToString("x16"));
            }

            public override bool MoveToNextMer()
            {
                if (closed)
                    return false;

                if (mersRead == length)
                {
                    countPair = 0;
                    opened = false;
                    closed = true;
                    flushedLowReps.Close();
                    File.Delete(flushedFN);
                    //Console.WriteLine("closed flushed lowreps " + flushedFN);
                    return false;
                }
                else
                {
                    lowestMer = flushedLowReps.ReadUInt64();
                    countPair = flushedLowReps.ReadUInt64();
                    mersRead++;
                    return true;
                }
            }
        }

        static void WriteExpandedMers(BinaryReader cbtFile, StreamWriter results, TextFormat wantedTextFormat)
        {
            bool EOF = false;
            ulong packedMer = 0;
            int merAsReadCount = 0;
            int merRCCount = 0;

            int merSize = cbtFile.ReadInt32();

            while (!EOF)
            {
                try
                {
                    packedMer = cbtFile.ReadUInt64();
                    merAsReadCount = cbtFile.ReadInt32();
                    merRCCount = cbtFile.ReadInt32();
                }
                catch
                {
                    EOF = true;
                    cbtFile.Close();
                }
                if (EOF)
                    break;

                string key = kMers.ExpandMer(packedMer, merSize);

                if (wantedTextFormat == TextFormat.tabPairs)
                    results.WriteLine(key + "\t" + merAsReadCount + "\t" + merRCCount);

                if (wantedTextFormat == TextFormat.tabSummed)
                    results.WriteLine(key + "\t" + (merAsReadCount + merRCCount));

                if (wantedTextFormat == TextFormat.faPairs)
                {
                    results.WriteLine(">" + merAsReadCount + " " + merRCCount);
                    results.WriteLine(key);
                }

                if (wantedTextFormat == TextFormat.faSummed)
                {
                    results.WriteLine(">" + (merAsReadCount + merRCCount));
                    results.WriteLine(key);
                }
            }
        }

        private static void WriteUsage()
        {
            Console.WriteLine("usage: Tessel -k k-merLength -g genomeLength [-t #threads] [-f fasta|fastq] [-tmp tempDir] [-min minCount] [-trim nn] [-trimQual nn] [-s] [-pairs] [-nopairs] [-canonical|asread] [-text textFN] [-textFormat pairs|sum|faPairs|faSum] cbtName readsFN... or readsPattern (" + version + ")");
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

        static void RateReporter()
        {
            long lastReadsRead = 0;
            bool startedPairPhase = false;

            while (!stopMonitor)
            {
                Thread.Sleep(reportInterval);

                if (programPhase == phases.tilingMers)
                {
                    long currentReadsRead = progressReadsRead;
                    Console.WriteLine("tiled " + currentReadsRead + " reads (+" + (currentReadsRead - lastReadsRead) + ")");
                    lastReadsRead = currentReadsRead;
                }
                if (programPhase == phases.writingMers)
                {
                    Console.WriteLine("wrote " + progressMersWritten + "/" + progressMersToWrite + " " + kLength + "-mers");
                    lastReadsRead = 0;  // prep for pairs
                }
                if (programPhase == phases.tilingPairs)
                {
                    if (!startedPairPhase)
                        lastReadsRead = 0;
                    startedPairPhase = true;
                    long currentReadsRead = progressReadsRead;
                    Console.WriteLine("tiled " + currentReadsRead + " reads (+" + (currentReadsRead - lastReadsRead) + ") for pairs");
                    lastReadsRead = currentReadsRead;
                }
                if (programPhase == phases.writingPairs)
                {
                    Console.WriteLine("wrote " + progressPairsWritten + "/" + progressPairsTiled + " kMer pairs");
                }
            }
        }
    }

    public class CountingThreadParams
    {
        public int threadNumber;
        public BufferedReader readsFile;
        public int merSize;
        public MerTables merTables;
        public int trimLength;
        public int trimQual;
        public int qualOffset;
        public ThreadStats threadStats;
    }

    public class SortingThreadParams
    {
        public int threadNumber;
        public EventWaitHandle sortingEWH;
        public EventWaitHandle writingEWH;
        public Queue<MerBuffer> buffersToBeSorted;
        public Queue<MerBuffer> buffersToBeWritten;
    }

    public class WritingThreadParams
    {
        public int threadNumber;
        public EventWaitHandle writingEWH;
        public Queue<MerBuffer> buffersToBeWritten;
        public EventWaitHandle freedEWH;
        public BinaryWriter bwFile;
        public Program.TypeOfCount typeOfCount;
        public int minCount;
        public ThreadStats tesselStats;
        public Dictionary<int, long> sumReps;
    }

    public class PairCountingThreadParams
    {
        public int threadNumber;
        public BufferedReader readsFile;
        public int merSize;
        public int minReadLength;
        public int pairGap;
        public int minCount;
        public MerCollections.PairTable pairsTable;
        public kMerTable kMersTable;
        public ThreadStats threadStats;
    }

    public class SorterParams
    {
        public Queue<SorterTask> waitingTasks;
    }

    public class SorterTask
    {
        public SingletonCollection merCollection;
        public MerDictionary merDictionary;
    }

    public class ThreadStats
    {
        public long noOfReadsRead = 0;              // no. of reads read from all files
        public long noOfTiledMers = 0;              // no. of mers tiled from these reads
        public long noOfValidMers = 0;              // no. of mers that are valid (do not contain Ns)
        public long noOfNMers = 0;                  // no. of mers not written because they contained Ns
        public long noOfPoorQuals = 0;              // no. of bases with poor quals trimmed from end of reads
        public long noOfDistinctMers = 0;           // no. of distinct mers found (including generated RC forms)
        public long noOfMersWritten = 0;            // no. of mers counted
        public long noOfDistinctMersDropped = 0;    // no. of distinct mers dropped (too few reps)
        public long noOfMersDropped = 0;            // no. of mers dropped         
        public long noOfTiledPairs = 0;             // no. of pairs generated/counted
        public long noOfDistinctPairs = 0;          // no of distinct pairs found
        public long noOfPairsWritten = 0;           // no. of kMer pairs written
    } 

    public class MerBuffer
    {
        public int bufferIdx;                                      // index within buffer array - tracing use only 
        public bool finishSorter;                                  // tells sorting thread to finish when it dequeues this buffer 
        public bool finishWriter;                                  // and this tells the Writer thread to finish 
        public int bufferNo;                                       // buffers are created and written in this order 
        public bool free;                                          // free/in-use flag 
        public int activeMers;                                     // no. of active mers in this buffer
        public ulong[] mers;                                       // mers waiting sorting/writing
        public ulong[] counts;                                     // and their count pairs
        public byte[] byteBuffer;                                  // byte buffer used when writing to file
        public const int byteBufferlen = 65536;                    // and size of this buffer  

        public MerBuffer(int bufferIdx, int bufferSize)
        {
            this.bufferIdx = bufferIdx;
            this.free = true;
            this.finishSorter = false;
            this.bufferNo = 0;
            this.activeMers = 0;
            this.mers = new ulong[bufferSize];                     // will be resized upwards as needed
            this.counts = new ulong[bufferSize];
            this.byteBuffer = new byte[byteBufferlen+16];          // allow a bit more in case what's being written doesn't go evenly into byteBufferLen
        }

        public void Sort()
        {
            MerArraySorter<ulong> merSorter = new MerArraySorter<ulong>(mers, counts);
            merSorter.Sort(0, activeMers - 1);
            //for (int i = 0; i < activeMers - 2; i++)
            //    if (mers[i] > mers[i + 1])
            //    {
            //        Console.WriteLine("sort failed");
            //        break;
            //    }
        }

    } // MerBuffer
}