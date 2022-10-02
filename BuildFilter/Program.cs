using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using WorkingDogsCore;

namespace BuildFilter
{
    // Builds a filter from a set of genes/contigs/... 
    // 
    // usage: BuildFilter [-v1|-v1c|-v2c|-v2a] -k merSize [+/-lcf] [-s] [-mindepth nn] [-minlength nn] genesFN 
    //
    // The set of genes/contigs etc are expected to be in FASTA format.
    //
    // Version 1. Always generate canonical kMers. Only one kMer set. Each entry is a ulong containing the packed canonical kMer.
    //            Start of file is a single int that holds the kMer size. [000000KK]
    //
    // Version 2. Steals the upper byte of the previous kMer size int as a version number. 0 == original; 1 == version 2
    //            The second byte is stolen as an as-read(0)/canonical(1) flag, and the low byte is the kMer size once again. [VVCCCxxKK]
    //
    // The rest of the file in both V1 and V2 are just the packed kMers - end of file. No length is kept because low-complexity kMers are dropped
    // during writing so the final length is unknown. This could be remedied in V3 but is currently not an issue.
    //

    class Program
    {
        //static StreamWriter trace = new StreamWriter("buildTrace.txt");

        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                WriteUsage();
                return;
            }

            int kMerSize = 0;
            bool filterLowComplexity = true;
            int minDepth = 1;
            int minLength = 0;
            List<string> targetsList = new List<string>();
            string mersFN = null;
            bool recursiveFileSearch = false;
            int version = 2;
            bool canonical = true;

            for (int p = 0; p < args.Length; p++)
            {
                if (args[p][0] == '-' || args[p][0] == '+')
                {
                    args[p] = args[p].ToLower();

                    if (args[p] == "-h" || args[p] == "-help")
                    {
                        WriteUsage();
                        return;
                    }

                    if (args[p] == "-v1" || args[p] == "-v1c")
                    {
                        version = 1;
                        canonical = true;
                        continue;
                    }

                    if (args[p] == "-v2a")
                    {
                        version = 2;
                        canonical = false;
                        continue;
                    }

                    if (args[p] == "-v2c")
                    {
                        version = 2;
                        canonical = true;
                        continue;
                    }

                    if (args[p] == "-lcf")
                    {
                        filterLowComplexity = false;
                        continue;
                    }

                    if (args[p] == "+lcf")
                    {
                        filterLowComplexity = true;
                        continue;
                    }

                    if (args[p] == "-s")
                    {
                        recursiveFileSearch = true;
                        continue;
                    }

                    if (args[p] == "-k")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -k"))
                            return;
                        try
                        {
                            kMerSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-md" || args[p] == "-mindepth")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -mindepth"))
                            return;
                        try
                        {
                            minDepth = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -mindepth parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-ml" || args[p] == "-minlength")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -minlength"))
                            return;
                        try
                        {
                            minLength = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -minlength parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    Console.WriteLine("unrecognised arg: " + args[p]);
                    return;

                } // arg starts with +/-

                targetsList.Add(args[p]);
            }

            if (targetsList == null)
            {
                Console.WriteLine("no target sequences specified");
                return;
            }

            // only one name specified - assume we'll just modify it to generate the mers file name
            if (targetsList.Count == 1)
                mersFN = targetsList[0].Substring(0, targetsList[0].LastIndexOf('.')) + "_" + kMerSize + ".mer";
            else
            {
                mersFN = targetsList[0];
                targetsList.RemoveAt(0);
            }

            List<string> targetFNs = new List<string>();
            SearchOption searchOption = recursiveFileSearch ? SearchOption.AllDirectories : SearchOption.TopDirectoryOnly;
            foreach (string targetName in targetsList)
            {
                string targetDir;
                string targetFNpart;
                SeqFiles.SplitFileName(targetName, out targetDir, out targetFNpart);
                string[] filesFromTarget = Directory.GetFiles(targetDir, targetFNpart, searchOption);
                foreach (string FN in filesFromTarget)
                    targetFNs.Add(FN);
            }

            if (targetFNs.Count == 0)
            {
                Console.WriteLine("no files found");
                return;
            }

            if (filterLowComplexity)
                Console.WriteLine("discarding low-complexity k-mers");
            else
                Console.WriteLine("retaining low-complexity k-mers");
            if (minDepth > 1)
                Console.WriteLine("discarding kMers found less than " + minDepth + " times");

            long cumulativeLength = 0;
            foreach (string FN in targetFNs)
            {
                FileInfo fi = new FileInfo(FN);
                long fileLength = fi.Length;
                cumulativeLength += fileLength;
            }

            MerDictionary<int> distinctMers = new MerDictionary<int>(2*cumulativeLength, kMerSize);

            bool EOF = false;
            ulong[] merSet = new ulong[2000];
            bool[] merValid = new bool[2000];
            DateTime start = DateTime.Now;
            int countSequences = 0;
            int tooShortCount = 0;

            foreach (string targetFN in targetFNs)
            {
                int fileFormat = SeqFiles.DetermineFileFormat(targetFN);
                StreamReader targets = new StreamReader(targetFN);

                while (!EOF)
                {
                    string seq = SeqFiles.ReadRead(targets, fileFormat);
                    if (seq == null)
                        break;
                    countSequences++;

                    if (seq.Length < minLength)
                    {
                        tooShortCount++;
                        continue;
                    }

                    //if (nextGene == "TCTCAAAGATTAAGCCATGCATGTCTAAGTACATGCCGTATTAAGGTGAAACCGCGAATGGCTCATTAAATCAGTTACGGTTCATTAGAACTTGAGCTAACTTACATGGATAACTGTGGTAATTCTAGAGCTAATACATGCACAAAAGCTTTGACCAAGCCGTGCTCGTCGCGGTGAAGGAAAAAGCGCATTTATTAGACCAAGACCAATGGGAATCATTGGGCTTTGTCTCGGTTGCCGTCAAAACGTAACCGGGCAAAGATTCCAGTTCCCTCAAACATTATGGTGACTCTAGATAACTGTAGCTAATCGCATGGCCAATGAGCCGGCGATAGATCTTTCAAGCGTCTGCCTTATCAACTGTCGATGGTAGGTTATGCGCCTACCATGGTTTTAACGGGTGACGAGGAATCAGGGTTCGATTTCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCACGCAAATTACCCAATGCCAGAACGGCGAGGTAGTGACGAAAAATAACAATACGGGACTCTAATGAGGCCCCGTAATTGGAATGAGAACAATCTAAATCCTTTAACGAGGATCTATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATCTCAGTTTTTAGTTGTTTGCGGTCCACTAACAAGTGGTTACTGCTCAACTGAACTCAACAATAATACCGATTTATTTTGGGTGGTTCTGTTGCCGGTTAGTTGACGGCGGTTTGGCCGCGATCAGTGGATTTCACCTCGTGTGGGCTGCTGGTCGTGTCCAACTGTCTCTCATGCTAGCCCAACGGCCATTCAAATGCTCATGGTGCTCTTAACCGGGTGTCATGCGGCGATCGGTACGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCGTTGAAAACGGATTTAAAACGTCTAAATTGCCCAGAATAATGTTGCATGGAATAATAAAATATGACCTCGGTTCTATTTTGTTGGTCTTTAGAACTTATTACCAATTAAGAGGTAATGATTAAGAGGGACAGACGGGGGCATTCGTATTGCGGCGCTAGAGGTGAAATTCTTGGACCGTCGCAAGACGAACTAAAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGTTCGAAGGCGATCAGATACCGCCCTAGTTCTAACCATAAACGATGCCAACCAGCAATCCGTCTGAGTTCCTTAAATGACTCGACGGGCGGCTTCCGGGAAACCAAAGTTTTTCGGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGGGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTTACCTGGCCCGGACACTAAAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTAGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTAGCCTACTAAATAGACTATGTTGGCTTGTGTAGAAATTTAACCACTAGCCATTTCATGTCGTTCGGRGCTCGTTTCGGCCTGTGTCGTTYKGGGCATCGGSARGTCGGTTCTGCCGGCTTGTCRRTGTTCTTGCGGCAYRGGTTTCGGAGCGGGTTTTCGGCGGCRTGATTTACTAGTGGCGTTTCAATAACGCCAACAGTGCTTCTTAGAGGGACAGGCGGCGATTCAGCCGCACGAAACAGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCAGGGCCGCACGCGCGCCACACTGAAGTGATCAGCGTGCTATTTTTGATGTCCTGCTCTGTTAAGAGTAGGTAACCCAATCAACCTTCTTCGTGATTGGGATAGGGGATTGTAATTATTCCCCTTGAACGAGGAATTCCCAGTAAGCGCGAGTCATAAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGATTTAGTGAGGTCTTCAGACCAGCCGATGCAGTTGTTTACTTGTTGAACAACCGCGTCTGTTGTTGGAAAGATGCCCAAACTTGATTATTTAGAGGAAGTAAAAGTCGTAA")
                    //    Debugger.Break();

                    if (kMers.GeneContainsAmbiguousBases(seq))
                        seq = kMers.CollapseAmbiguousBases(seq);

                    int mersInSeq = kMers.GenerateMersFromRead(seq, kMerSize, ref merSet, ref merValid);

                    for (int m = 0; m < mersInSeq; m++)
                    {
                        if (merValid[m])
                        {
                            ulong mer = merSet[m];
                            if (canonical)
                            {
                                ulong rcMer = kMers.ReverseComplement(mer, kMerSize);
                                if (rcMer < mer)
                                    mer = rcMer;
                            }

                            if (distinctMers.ContainsKey(mer))
                                distinctMers[mer]++;
                            else
                                distinctMers.Add(mer, 1);
                        }
                    } // all mers in ref sequence
                }

                targets.Close();
            }

            Console.WriteLine("Loaded " + distinctMers.Count + " distinct " + (canonical ? "canonical " : "as-read ") + kMerSize + "-mers from " + countSequences + " sequences in " +
                    (DateTime.Now - start).TotalSeconds.ToString("#.0") + "s. " + (tooShortCount > 0 ? (tooShortCount + " too short.") : ""));


            BinaryWriter tiledSeqs = new BinaryWriter(File.Open(mersFN, FileMode.Create, FileAccess.Write));
            int mersWritten = 0;
            int skippedLC = 0;
            int skippedLowCount = 0;

            if (version == 1)
                tiledSeqs.Write(kMerSize);                         // first int is just kMer size with zero first byte
            if (version == 2)
            {
                uint v2Header = 0x01000000 | (uint)(canonical ? 1 : 0) << 16 | (uint)kMerSize;
                tiledSeqs.Write(v2Header);
            }

            foreach (KeyValuePair<ulong, int> kvp in distinctMers)
            {
                ulong mer = kvp.Key;
                int count = kvp.Value;

                //string kMer = kMers.ExpandMer(mer, merSize);
                //if (mer == 0x25c70321ce000000 || kMer == "AAAAAAAACCATTTTTTAAC")
                //    Debugger.Break();

                // remove any low complexity k-mers from the filter on the fly
                if (filterLowComplexity)
                {
                    // just don't write out low complexity mers... 
                    if (LowComplexityMer(kMerSize, mer) || LowComplexityMer(kMerSize, kMers.ReverseComplement(mer, kMerSize)))
                    {
                        skippedLC++;
                        //trace.WriteLine("- " + kMers.ExpandMer(mer, merSize));
                        continue;
                    }
                }

                if (count < minDepth)
                {
                    skippedLowCount++;
                    continue;
                }

                tiledSeqs.Write(mer);                        // ulong - packed mer - max 32 bases
                mersWritten++;
                //trace.WriteLine("+ " + kMers.ExpandMer(mer, merSize));
            }

            tiledSeqs.Close();
            //trace.Close();

            Console.WriteLine("Wrote " + mersWritten + " " + kMerSize + "-mers (skipped " + skippedLC + " low complexity, skipped " + skippedLowCount + " low reps) in " +
                    (DateTime.Now - start).TotalSeconds.ToString("#.0") + "s");

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

        private static void WriteUsage()
        {
            Console.WriteLine("usage: BuildFilter [-v1|-v2]  -k merSize [+/-lcf] [-mindepth nn] [-minlength nn] [-s] [filterFN] genesFN");
        }

        private static bool LowComplexityMer(int merSize, ulong mer)
        {
            // check if any pairs or triplets of bases is dominant. 

            Dictionary<ulong, int> distinctSet = new Dictionary<ulong, int>(merSize);
            if (CheckForDiversity(merSize, mer, distinctSet, 2))
                return true;
            if (CheckForDiversity(merSize, mer, distinctSet, 3))
                return true;

            return false;

        }

        private static bool CheckForDiversity(int merSize, ulong mer, Dictionary<ulong, int> distinctSet, int ptqSize)
        {
            // 'ptqSize' is the size of the repeat targets - pairs, triplets or quads
            ulong ptqMask = 0xffffffffffffffff << (64 - (ptqSize * 2));     // mask for pair/triplet/quad at top ok
            bool lacksDiversity = false;
            distinctSet.Clear();
            //bool trace = true;

            // number of complete pairs/triplets in the kMer (overlapping)
            int ptqInMer = merSize - ptqSize + 1;
            ulong tempMer = mer;

            // generate all 'ptq's and save their counts
            for (int b = 0; b < ptqInMer; b++)
            {
                ulong topBits = tempMer & ptqMask;
                tempMer = tempMer << 2;
                if (distinctSet.ContainsKey(topBits))
                    distinctSet[topBits]++;
                else
                    distinctSet.Add(topBits, 1);
            }

            int maxCount = 0;
            ulong maxKey = 0;
            foreach (KeyValuePair<ulong, int> kvp in distinctSet)
            {
                int count = kvp.Value;
                if (count > maxCount)
                {
                    maxKey = kvp.Key;
                    maxCount = count;
                }
            }
            int secondCount = 0;
            foreach (KeyValuePair<ulong, int> kvp in distinctSet)
            {
                int count = kvp.Value;
                if (kvp.Key != maxKey && count > secondCount)
                    secondCount = count;
            }

            int topTwoCount = maxCount + secondCount;
            int cutOff = ptqInMer * 55/100;
            if (secondCount == 1)
                cutOff = ptqInMer / 2;
            if (topTwoCount > cutOff)
                lacksDiversity = true;

            //if (trace && lacksDiversity)
            //{
            //    Console.WriteLine(kMers.ExpandMer(mer, merSize));
            //    foreach (KeyValuePair<ulong, int> kvp in distinctSet)
            //        Console.WriteLine(kMers.ExpandMer(kvp.Key, ptqSize) + "=" + kvp.Value);
            //    Console.WriteLine();
            //}

            return lacksDiversity;
        }
    }
}


