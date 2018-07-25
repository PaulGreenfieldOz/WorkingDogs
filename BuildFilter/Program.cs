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
    // usage: BuildFilter -k merSize [+/-lcf] [-s] genesFN 
    //
    // The set of genes/contigs etc are expected to be in FASTA format.

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

            int merSize = 0;
            bool filterLowComplexity = true;
            int min = 1;
            List<string> targetsList = new List<string>();
            string mersFN = null;
            bool recursiveFileSearch = false;

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
                            merSize = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -k parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }

                    if (args[p] == "-m" || args[p] == "-min")
                    {
                        if (!CheckForParamValue(p, args.Length, "number expected after -min"))
                            return;
                        try
                        {
                            min = Convert.ToInt32(args[p + 1]);
                        }
                        catch
                        {
                            Console.WriteLine("expected a number for the -min parameter: " + args[p + 1]);
                            return;
                        }
                        p++;
                        continue;
                    }


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
                mersFN = targetsList[0].Substring(0, targetsList[0].LastIndexOf('.')) + "_" + merSize + ".mer";
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
            if (min > 1)
                Console.WriteLine("discarding kMers found less than " + min + " times");

            long cumulativeLength = 0;
            foreach (string FN in targetFNs)
            {
                FileInfo fi = new FileInfo(FN);
                long fileLength = fi.Length;
                cumulativeLength += fileLength;
            }

            MerDictionary<int> allFilterMers = new MerDictionary<int>(cumulativeLength/8, merSize);

            bool EOF = false;
            ulong[] merSet = new ulong[2000];
            bool[] merValid = new bool[2000];
            DateTime start = DateTime.Now;
            int countSequences = 0;

            foreach (string targetFN in targetFNs)
            {
                int fileFormat = SeqFiles.DetermineFileFormat(targetFN);
                StreamReader targets = new StreamReader(targetFN);

                while (!EOF)
                {
                    string nextGene = SeqFiles.ReadRead(targets, fileFormat);
                    if (nextGene == null)
                        break;
                    countSequences++;

                    //if (nextGene == "TCTCAAAGATTAAGCCATGCATGTCTAAGTACATGCCGTATTAAGGTGAAACCGCGAATGGCTCATTAAATCAGTTACGGTTCATTAGAACTTGAGCTAACTTACATGGATAACTGTGGTAATTCTAGAGCTAATACATGCACAAAAGCTTTGACCAAGCCGTGCTCGTCGCGGTGAAGGAAAAAGCGCATTTATTAGACCAAGACCAATGGGAATCATTGGGCTTTGTCTCGGTTGCCGTCAAAACGTAACCGGGCAAAGATTCCAGTTCCCTCAAACATTATGGTGACTCTAGATAACTGTAGCTAATCGCATGGCCAATGAGCCGGCGATAGATCTTTCAAGCGTCTGCCTTATCAACTGTCGATGGTAGGTTATGCGCCTACCATGGTTTTAACGGGTGACGAGGAATCAGGGTTCGATTTCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCACGCAAATTACCCAATGCCAGAACGGCGAGGTAGTGACGAAAAATAACAATACGGGACTCTAATGAGGCCCCGTAATTGGAATGAGAACAATCTAAATCCTTTAACGAGGATCTATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGATCTCAGTTTTTAGTTGTTTGCGGTCCACTAACAAGTGGTTACTGCTCAACTGAACTCAACAATAATACCGATTTATTTTGGGTGGTTCTGTTGCCGGTTAGTTGACGGCGGTTTGGCCGCGATCAGTGGATTTCACCTCGTGTGGGCTGCTGGTCGTGTCCAACTGTCTCTCATGCTAGCCCAACGGCCATTCAAATGCTCATGGTGCTCTTAACCGGGTGTCATGCGGCGATCGGTACGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCGTTGAAAACGGATTTAAAACGTCTAAATTGCCCAGAATAATGTTGCATGGAATAATAAAATATGACCTCGGTTCTATTTTGTTGGTCTTTAGAACTTATTACCAATTAAGAGGTAATGATTAAGAGGGACAGACGGGGGCATTCGTATTGCGGCGCTAGAGGTGAAATTCTTGGACCGTCGCAAGACGAACTAAAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTTAGAGGTTCGAAGGCGATCAGATACCGCCCTAGTTCTAACCATAAACGATGCCAACCAGCAATCCGTCTGAGTTCCTTAAATGACTCGACGGGCGGCTTCCGGGAAACCAAAGTTTTTCGGTTCCGGGGGAAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGGGCCTGCGGCTTAATTTGACTCAACACGGGAAAACTTACCTGGCCCGGACACTAAAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTAGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTAGCCTACTAAATAGACTATGTTGGCTTGTGTAGAAATTTAACCACTAGCCATTTCATGTCGTTCGGRGCTCGTTTCGGCCTGTGTCGTTYKGGGCATCGGSARGTCGGTTCTGCCGGCTTGTCRRTGTTCTTGCGGCAYRGGTTTCGGAGCGGGTTTTCGGCGGCRTGATTTACTAGTGGCGTTTCAATAACGCCAACAGTGCTTCTTAGAGGGACAGGCGGCGATTCAGCCGCACGAAACAGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCAGGGCCGCACGCGCGCCACACTGAAGTGATCAGCGTGCTATTTTTGATGTCCTGCTCTGTTAAGAGTAGGTAACCCAATCAACCTTCTTCGTGATTGGGATAGGGGATTGTAATTATTCCCCTTGAACGAGGAATTCCCAGTAAGCGCGAGTCATAAGCTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGAATGATTTAGTGAGGTCTTCAGACCAGCCGATGCAGTTGTTTACTTGTTGAACAACCGCGTCTGTTGTTGGAAAGATGCCCAAACTTGATTATTTAGAGGAAGTAAAAGTCGTAA")
                    //    Debugger.Break();

                    if (kMers.GeneContainsAmbiguousBases(nextGene))
                        nextGene = kMers.CollapseAmbiguousBases(nextGene);

                    int mersInRefGene = kMers.GenerateMersFromRead(nextGene, merSize, ref merSet, ref merValid);

                    for (int m = 0; m < mersInRefGene; m++)
                    {
                        if (merValid[m])
                        {
                            ulong mer = merSet[m];
                            // use canonical form to reduce hashSet size (may only have small effect due to use of strand-specific primers)
                            ulong rcMer = kMers.ReverseComplement(mer, merSize);
                            if (rcMer < mer)
                                mer = rcMer;

                            if (allFilterMers.ContainsKey(mer))
                                allFilterMers[mer]++;
                            else
                                // add the canonical mer to the hashSet as we know it's not already there
                                allFilterMers.Add(mer, 1);
                        }
                    } // all mers in ref sequence
                }

                targets.Close();
            }

            Console.WriteLine("Loaded " + allFilterMers.Count + " distinct " + merSize + "-mers from " + countSequences + " sequences in " +
                    (DateTime.Now - start).TotalSeconds.ToString("#.0") + "s");


            BinaryWriter tiledGenes = new BinaryWriter(File.Open(mersFN, FileMode.Create, FileAccess.Write));
            int mersWritten = 0;
            int skippedLC = 0;
            int skippedLowCount = 0;

            tiledGenes.Write(merSize);                         // tag file with mer size for safety

            foreach (KeyValuePair<ulong, int> kvp in allFilterMers)
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
                    if (LowComplexityMer(merSize, mer) || LowComplexityMer(merSize, kMers.ReverseComplement(mer, merSize)))
                    {
                        skippedLC++;
                        //trace.WriteLine("- " + kMers.ExpandMer(mer, merSize));
                        continue;
                    }
                }

                if (count < min)
                {
                    skippedLowCount++;
                    continue;
                }

                tiledGenes.Write(mer);                        // ulong - packed mer - max 32 bases
                mersWritten++;
                //trace.WriteLine("+ " + kMers.ExpandMer(mer, merSize));
            }

            tiledGenes.Close();
            //trace.Close();

            Console.WriteLine("Wrote " + mersWritten + " " + merSize + "-mers (skipped " + skippedLC + " low complexity, skipped " + skippedLowCount + " low reps)  in " +
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
            Console.WriteLine("usage: BuildFilter -k merSize [+/-lcf] [-min min] [-s] [filterFN] genesFN");
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


