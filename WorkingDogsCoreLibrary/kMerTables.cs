using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;

namespace WorkingDogsCore
{
    // a collection of canonical kMers + count pairs loaded from a .cbt file
    // ---------------------------------------------------------------------
    public class kMerTable
    {
        const int bytesPerCBTMer = 8 * 2;

        public MerTable<ulong> kMersTable = null;

        public long distinctMersLoaded = 0;
        public long totalMersLoaded = 0;
        public int averageDepthLoaded = 0;

        public int merSize = 0;

        List<ulong> hdubMerList;
        public HashSet<ulong> hdubFilter;

        public kMerTable(string cbtFN, int minLoadDepth, int deepThreshold)
        {
            Console.WriteLine("Loading kMers from " + cbtFN + " (min=" + minLoadDepth + ")");
            DateTime loadStart = DateTime.Now;

            FileInfo tileFI = new FileInfo(cbtFN);
            long merArraySize = 0;
            long cbtFileLength = tileFI.Length;
            merArraySize = cbtFileLength / bytesPerCBTMer;        // file size in kmers

            if (deepThreshold > 0)
                hdubMerList = new List<ulong>(100);
            else
                hdubMerList = null;

            BinaryReader cbtFile = null;
            cbtFile = new BinaryReader(new FileStream(cbtFN, FileMode.Open, FileAccess.Read, FileShare.Read, 1000000, FileOptions.SequentialScan));

            merSize = cbtFile.ReadInt32();
            // allocate the kMer table
            kMersTable = new MerTable<ulong>(merArraySize, merSize);

            LoadCBTTiles(cbtFile, minLoadDepth, deepThreshold);

            if (merSize <= 0 || merSize > 32)
                return;

            cbtFile.Close();

            averageDepthLoaded = (int)(totalMersLoaded / distinctMersLoaded);

            hdubFilter = GenerateHDUBFilter(hdubMerList);
            hdubMerList = null;

            Console.WriteLine("Finished loading kMer table. " +
                                  distinctMersLoaded + "/" + totalMersLoaded + " " + merSize + "-mers loaded, average depth " + averageDepthLoaded +
                                  " in " + (DateTime.Now - loadStart).TotalSeconds.ToString("#.0") + "s");
        }

        private void LoadCBTTiles(BinaryReader cbtFile, int minLoadDepth, int minDeepDepth)
        {
            ulong kMer = 0;                                         // packed form of next read kMer
            ulong countPair = 0;                                    // a packed count pair for this kMer
            int asReadCount = 0;                                    // as-read count for the next read mer
            int rcCount = 0;                                        // RC count for the next read mer
            int sumCount = 0;                                       // sum of as-read and rc counts

            List<ulong> hdubMers = new List<ulong>();               // accumulating set of high-depth, unbalanced kMers
            List<int> hdubCounts = new List<int>();                 // and their summed counts
            int deepestDepthSoFar = 0;                              // the depest depth in the above lists
            int deepThreshold = minDeepDepth;                       // and the threshold for admission to the 'hdub' list

            bool EOF = false;

            while (!EOF)
            {
                try
                {
                    kMer = cbtFile.ReadUInt64();
                    asReadCount = cbtFile.ReadInt32();
                    rcCount = cbtFile.ReadInt32();
                    countPair = ((ulong)asReadCount << 32) + (ulong)rcCount;
                    sumCount = asReadCount + rcCount;
                }
                catch
                {
                    EOF = true;
                }

                if (EOF)
                    break;

                // add all mers tiled from the reads into our hash table
                if (sumCount >= minLoadDepth)
                {
                    kMersTable.Add(kMer, countPair);
                    distinctMersLoaded++;
                    totalMersLoaded += sumCount;
                }

                // if we're accumulating hdubs...
                if (hdubMerList != null)
                { 
                    int highCount = asReadCount;
                    int lowCount = rcCount;
                    bool rcAsRead = false;

                    if (lowCount > highCount)
                    {
                        int temp = highCount;
                        highCount = lowCount;
                        lowCount = temp;
                        rcAsRead = true;
                    }

                    // does this mer look like a candidate for the HDUB list?
                    if (sumCount > deepThreshold && lowCount <= sumCount / 100)
                    {
                        // if the new 'deep' kMer is deeper than anything in the list already, adjust the threshold 
                        if (sumCount > deepestDepthSoFar)
                        {
                            deepestDepthSoFar = sumCount;
                            deepThreshold = sumCount / 4;
                        }

                        // add the as-read form of the kMer to the lists
                        hdubMers.Add(rcAsRead ? kMers.ReverseComplement(kMer, merSize) : kMer);
                        hdubCounts.Add(sumCount);

                        // trim down the list if it's getting a bit big (just want the top few HDUBs)
                        if (hdubMers.Count > 100)
                        {
                            int idx = 0;
                            while (hdubCounts[idx] < deepThreshold)
                            {
                                hdubMers.RemoveAt(idx);
                                hdubCounts.RemoveAt(idx);
                                idx++;
                                if (idx == hdubCounts.Count)
                                    break;
                            }
                        }
                    }
                }

            } // all cbt tuples in the file

            cbtFile.Close();

            // flush final buffer
            bool mersLoaded = kMersTable.LoadFinished();
            GC.Collect();

            if (!mersLoaded)
            {
                Console.WriteLine(".cbt load failed - check file is sorted");
                distinctMersLoaded = 0;
                totalMersLoaded = 0;
                kMersTable = null;
            }

            if (hdubMerList != null)
            {
                // return the HDUBs in sorted order (highest count first)
                ulong[] sortedMers = hdubMers.ToArray();
                int[] sortedCounts = hdubCounts.ToArray();
                Array.Sort<int, ulong>(sortedCounts, sortedMers);
                Array.Reverse(sortedMers);
                hdubMers.Clear();
                foreach (ulong mer in sortedMers)
                    hdubMerList.Add(mer);
            }
        }

        // turns the collected set of HDUB kMers into the HDUB filter
        private HashSet<ulong> GenerateHDUBFilter(List<ulong> hdubMers)
        {
            HashSet<ulong> hdubFilterSet = new HashSet<ulong>();

            if (hdubMers.Count == 0)
                return hdubFilterSet;

            //foreach (ulong mer in hdubMers)
            //{
            //    int plusCount;
            //    int rcCount;
            //    GetDepthPair(mer, out plusCount, out rcCount);
            //    Console.WriteLine(kMers.ExpandMer(mer, merSize) + " " + plusCount + " " + rcCount);
            //}

            List<ulong> hdubMerList = new List<ulong>();
            List<ulong> hdubMerTopFew = new List<ulong>();
            int hdubThreshhold = GetDepthSum(hdubMers[0]) / 2;

            // takes the deepest remaining kMer and extend it - capturing the extended kMers and removing them from the hdub set
            while (hdubMers.Count > 0 && GetDepthSum(hdubMers[0]) > hdubThreshhold)
            {
                ulong startingHDUBMer = hdubMers[0];
                hdubMerList.Add(startingHDUBMer);
                hdubMers.RemoveAt(0);
                ExtendHDUBMers(startingHDUBMer, hdubMers, hdubMerList);

                for (int i = 0; i < Math.Min(hdubMerList.Count, 10); i++)
                {
                    hdubMerTopFew.Add(hdubMerList[i]);
                    // add all of the HDUBs to the filter up front in an attenpt to put them together and have them cache better
                    hdubFilterSet.Add(hdubMerList[i]);
                }
            }

            // generate single-base variants of these adapter-like kMers
            foreach (ulong hdubMer in hdubMerTopFew)
            {
                if (!hdubFilterSet.Contains(hdubMer))
                    hdubFilterSet.Add(hdubMer);

                List<ulong> hdubMerVariants = new List<ulong>();
                kMers.GenerateMerSubVariants(hdubMer, hdubMerVariants, merSize);
                foreach (ulong hdubMerVariant in hdubMerVariants)
                {
                    // only need to add variants that appear in the seq file somewhere
                    int plusCount = 0;
                    int rcCount = 0;
                    GetDepthCounts(hdubMerVariant, out plusCount, out rcCount);
                    if (plusCount == 0 && rcCount == 0)
                        continue;

                    // and have one zero or near zero count
                    int lowCount = plusCount;
                    int highCount = rcCount;
                    if (lowCount > highCount)
                    {
                        lowCount = rcCount;
                        highCount = plusCount;
                    }
                    if (lowCount == 0 || lowCount * 10 < highCount)
                        if (!hdubFilterSet.Contains(hdubMerVariant))
                            hdubFilterSet.Add(hdubMerVariant);

                }
            }

            return hdubFilterSet;
        }

        // extends the starting kMer from those in the hdubMers list, adding the selected kMers to the set and removing them from the list
        // The kMers in the starting list are as-read and the extension has to be done in as-read form 
        private void ExtendHDUBMers(ulong startingHDUBMer, List<ulong> hdubMers, List<ulong> hdubMerSet)
        {
            hdubMerSet.Clear();
            bool extending = true;
            ulong nextAsReadHDUBMer = startingHDUBMer;

            //Console.WriteLine("finding followers of " + kMers.ExpandMer(startingHDUBMer, merSize) + "\t" + kMers.ExpandMer(nextAsReadHDUBMer, merSize));

            int fillShift = 64 - merSize * 2;   // unused RHS bits 
            int startShift = fillShift + 2;     // 2 bits more to remove the last base
            ulong startMask = 0xFFFFFFFFFFFFFFFF >> startShift;

            // get the starting k-1 bases of the remaining kMers
            Dictionary<ulong, ulong> startsOfMers = new Dictionary<ulong, ulong>();
            foreach (ulong merAsRead in hdubMers)
            {
                ulong merFragment = merAsRead >> startShift;
                if (!startsOfMers.ContainsKey(merFragment))
                    startsOfMers.Add(merAsRead >> startShift, merAsRead);
            }
            //foreach (KeyValuePair<ulong, ulong> kvp in startsOfMers)
            //    Console.WriteLine("start: " + kMers.ExpandMer(kvp.Key, 32) + "\tfull " + kMers.ExpandMer(kvp.Value, 32));

            while (extending)
            {
                ulong endOfMer = nextAsReadHDUBMer >> fillShift & startMask;
                //Console.WriteLine("e: " + kMers.ExpandMer(endOfMer, 32));
                if (startsOfMers.ContainsKey(endOfMer))
                {
                    //Console.WriteLine("m: " + kMers.ExpandMer(endOfMer, 32));
                    hdubMers.Remove(nextAsReadHDUBMer);
                    hdubMerSet.Add(nextAsReadHDUBMer);
                    nextAsReadHDUBMer = startsOfMers[endOfMer];
                    startsOfMers.Remove(endOfMer);
                    //Console.WriteLine("next: " + kMers.ExpandMer(nextAsReadHDUBMer, 32));
                }
                else
                {
                    break;
                }
            }
        }

        // gets total depth of a kMer that *may* be in the table
        public int GetDepthSum(ulong mer)
        {
            ulong rcMer = kMers.ReverseComplement(mer, merSize);
            if (rcMer < mer)
                mer = rcMer;

            ulong countPair = 0;
            kMersTable.TryGetValue(mer, out countPair);
            int plusCount = (int)(countPair >> 32);
            int rcCount = (int)(countPair & 0xFFFFFFFF);

            return plusCount + rcCount;
        }

        // gets total depth of a kMer that *may* be in the table (and remember if the depths are badly unbalanced, or somewhat unbalanced)
        public int GetDepthSum(ulong mer, int minDepth, out bool unbalanced, out bool tilted)
        {
            ulong rcMer = kMers.ReverseComplement(mer, merSize);
            if (rcMer < mer)
                mer = rcMer;

            ulong countPair = 0;
            kMersTable.TryGetValue(mer, out countPair);
            int plusCount = (int)(countPair >> 32);
            int rcCount = (int)(countPair & 0xFFFFFFFF);
            int lowerDepth = plusCount;
            int higherDepth = rcCount;
            if (rcCount < lowerDepth)
            {
                lowerDepth = rcCount;
                higherDepth = plusCount;
            }

            unbalanced = false;
            int summedDepth = plusCount + rcCount;
            unbalanced = lowerDepth < minDepth &&
                            summedDepth > 20 ? (lowerDepth * 10 < higherDepth) : lowerDepth <= 1;
            tilted = !unbalanced && summedDepth > minDepth && (plusCount * 2 < rcCount | rcCount * 2 < plusCount);
            return summedDepth;
        }

        // gets depths for a kMer that *may* be in the table
        public bool GetDepthCounts(ulong mer, out int plusCount, out int rcCount)
        {
            bool foundMer = true;
            ulong rcMer = kMers.ReverseComplement(mer, merSize);
            ulong countPair;
            bool rcMerWasCanonical = false;

            rcMerWasCanonical = rcMer < mer;
            if (rcMerWasCanonical)
                mer = rcMer;

            if (!kMersTable.TryGetValue(mer, out countPair))
            {
                //string missingMer = kMers.ExpandMer(packedMer);
                plusCount = 0;
                rcCount = 0;
                foundMer = false;
            }

            // extract the plus, RC and qual values from the packed ulong value
            if (rcMerWasCanonical)
            {
                rcCount = (int)(countPair >> 32);
                plusCount = (int)(countPair & 0xFFFFFFFF);
            }
            else
            {
                plusCount = (int)(countPair >> 32);
                rcCount = (int)(countPair & 0xFFFFFFFF);
            }

            return foundMer;
        }


    }

    // a collection of canonical kMers + count pairs loaded from a .cbt file
    // ---------------------------------------------------------------------
    public class PairTable
    {
        const int bytesPerPair = (64 + 32) / 8;

        public MerTable<int> pairsTable = null;
        public int pairGap = 0;
        public int pairFullLength = 0;
        public long distinctPairsLoaded = 0;
        public long totalPairsLoaded = 0;
        public int averageDepthLoaded = 0;

        public long gpdCalls = 0;

        public PairTable(string pairsFN, int minLoadDepth)
        {
            Console.WriteLine("Loading kMer pairs from " + pairsFN + " (min=" + minLoadDepth + ")");
            DateTime loadStart = DateTime.Now;

            FileInfo pairsFI = new FileInfo(pairsFN);

            long pairsFileLength = pairsFI.Length - 4;
            long pairsArrayLength = pairsFileLength / bytesPerPair + 1;

            // allocate the pairs table (32-mers - 2x16)
            pairsTable = new MerTable<int>(pairsArrayLength, 32);

            BinaryReader pairsFile = new BinaryReader(new FileStream(pairsFN, FileMode.Open, FileAccess.Read, FileShare.Read, 1000000, FileOptions.SequentialScan));
            pairGap = pairsFile.ReadInt32();

            pairFullLength = 2 * kMerPairs.pairFragmentSize + pairGap;

            bool EOF = false;

            while (!EOF)
            {
                try
                {
                    ulong pair = pairsFile.ReadUInt64();
                    int pairDepth = pairsFile.ReadInt32();

                    if (pairDepth > minLoadDepth)
                    {
                        pairsTable.Add(pair, pairDepth);
                        distinctPairsLoaded++;
                        totalPairsLoaded += pairDepth;
                    }
                }
                catch (EndOfStreamException)
                {
                    EOF = true;
                }
                if (EOF)
                    break;
            }

            pairsFile.Close();
            bool prsLoaded = pairsTable.LoadFinished();

            averageDepthLoaded = (int)(totalPairsLoaded / distinctPairsLoaded);

            if (prsLoaded)
                Console.WriteLine("Finished loading kMer pairs table. " +
                      distinctPairsLoaded + "/" + totalPairsLoaded + " " + "pairs loaded, average depth " + averageDepthLoaded +
                      " in " + (DateTime.Now - loadStart).TotalSeconds.ToString("#.0") + "s");
            else
            {
                Console.WriteLine(".prs load failed - check file is sorted");
                totalPairsLoaded = 0;
                pairsTable = null;
            }

            GC.Collect();
        }

        public int GetPairDepth(ulong pair)
        {
            ulong rcPair = kMers.ReverseComplement(pair, kMerPairs.pairSize);

            if (rcPair < pair)
                pair = rcPair;

            int pairDepth = 0;
            pairsTable.TryGetValue(pair, out pairDepth);

            return pairDepth;
        }
    }
}

