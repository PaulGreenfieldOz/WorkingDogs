using System;
using System.Collections.Generic;
using System.IO;

namespace WorkingDogsCore
{

    public static class TrimExtend
    {
        const int contigSmallRegion = 20;

        // Trim any high depth, unbalanced kMers - likely to be adapters. These can be at either end of a read.
        public static int TrimUnbalancedHighDepth(kMerTable kMersTable,  Sequence read, Sequence quals, int minDepth, ulong[] mers, out int startingBase, out int trimmedLength)
        {
            int merSize = kMersTable.merSize;
            HashSet<ulong> hdubFilter = kMersTable.hdubFilter;

            startingBase = 0;
            trimmedLength = read.Length;

            int merIdx = 0;
            int plusCount = 0;
            int rcCount = 0;
            bool startsWithDeepUnbalanced = false;
            int startingLength = read.Length;
            int mersInRead = read.Length - merSize + 1;

            for (int i = 0; i < Math.Min(mersInRead, merSize + 1); i++)
            {
                if (hdubFilter.Contains(mers[i]))
                {
                    startsWithDeepUnbalanced = true;
                    break;
                }
            }

            // does the read start with a suspect kMer?
            if (startsWithDeepUnbalanced)
            {
                // Skip until we find a normal looking kMer 
                while (merIdx < mersInRead)
                {
                    kMersTable.GetDepthCounts(mers[merIdx], out plusCount, out rcCount);

                    // found a possibly reasonable kMer?
                    if (plusCount >= minDepth && rcCount >= minDepth && !hdubFilter.Contains(mers[merIdx]))
                        break;

                    merIdx++;
                }

                startingBase = merIdx;
                trimmedLength = read.Length - merIdx;
            }

            // now have found a sensible-looking part of the read, so search for any (other) suspect kMers and end-trim if we find any.
            // this will be the normal trimming path for reads where the adapters appear at the end of the read
            if (trimmedLength > merSize)
            {
                while (merIdx < mersInRead)
                {
                    //FindMerInUniqueMers(mers[merIdx], out plusCount, out rcCount);
                    //if (DeepUnbalanced(plusCount, rcCount))
                    //{
                    //    trimmedLength = trimmedLength - (mersInRead - merIdx - 1 + merSize);
                    //    break;
                    //}
                    if (hdubFilter.Contains(mers[merIdx]))
                    {
                        trimmedLength = trimmedLength - (mersInRead - merIdx - 1 + merSize);
                        break;
                    }
                    merIdx++;
                }
            }

            return startingLength - trimmedLength;
        }

        public static int TrimReadStart(Sequence read, Sequence quals, int mersInRead, int pairsInRead, bool contigs, 
                                        int[] depths, int minDepth, int[] pairDepths, int minPairDepth, 
                                        int startingBase, int readLength)
        {
            int firstMerIdx = startingBase;
            int basesTrimmed = 0;

            // trim not-good-looking bases from the start of a read
            // always check for at least two consecutive good bases before deciding trimming is done
            while (firstMerIdx < mersInRead-1)
            {
                // if trimming contigs, recalculate the minimum depth each time to better handle the coverage depth droop that can occur at the end of contigs
                if (contigs)
                {
                    int lengthForMean = contigSmallRegion;
                    if (firstMerIdx + contigSmallRegion > mersInRead)
                        lengthForMean = mersInRead - firstMerIdx;
                    if (lengthForMean != 0)
                    {
                        minDepth = HarmonicMean(depths, firstMerIdx, lengthForMean) / 3;
                        minPairDepth = HarmonicMean(pairDepths, firstMerIdx, lengthForMean) / 3;
                    }
                }

                bool foundPoorMer = false;

                // check that the kMer depth is satisfactory
                int depth = depths[firstMerIdx];
                if (depths[firstMerIdx] < minDepth || depths[firstMerIdx+1] < minDepth)
                    foundPoorMer = true;

                // check that we have a satisfactory pair (if we can)
                if (!foundPoorMer && firstMerIdx < pairsInRead-1)
                {
                    if (pairDepths[firstMerIdx] < minPairDepth || pairDepths[firstMerIdx+1] < minPairDepth)
                        foundPoorMer = true;
                }

                if (foundPoorMer)
                {
                    firstMerIdx++;
                    basesTrimmed++;
                }
                else
                {
                    break;
                }

            } // scan for a good starting kMer

            return basesTrimmed;
        }

        public static int TrimReadEnd(Sequence read, Sequence quals, int mersInRead, int pairsInRead, bool contigs,
                                      ulong[] mers, int[] depths, int minDepth, int[] pairDepths, int minPairDepth, 
                                      int startingBase, int readLength, int merSize, kMerTable kMersTable)
        {
            int lastMerIdx = -1;
            int basesTrimmed = 0;
            int mersToPairs = mersInRead - pairsInRead;
            bool lowComplexityFound = false;

            // find the first bad kMer in the read (could be in the middle rather than at the end)
            for (int mi = startingBase; mi < mersInRead; mi++)
            {
                if (!lowComplexityFound)
                    lowComplexityFound = kMers.LowComplexity(mers[mi], merSize);
                if (depths[mi] == 0)
                {
                    lastMerIdx = mi;
                    break;
                }
            }
            // or the first bad pair if we didn't find a bad kMer
            if (lastMerIdx < 0)
                for (int pi = startingBase; pi < pairsInRead; pi++)
                {
                    if (pairDepths[pi] == 0)
                    {
                        // convert pairdIdx to merIdx
                        lastMerIdx = pi + mersToPairs;
                        break;
                    }
                }

            // didn't find a bad depth/pair so don't trim anything
            if (lastMerIdx < 0)
                return 0;

            // found a bad depth/pair go backwards looking for a good kMer
            basesTrimmed = mersInRead - lastMerIdx;

            while (lastMerIdx > startingBase)
            {
                // if trimming contigs, recalculate the minimum depth each time to better handle the coverage depth droop that can occur at the ends of contigs
                if (contigs)
                {
                    int startForMean = lastMerIdx - contigSmallRegion;
                    if (startForMean >= 0)
                        minDepth = HarmonicMean(depths, startForMean, contigSmallRegion) / 3;
                }

                bool trimThisMer = false;

                // trim kMer if it isn't deep enough
                if (depths[lastMerIdx] < minDepth)
                    trimThisMer = true;

                // trim even harder if we found a low-complexity kMer. These will typically have one zero depth and we'll trim back to the last balanced kMer
                if (!trimThisMer && lowComplexityFound)
                {
                    int plusCount = 0;
                    int rcCount = 0;

                    kMersTable.GetDepthCounts(mers[lastMerIdx], out plusCount, out rcCount);

                    if (plusCount * 10 < rcCount | rcCount * 10 < plusCount)
                        trimThisMer = true;
                }

                if (!trimThisMer)
                {
                    int pairIdx = lastMerIdx - mersToPairs;
                    if (pairIdx > startingBase)
                    {
                        int pairDepth = pairDepths[pairIdx];
                        if (pairDepth <= minPairDepth)
                            trimThisMer = true;
                    }
                }

                if (trimThisMer)
                {
                    lastMerIdx--;
                    basesTrimmed++;
                    continue;
                }

                // stop scanning once we found a deep-enough kMer with a good pair (or no pair was available)
                break;
            }

            return basesTrimmed;
        }

        public static int ExtendRead(Sequence read, Sequence quals, int mersInRead, int pairsInRead, 
                                     ulong[] mers, int[] depths, int minDepth, int[] pairDepths, int minPairDepth,
                                     kMerTable kMersTable, PairTable pairsTable, int wantedExtension)
        {
            int merSize = kMersTable.merSize;
            int pairGap = 0;
            if (pairsInRead > 0)
                pairGap = pairsTable.pairGap;
            int lastMerIdx = mersInRead - 1;
            int merFill = 64 - merSize * 2;                         // unused bits at RHS of ulong kMer
            int basesAdded = 0;
            char[] bases = new char[] { 'A', 'C', 'G', 'T' };

            while (basesAdded < wantedExtension)
            {
                ulong lastMer = mers[lastMerIdx];
                ulong rshiftedMer = (lastMer >> merFill) << 2;      // right-shifted kMer with a hole at the last base
                int bestBase = -1;                                  // assume no good base found
                int bestBasesFound = 0;                             // how many acceptable extensions were found - only alllowed answer is '1'
                ulong bestMer = 0;                                  // the kMer that was accepted
                int bestDepth = 0;                                  // and its depth
                int bestPairDepth = 0;                              // and its pair depth

                bool lastGoodMerUnbalanced = false;
                int plusCount = 0;
                int rcCount = 0;
                int pairDepth = 0;

                kMersTable.GetDepthCounts(lastMer, out plusCount, out rcCount);
                lastGoodMerUnbalanced = (plusCount + rcCount > minDepth) && (plusCount <= 1 || rcCount <= 1);

                // try each of the possible bases
                for (ulong b = 0; b < 4; b++)
                {
                    ulong possibleMer = (rshiftedMer | b) << merFill;
                    kMersTable.GetDepthCounts(possibleMer, out plusCount, out rcCount);

                    // skip alternative if not deep enough
                    if (plusCount + rcCount < minDepth)
                        continue;

                    // and be very wary of kMers with 0/1 on one strand
                    if (!lastGoodMerUnbalanced && (plusCount <= 1 || rcCount <= 1))
                        continue;

                    // and also skip if it doesn't have a pair match
                    if (pairsInRead > 0)
                    {
                        int firstMerIdx = lastMerIdx + 1 + merSize - pairGap - 2 * kMerPairs.pairFragmentSize;
                        if (firstMerIdx >= 0)
                        {
                            ulong pair;
                            bool pairValid = kMerPairs.ConstructPair(read, firstMerIdx, pairGap, out pair);

                            pairDepth = pairsTable.GetPairDepth(pair);
                            if (pairDepth < minPairDepth)
                                continue;
                        }
                    }

                    // have an alternative that is both deep enough and has a pair
                    bestBase = (int)b;
                    bestMer = possibleMer;
                    bestDepth = plusCount + rcCount;
                    bestPairDepth = pairDepth;
                    bestBasesFound++;
                }

                // and if we found just one acceptable extension, accept it and move on
                if (bestBasesFound == 1)
                {
                    read.Append(bases[bestBase]);
                    if (quals.Length > 0)
                        quals.Append(quals.Bases[quals.Length - 1]);    // copy the qual from the previous base
                    lastMerIdx++;
                    mers[lastMerIdx] = bestMer;
                    depths[lastMerIdx] = bestDepth;
                    pairDepths[lastMerIdx] = bestPairDepth;
                    basesAdded++;
                }
                else
                    break;
            }

            return basesAdded;
        }

        private static int FindPoorMer(kMerTable kMersTable, PairTable pairsTable, Sequence read, ulong[] mers, int[] depths)
        {
            int merSize = kMersTable.merSize;
            int pairGap = pairsTable.pairGap;

            int firstMerIdx = 0;
            int mersInRead = read.Length - merSize + 1;
            int firstPoorMer = -1;
            int minDepthForRead = HarmonicMean(depths, 0, mersInRead) / 3;
            if (minDepthForRead == 0)
                minDepthForRead = kMersTable.averageDepthLoaded / 10;

            while (firstMerIdx < mersInRead)
            {
                // check every kMer in the trimmed/extended read
                int plusCount = 0;
                int rcCount = 0;
                kMersTable.GetDepthCounts(mers[firstMerIdx], out plusCount, out rcCount);

                if (plusCount + rcCount < minDepthForRead)
                {
                    firstPoorMer = firstMerIdx;
                    break;
                }

                // finally check that we have a pair (if we can)
                if (pairsTable != null)
                {
                    int lastMerIdx = firstMerIdx + pairGap + 2 * kMerPairs.pairFragmentSize - merSize;
                    if (lastMerIdx < mersInRead)
                    {
                        ulong firstMerInPair = mers[firstMerIdx];
                        ulong lastMerInPair = mers[lastMerIdx];
                        ulong pair = (firstMerInPair & kMerPairs.firstFragmentMask) | ((lastMerInPair >> (64 - merSize * 2)) & kMerPairs.lastFragmentMask);
                        ulong rcPair = kMers.ReverseComplement(pair, kMerPairs.pairSize);

                        if (rcPair < pair)
                            pair = rcPair;

                        int pairDepth = 0;
                        pairDepth = pairsTable.GetPairDepth(pair);
                        if (pairDepth < minDepthForRead)
                        {
                            firstPoorMer = firstMerIdx;
                            break;
                        }
                    }
                }

                firstMerIdx++;


            } // scan for a good starting kMer

            return firstPoorMer;
        }

        private static int HarmonicMean(int[] depths, int startingMer, int merCount)
        {
            int noOKMers = 0;                       // no. of mers that are better than OK depth
            int HM = 0;

            // calculate the average depth for the passable mers
            double invSum = 0.0f;
            for (int m = 0; m < merCount; m++)
            {
                int depth = depths[startingMer + m];
                if (depth > 0)
                {
                    noOKMers++;
                    invSum += 1.0f / (double)depth;
                }
            }

            if (noOKMers > 0)
                HM = (int)((double)noOKMers / invSum);      // harmonic mean of just the 'normal' k-mers

            return HM;
        }

        private static int MinDepth(int[] depths, int startingMer, int merCount)
        {
            int minDepth = int.MaxValue;

            for (int m = 0; m < merCount; m++)
            {
                int depth = depths[startingMer + m];
                if (depth < minDepth && depth > 0)
                    minDepth = depth;
            }

            return minDepth;
        }
    }
}