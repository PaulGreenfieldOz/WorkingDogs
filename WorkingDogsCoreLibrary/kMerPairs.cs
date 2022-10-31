using System;
using System.Collections.Generic;
using System.Text;

namespace WorkingDogsCore
{
    public class kMerPairs
    {
        public const int pairFragmentSize = 16;               // in bases
        public const ulong firstFragmentMask = 0xFFFFFFFFFFFFFFFF << (64 - pairFragmentSize * 2);
        public const ulong lastFragmentMask = 0xFFFFFFFFFFFFFFFF >> (64 - pairFragmentSize * 2);
        public const int pairSize = 2 * pairFragmentSize;     // in bases (will always be 32)

        public static bool ConstructPair(Sequence seq, int m, int pairGap, out ulong pair)
        {
            pair = 0;

            if (m + kMerPairs.pairFragmentSize * 2 + pairGap > seq.Length)
                return false;

            ulong firstFragment;
            bool firstFragmentOK = Sequence.CondenseMer(seq, m, pairFragmentSize, out firstFragment);
            ulong lastFragment;
            bool lastFragmentOK = Sequence.CondenseMer(seq, (m + pairFragmentSize + pairGap), pairFragmentSize, out lastFragment);

            pair = firstFragment | (lastFragment >> (64 - pairFragmentSize * 2));

            return firstFragmentOK && lastFragmentOK;
        }

        // generate a pair going forward into a read/contig. The current kMer is at the start of the pair region, m --> first base of starting kMer
        // XXXXXXXXXXXXXXXXxxxxxxxx----------------yyyyyyyyYYYYYYYYYYYYYYYY (24-mers & pair region length = 64)
        // ^m              gggggggggggggggggggggggggggggggg
        // --> XXXXXXXXXXXXXXXXYYYYYYYYYYYYYYYY
        public static bool ConstructPairFwd(Sequence seq, int m, int pairGap, out ulong pair)
        {
            pair = 0;

            // return if there's not enough room for a pair region
            if (m + kMerPairs.pairFragmentSize * 2 + pairGap > seq.Length)
                return false;

            ulong firstFragment;
            bool firstFragmentOK = Sequence.CondenseMer(seq, m, pairFragmentSize, out firstFragment);
            ulong lastFragment;
            bool lastFragmentOK = Sequence.CondenseMer(seq, (m + pairFragmentSize + pairGap), pairFragmentSize, out lastFragment);

            pair = firstFragment | (lastFragment >> (64 - pairFragmentSize * 2));

            return firstFragmentOK && lastFragmentOK;
        }

        // generate a pair going backwards into a read/contig. The current kMer is at the end of the pair region, m --> first base of last kMer
        //                 gggggggggggggggggggggggggggggggg
        // XXXXXXXXXXXXXXXXxxxxxxxx----------------yyyyyyyyYYYYYYYYYYYYYYYY (24-mers & pair length = 72)
        //                                         ^m
        // --> XXXXXXXXXXXXXXXXYYYYYYYYYYYYYYYY
        public static bool ConstructPairRvs(Sequence seq, int m, int kMerSize, int pairGap, out ulong pair)
        {
            pair = 0;

            int startOfPairRegion = (m + kMerSize) - (pairFragmentSize * 2 + pairGap);
            // return if there's not enough room for a pair region
            if (startOfPairRegion < 0)
                return false;

            ulong firstFragment;
            bool firstFragmentOK = Sequence.CondenseMer(seq, startOfPairRegion, pairFragmentSize, out firstFragment);
            ulong lastFragment;
            bool lastFragmentOK = Sequence.CondenseMer(seq, startOfPairRegion+pairFragmentSize+pairGap, pairFragmentSize, out lastFragment);


            pair = firstFragment | (lastFragment >> (64 - pairFragmentSize * 2));

            return firstFragmentOK && lastFragmentOK;
        }

        public static bool ConstructPair(Sequence seq, int m, int merSize, int pairGap, out ulong pair, out ulong firstMer, out ulong secondMer)
        {
            pair = 0;
            firstMer = 0;
            secondMer = 0;

            if (m + kMerPairs.pairFragmentSize * 2 + pairGap > seq.Length)
                return false;

            bool firstMerOK = Sequence.CondenseMer(seq, m, merSize, out firstMer);
            bool secondMerOK = Sequence.CondenseMer(seq, (m + 2 * kMerPairs.pairFragmentSize + pairGap - merSize), merSize, out secondMer);

            pair = (firstMer & kMerPairs.firstFragmentMask) | ((secondMer >> (64 - merSize * 2)) & kMerPairs.lastFragmentMask);

            return firstMerOK && secondMerOK;
        }

        public static bool ConstructPairIncremental(Sequence seq, int m, int pairGap, ulong previousPair, out ulong pair)
        {
            pair = 0;

            if (m + kMerPairs.pairFragmentSize * 2 + pairGap > seq.Length)
                return false;

            bool firstFragmentOK = true;
            ulong firstFragment = (previousPair & kMerPairs.firstFragmentMask) << 2;
            long newBase = kMers.BaseCharToInt(seq.Bases[m + kMerPairs.pairFragmentSize - 1]);
            if (newBase < 0)
                firstFragmentOK = false;
            firstFragment = firstFragment | (ulong)newBase << kMerPairs.pairFragmentSize * 2;

            bool lastFragmentOK = true;
            ulong lastFragment = (previousPair << 2) & kMerPairs.lastFragmentMask;
            newBase = kMers.BaseCharToInt(seq.Bases[m + kMerPairs.pairFragmentSize * 2 + pairGap - 1]);
            if (newBase < 0)
                lastFragmentOK = false;
            lastFragment = lastFragment | (ulong)newBase;

            pair = firstFragment | lastFragment;

            return firstFragmentOK && lastFragmentOK;
        }


        public static int GeneratePairsFromRead(Sequence read, int pairGap, ref ulong[] pairs, ref bool[] pairValid)
        {
            int pairsInRead = read.Length - (kMerPairs.pairSize + pairGap) + 1;
            bool pairIsValid = false;
            ulong pair = 0;

            if (pairsInRead < 1)
                return 0;

            if (pairs.Length < pairsInRead)
            {
                Array.Resize<ulong>(ref pairs, pairsInRead + 100);
                if (pairValid != null)
                    Array.Resize<bool>(ref pairValid, pairsInRead + 100);
            }

            for (int i = 0; i < pairsInRead; i++)
            {
                if (pairIsValid)
                {
                    pairIsValid = kMerPairs.ConstructPairIncremental(read, i, pairGap, pair, out pair);
                    if (pairValid != null)
                        pairValid[i] = pairIsValid;
                    pairs[i] = pair;
                }
                else
                {
                    pairIsValid = kMerPairs.ConstructPair(read, i, pairGap, out pair);
                    if (pairValid != null)
                        pairValid[i] = pairIsValid;
                    pairs[i] = pair;
                }
            }

            return pairsInRead;
        }

        public static int GeneratePairsFromRead(Sequence read, int pairGap, ref ulong[] pairs)
        {
            bool[] pairValid = null;
            return GeneratePairsFromRead(read, pairGap, ref pairs, ref pairValid);
        }

        public static int GeneratePairsFromRead(string read, int pairGap, ref ulong[] pairs)
        {
            bool[] pairValid = null;
            Sequence seq = new Sequence(read);
            return GeneratePairsFromRead(seq, pairGap, ref pairs, ref pairValid);
        }

        public static int GeneratePairsFromRead(string read, int pairGap, ref ulong[] pairs, ref bool[] pairsValid)
        {
            Sequence seq = new Sequence(read);
            return GeneratePairsFromRead(seq, pairGap, ref pairs, ref pairsValid);
        }
    }
}
