using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace WorkingDogsCore
{
    //
    // A collection of methods used to create and manipulate k-mers (k<=32).
    //
    // k-mers are normally handled in 'packed' form, with each base taking 2 bits, and the bases are left-adjusted within a ulong (64 bits).
    //
    public static class kMers
    {
        // Mapping low 5 bits of ASCII chars to 2-bit base encoding (A=00, C=01, G=10, T=11). Same mapping works for lower case bases.
        // Used by BaseCharToInt(char baseChar) to convert char bases to 2-bit longs
        private static long[] baseToInt = new long[] { -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
        //                                              @, A,  B, C,  D,  E,  F, G,  H,  I,  J,  K,  L,  M,  N,  O,  P,  Q,  R,  S, T,  U,  V,  W,  X,  Y,  Z,  [, \\,  ],  ^,  _
        // mapping 2-bit bases back to char form
        public static char[] baseToChar = new char[] { 'A', 'C', 'G', 'T' };
        //public static char[] baseToChar = new char[] { 'a', 'c', 'g', 't' };

        // Mapping low 5 bits of ASCII chars to ACGT or not. Used to detect degenerate bases. Summing over all translated bases will produce 0 if there are no ambiguous bases
        private static int[] baseToNotACGT = new int[] { 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        //                                               @, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z, [,\\, ], ^, _

        // Translate bases to ACGT - collapsing degenerate bases. Indexed by low 5 bits of char. Always map to first base of set, default is A
        private static char[] baseToACGT = new char[] { '@', 'A', 'B', 'C', 'G', 'E', 'F', 'G', 'A', 'I', 'J', 'G', 'L', 'A', 'A', 'O', 'P', 'Q', 'G', 'G', 'T', 'U', 'G', 'A', 'X', 'T', 'Z', '[', '\\', ']', '^', '_' };
        //                                               @,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,   [,  \\,   ],   ^,   _

        // Complement bases - including degenerate bases. Indexed by low 5 bits of char. Always map to first base in following table. Unknown are mapped to N
        public static char[] baseToComplement = new char[] { '@', 'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', 'N', 'N', 'Y', 'S', 'A', 'N', 'B', 'W', 'N', 'R', 'N', '[', '\\', ']', '^', '_' };
        //                                                    @,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,   [,  \\,   ],   ^,   _

        // valid base codes (including degenerate bases)
        public static HashSet<char> validBaseCodes = new HashSet<char> { 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N', 'I' };

        //  b   meaning         reason                                  rcm     rcb
        //
        //  G   G               Guanine                                 C       C
        //  A   A               Adenine                                 T       T
        //  T   T               Thymine                                 A       A
        //  C   C               Cytosine                                G       G
        //  R   G or A          puRine                                  C T     Y
        //  Y   T or C          pYrimidine                              A G     R
        //  M   A or C          aMino                                   T G     K
        //  K   G or T          Keto                                    C A     M
        //  S   G or C          Strong interaction(3 H bonds)           C G     S
        //  W   A or T          Weak interaction(2 H bonds)             T A     W
        //  H   A or C or T     not - G, H follows G in the alphabet    T G A   D
        //  B   G or T or C     not - A, B follows A                    C A G   V
        //  V   G or C or A     not-T(not - U), V follows U             C G T   B
        //  D   G or A or T     not-C, D follows C                      C T A   H
        //  N   any             aNy                                     N       N
        //  I   any             Inosine (a match-anything base)         N       N

        public static string ExpandMer(ulong packedMer, int merSize)
        {
            char[] charMer = new char[merSize];             // char mer string under construction
            ulong tempMer = packedMer;                      // copy of packedMer for shifting
            ulong intBase = 0;                              // current base in binary form

            for (int i = 0; i < merSize; i++)
            {
                intBase = (ulong)tempMer >> 62;
                tempMer = (ulong)(tempMer << 2);

                charMer[i] = baseToChar[intBase];
            }

            return new string(charMer, 0, merSize);
        }

        public static ulong CondenseMer(string mer, int merSize)
        {
            ulong packedMer;
            CondenseMer(mer, merSize, out packedMer);
            return packedMer;
        }

        public static ulong CondenseMer(string mer)
        {
            ulong packedMer;
            int merSize = mer.Length;
            CondenseMer(mer, merSize, out packedMer);
            return packedMer;
        }

        // Convert a text-form k-mer into its 2-bit binary form. 
        // The function will return false if any non-ACGT bases are found. Lower case acgt are also accepted.
        // 
        public static bool CondenseMer(string seq, int merSize, out ulong packedMer)
        {
            return CondenseMer(seq, 0, merSize, out packedMer);
        }

        public static bool CondenseMer(string seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq[m + start];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return true;
        }

        public static bool CondenseMer(StringBuilder seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq[m + start];
                long packedBase = BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedMer = 0;
                    return false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return true;
        }


        public static bool CondenseMerIncremental(int merSize, ulong previousMer, string line, int m, out ulong nextMer)
        {
            char nextBase = line[m + merSize - 1];
            long packedBase = BaseCharToInt(nextBase);
            nextMer = (previousMer << 2) | (ulong)packedBase << (64 - merSize * 2);
            if (packedBase < 0)
            {
                nextMer = 0;
                return false;
            }
            return true;
        }

        public static bool CondenseMerIncremental(int merSize, ulong previousMer, char nextBase, out ulong nextMer)
        {
            long packedBase = BaseCharToInt(nextBase);
            nextMer = (previousMer << 2) | (ulong)packedBase << (64 - merSize * 2);
            if (packedBase < 0)
            {
                nextMer = 0;
                return false;
            }
            return true;
        }

        public static long BaseCharToInt(char baseChar)
        {
            return baseToInt[baseChar & 0x1f];
        }

        public static string CollapseAmbiguousBases(string seq)
        {
            char[] seqChars = seq.ToCharArray();
            for (int i = 0; i < seq.Length; i++)
                seqChars[i] = baseToACGT[seq[i] & 0x1f];

            return new string(seqChars);
        }

        public static bool GeneContainsAmbiguousBases(string seq)
        {
            int nonACGTCount = 0;
            for (int i = 0; i < seq.Length; i++)
                nonACGTCount += baseToNotACGT[seq[i] & 0x1f];

            return nonACGTCount > 0;
        }

        // Takes a sequence (such as a primer) that can contain degenerate bases, and returns all possible expansions of this seq.
        // Degenerate bases are expanded, and variants containing up to subsWanted substitution differences are generated. 
        public static void GenerateSeqVariants(string seq, int variedLength, char[] degenerateBases, Dictionary<char, List<char>> degenerateBaseExpansions, HashSet<string> variantsSet, int subsWanted)
        {
            GenerateSeqVariants(seq, 0, variedLength, degenerateBases, degenerateBaseExpansions, variantsSet, subsWanted);
        }

        public static void GenerateSeqVariants(string seq, int startIdx, int variedLength, char[] degenerateBases, Dictionary<char, List<char>> degenerateBaseExpansions, HashSet<string> variantsSet, int subsWanted)
        {
            // trim seq down to common length if necessary (never used)
            //if (seq.Length > seqLength)
            //    seq = seq.Substring(0, seqLength);

            HashSet<string> expandedVariantsSet = new HashSet<string>();

            // does the seq contain any degenerate bases?
            if (seq.IndexOfAny(degenerateBases) >= 0)
            {
                List<string> expandedSeqs = new List<string>();
                expandedSeqs.Add(seq);
                foreach (KeyValuePair<char, List<char>> expansion in degenerateBaseExpansions)
                    expandedSeqs = ExpandDegenerateBase(expandedSeqs, expansion.Key, expansion.Value);
                foreach (string expandedSeq in expandedSeqs)
                    expandedVariantsSet.Add(expandedSeq);
            }
            else
                expandedVariantsSet.Add(seq);

            while (subsWanted > 0)
            {
                HashSet<string> startingVariants = new HashSet<string>(expandedVariantsSet);
                expandedVariantsSet.Clear();

                foreach (string startingVariant in startingVariants)
                {
                    if (!expandedVariantsSet.Contains(startingVariant))
                        expandedVariantsSet.Add(startingVariant);

                    List<ulong> seqVariantsPacked = new List<ulong>();
                    ulong packedVariant = kMers.CondenseMer(startingVariant);
                    kMers.GenerateMerSubVariants(packedVariant, seqVariantsPacked, startIdx, variedLength);
                    foreach (ulong seqVariantPacked in seqVariantsPacked)
                    {
                        string seqVariant = kMers.ExpandMer(seqVariantPacked, seq.Length);
                        if (!expandedVariantsSet.Contains(seqVariant))
                            expandedVariantsSet.Add(seqVariant);
                    }
                }

                subsWanted--;
            }

            variantsSet.UnionWith(expandedVariantsSet);
        }

        // simple GenerateSeqVariants (without degenerate bases)
        public static void GenerateSeqVariants(string seq, int variedLength, HashSet<string> variantsSet, int subsWanted)
        {
            // trim seq down to specified length if necessary (never used)
            //if (seq.Length > seqLength)
            //    seq = seq.Substring(0, seqLength);

            HashSet<string> expandedVariantsSet = new HashSet<string>();

            expandedVariantsSet.Add(seq);

            while (subsWanted > 0)
            {
                HashSet<string> startingVariants = new HashSet<string>(expandedVariantsSet);
                expandedVariantsSet.Clear();

                foreach (string startingVariant in startingVariants)
                {
                    if (!expandedVariantsSet.Contains(startingVariant))
                        expandedVariantsSet.Add(startingVariant);

                    List<ulong> seqVariantsPacked = new List<ulong>();
                    ulong packedVariant = kMers.CondenseMer(startingVariant);
                    kMers.GenerateMerSubVariants(packedVariant, seqVariantsPacked, variedLength);
                    foreach (ulong seqVariantPacked in seqVariantsPacked)
                    {
                        string seqVariant = kMers.ExpandMer(seqVariantPacked, seq.Length);
                        if (!expandedVariantsSet.Contains(seqVariant))
                            expandedVariantsSet.Add(seqVariant);
                    }
                }

                subsWanted--;
            }

            variantsSet.UnionWith(expandedVariantsSet);
        }

        public static void InitialiseDegenerateBaseTables(out char[] degenerateBases, out Dictionary<char, List<char>> expansions)
        {
            //  G       G               Guanine
            //  A       A               Adenine
            //  T       T               Thymine
            //  C       C               Cytosine
            //  R       G or A          puRine
            //  Y       T or C          pYrimidine
            //  M       A or C          aMino
            //  K       G or T          Keto
            //  S       G or C          Strong interaction(3 H bonds)
            //  W       A or T          Weak interaction (2 H bonds) 
            //  H       A or C or T     not - G, H follows G in the alphabet
            //  B       G or T or C     not - A, B follows A
            //  V       G or C or A     not-T(not - U), V follows U
            //  D       G or A or T     not-C, D follows C
            //  N       any             aNy
            //  I       any             Inosine

            degenerateBases = new char[] { 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N', 'I' };
            expansions = new Dictionary<char, List<char>>();

            expansions.Add('R', new List<char>());
            expansions['R'].Add('G');
            expansions['R'].Add('A');
            expansions.Add('Y', new List<char>());
            expansions['Y'].Add('T');
            expansions['Y'].Add('C');
            expansions.Add('M', new List<char>());
            expansions['M'].Add('A');
            expansions['M'].Add('C');
            expansions.Add('K', new List<char>());
            expansions['K'].Add('G');
            expansions['K'].Add('T');
            expansions.Add('S', new List<char>());
            expansions['S'].Add('G');
            expansions['S'].Add('C');
            expansions.Add('W', new List<char>());
            expansions['W'].Add('A');
            expansions['W'].Add('T');
            expansions.Add('H', new List<char>());
            expansions['H'].Add('A');
            expansions['H'].Add('C');
            expansions['H'].Add('T');
            expansions.Add('B', new List<char>());
            expansions['B'].Add('G');
            expansions['B'].Add('T');
            expansions['B'].Add('C');
            expansions.Add('V', new List<char>());
            expansions['V'].Add('G');
            expansions['V'].Add('C');
            expansions['V'].Add('A');
            expansions.Add('D', new List<char>());
            expansions['D'].Add('G');
            expansions['D'].Add('A');
            expansions['D'].Add('T');
            expansions.Add('N', new List<char>());
            expansions['N'].Add('A');
            expansions['N'].Add('C');
            expansions['N'].Add('G');
            expansions['N'].Add('T');
            expansions.Add('I', new List<char>());
            expansions['I'].Add('A');
            expansions['I'].Add('C');
            expansions['I'].Add('G');
            expansions['I'].Add('T');
        }

        // takes a list of seqs and returns the same list with one single degenerate base expanded (in any seqs in which it occurs)
        private static List<string> ExpandDegenerateBase(List<string> seqs, char degenerateBase, List<char> baseExpansion)
        {
            List<string> expandedSeqs = new List<string>();

            foreach (string seq in seqs)
            {
                // does this barcode contain this degenerate base?
                if (seq.IndexOf(degenerateBase) >= 0)
                {
                    List<string> seqWithExpandedBase = ExpandBase(seq, degenerateBase, baseExpansion);
                    foreach (string expandedSeq in seqWithExpandedBase)
                        expandedSeqs.Add(expandedSeq);
                }
                else
                    expandedSeqs.Add(seq);
            }

            return expandedSeqs;
        }

        // takes a single sequence containing an degenerate base and returns a list of sequences with this base replaced
        // As the degenerate base can occur multiple times, we loop over the base replacement code until no more replacements are made
        private static List<string> ExpandBase(string seq, char degenerateBase, List<char> baseExpansion)
        {
            List<string> startingSeqs = new List<string>();
            List<string> expandedSeqs = new List<string>();
            bool baseReplaced = true;

            startingSeqs.Add(seq);

            while (baseReplaced)
            {
                baseReplaced = false;
                foreach (string nextSeq in startingSeqs)
                {
                    int abIdx = nextSeq.IndexOf(degenerateBase);
                    if (abIdx >= 0)
                    {
                        StringBuilder nextSeqSB = new StringBuilder(nextSeq);
                        foreach (char replacementBase in baseExpansion)
                        {
                            nextSeqSB[abIdx] = replacementBase;
                            string seqWithReplacementBase = nextSeqSB.ToString();
                            expandedSeqs.Add(seqWithReplacementBase);
                            baseReplaced = true;
                        }
                    }
                    else
                        expandedSeqs.Add(nextSeq);
                }

                if (baseReplaced)
                {
                    startingSeqs = expandedSeqs;
                    expandedSeqs = new List<string>();
                }
            }

            return expandedSeqs;
        }

        public static string ReverseComplement(string s)
        {
            return ReverseComplement(s, new char[s.Length]);   
        }

        public static string ReverseComplement(string s, char[] rcs)
        {
            int sl = s.Length;
            if (rcs.Length < sl)
                Array.Resize<char>(ref rcs, sl);

            for (int i = 0; i < sl; i++)
                rcs[sl - i - 1] = baseToComplement[s[i] & 0x1f];

            return new string(rcs);
        }

        public static string ReverseComplement(string s, int startIdx, int length, char[] rcs)
        {
            if (rcs.Length < length)
                Array.Resize<char>(ref rcs, length);

            for (int i = startIdx; i < startIdx+length; i++)
                rcs[length - i - 1] = baseToComplement[s[i] & 0x1f];

            return new string(rcs);
        }

        public static string ReverseComplement(string s, int startIdx, int length)
        {
            return ReverseComplement(s, startIdx, length, new char[length]);
        }

        public static StringBuilder ReverseComplement(StringBuilder s)
        {
            int sl = s.Length;
            StringBuilder rcs = new StringBuilder(sl);
            rcs.Length = sl;

            for (int i = 0; i < sl; i++)
                rcs[sl - i - 1] = baseToComplement[s[i] & 0x1f];

            return rcs;
        }

        // fast reverse complement for a packed k-mer 
        public static ulong ReverseComplement(ulong mer, int merSize)
        {
            mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
            mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
            mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
            mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
            mer = (mer >> 32) | (mer << 32);  // reversed and right-shifted after this line. Need to complement then left-shift.

            if (merSize < 32)
                return ~mer << (64 - merSize * 2);
            else
                return ~mer;
        }

        public static ulong Canonical(ulong mer, int merSize)
        {
            ulong topBase = mer >> 62;
            ulong bottomBaseComplemented = (~mer >> (32 - merSize) * 2) & 0x0000000000000003;
            if (topBase < bottomBaseComplemented)
                return mer;

            ulong startingMer = mer;
            mer = ((mer >> 2) & 0x3333333333333333) | ((mer & 0x3333333333333333) << 2);
            mer = ((mer >> 4) & 0x0F0F0F0F0F0F0F0F) | ((mer & 0x0F0F0F0F0F0F0F0F) << 4);
            mer = ((mer >> 8) & 0x00FF00FF00FF00FF) | ((mer & 0x00FF00FF00FF00FF) << 8);
            mer = ((mer >> 16) & 0x0000FFFF0000FFFF) | ((mer & 0x0000FFFF0000FFFF) << 16);
            mer = (mer >> 32) | (mer << 32);  // reversed and right-shifted after this line. Need to complement then left-shift.

            if (merSize < 32)
                mer = ~mer << (64 - merSize * 2);
            else
                mer = ~mer;

            if (startingMer < mer)
                return startingMer;
            else
                return mer;
        }

        public static int GenerateMerSubVariants(ulong mer, List<ulong> merVariants, int variedLength)
        {
            return GenerateMerSubVariants(mer, merVariants, 0, variedLength);
        }

        public static int GenerateMerSubVariants(ulong mer, List<ulong> merVariants, int startIdx, int variedLength)
        {
            int variantsAdded = 0;

            ulong baseMask = 0xc000000000000000;

            for (int m = startIdx; m < startIdx+variedLength; m++)
            {
                ulong merWithHole = mer & ~(baseMask >> (m * 2));
                for (ulong b = 0; b <= 3; b++)
                {
                    ulong newBase = b << (64 - (m + 1) * 2);
                    ulong merVariant = merWithHole | newBase;
                    if (merVariant == mer)
                        continue;
                    merVariants.Add(merVariant);
                    variantsAdded++;
                }
            }
            return variantsAdded;
        }

        public static int GenerateMerLastSubVariants(ulong mer, List<ulong> merVariants, int merSize)
        {
            int variantsAdded = 0;

            //ulong baseMask = 0xffffffffffffffff << (32 - merSize + 1) * 2;
            ulong baseMask = 0xffffffffffffffff << (merSize + 1) * 2;
            ulong merWithHole = mer & baseMask;
            for (ulong b = 0; b <= 3; b++)
            {
                ulong newBase = b << (64 - merSize * 2);
                ulong merVariant = merWithHole | newBase;
                if (merVariant == mer)
                    continue;
                merVariants.Add(merVariant);
                variantsAdded++;
            }

            return variantsAdded;
        }

        public static bool LowComplexity(ulong mer, int merSize)
        {
            ulong previousTopBasePair = 0xffffffffffffffff;
            ulong topBasePair;
            ulong shiftedMer = mer;
            int hpBases = 0;
            bool inHPRun = false;

            for (int b = 0; b < merSize; b++)
            {
                topBasePair = shiftedMer & 0xf000000000000000;
                if (topBasePair == previousTopBasePair)
                    if (inHPRun)
                        hpBases++;
                    else
                    {
                        inHPRun = true;
                        hpBases += 2;
                    }
                previousTopBasePair = topBasePair;
                shiftedMer = shiftedMer << 2;
            }

            return hpBases > 6;
        }

        public static int GenerateMersFromRead(string read, int merSize, ref ulong[] merSet, ref bool[] merValid)
        {
            int readLength = read.Length;
            int mersInRead = readLength - merSize + 1;
            bool merIsValid = false;
            ulong lastMer = 0;

            if (mersInRead < 1)
                return 0;

            if (merSet.Length < mersInRead)
            {
                Array.Resize<ulong>(ref merSet, mersInRead + 100);
                Array.Resize<bool>(ref merValid, mersInRead + 100);
            }

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                {
                    merIsValid = CondenseMerIncremental(merSize, lastMer, read, i, out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
                else
                {
                    merIsValid = CondenseMer(read.Substring(i, merSize), merSize, out lastMer);
                    merValid[i] = merIsValid;
                    merSet[i] = lastMer;
                }
            }

            return mersInRead;
        }

        public static int GenerateExpandedMersFromRead(string read, int kMerSize, List<ulong> kMerSet, char[] degenerateBases, Dictionary<char, List<char>> degenerateBaseExpansions)
        {
            int readLength = read.Length;
            int kMersInRead = readLength - kMerSize + 1;
            bool kMerWasNormal = false;
            ulong previousMer = 0;
            HashSet<String> expandedMers = new HashSet<string>();

            if (kMersInRead < 1)
                return 0;

            for (int i = 0; i < kMersInRead; i++)
            {
                // normal, fast path
                if (kMerWasNormal)
                {
                    kMerWasNormal = CondenseMerIncremental(kMerSize, previousMer, read, i, out ulong currentMer);

                    if (kMerWasNormal)
                    {
                        kMerSet.Add(currentMer);
                        previousMer = currentMer;
                    }
                }

                // kmers with degenerate bases (and first kMer in seq as well)
                if (!kMerWasNormal)
                {
                    string kMer = read.Substring(i, kMerSize);
                    expandedMers.Clear();
                    ulong currentMer = 0;
                    // expand any degenerate bases in the kMer (may be none)
                    GenerateSeqVariants(kMer, kMerSize, degenerateBases, degenerateBaseExpansions, expandedMers, 0);
                    foreach (string kMerExpansion in expandedMers)
                    {
                        CondenseMer(kMerExpansion, kMerSize, out currentMer);
                        kMerSet.Add(currentMer);
                    }
                    // no degenerate bases (just ACGT) so move to fast path 
                    if (expandedMers.Count == 1)
                    {
                        previousMer = currentMer;
                        kMerWasNormal = true;
                    }
                }
            }

            return kMerSet.Count;
        }
    }
}
