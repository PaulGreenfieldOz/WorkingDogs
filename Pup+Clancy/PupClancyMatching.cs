using System.Collections.Frozen;
using System.Diagnostics;
using System.Reflection.PortableExecutable;
using System.Threading.Tasks;
using WorkingDogsCore;

namespace PupClancy
{
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
        costTwoSubs = 2,
        costGapOpen = 2,
        costGapExtend = 3
    }

    class SharedCode
    {
        // PERF constants
        const int reducedMatch = 0;
        const int reducedX1Sub = 1;
        const int reducedX1InDel = 2;
        const int reducedX2Sub = 1;
        const int reducedX2Indel = 2;
        const int reducedGapIndel = 3;
        const int reducedGapInversion = 0;

        // direction of scanning
        const int plusStrand = 0;
        const int minusStrand = 1;

        // global constants (most settable by parameters)
        public static int largeVariantSetSize = 0;
        public static int smallVariantSetSize = 0;

        public static int maxMismatches = 5;
        public static int maxCost = 10;
        public static int maxDownstream = 20;
        public static int smallGap = 2;
        public static int biggerGap = 4;
        public static int longGap = 6;
        public static int maxConsecutiveFixes = 5;
        public static int trailingCostsLength = 10;
        public static int scanCheckInterval = 10;
        public static bool skipping = true;
        public static bool doubleSubs = false;

        static long seedsCheckedCounts = 0;
        static long matchingSeedSingleCounts = 0;
        static long matchingSeedMultipleCounts = 0;
        static long extendedSeedCounts = 0;

        public static bool useSeeds = true;
        public static int seedLength = 11;
        public static ulong seedMask = 0;
        public static int maxSeedMismatches = 5;

        const int noViableVariants = -1;
        const int maxDegenerate = 5;

        public static StreamWriter trace = null;
        public static bool dumpSeeds = true;

#if PERF
        public static long[][] perfCounts;
        public const int noOfCounts = 25;
        public const string perfHeaders = "Exact\tScanTests\tScanSingleMatches\tScanMultMatches\tExt1Tests\tExt1SingleMatches\tExt1MultMatches\tExt2Tests\tExt2Matches\tExt3Tests\tExt3Matches\tMatched\tCostly\tNoMatch\tSeedTest\tSeedViable\tSeedGrown\tScanMTests\tScanMMatches\tScanMMatchedM\tScanBackTracks\tSkipbacks\tRvsExact\tRvsFuzzy\tSkipbackProbs";
        //                          0     1          2                  3                4          5                  6                7          8            9          10           11       12      13       14        15          16         17          18            19             20              21         22        23
        enum PC //  ptGrowSeed
        {
            pcExact, pcScanTests, pcScanSingleMatches, pcScanMultMatches, pcExt1Tests, pcExt1SingleMatches, pcExt1MultMatches, pcExt2Tests, pcExt2Matches, pcExt3Tests, pcExt3Matches, pcMatched, pcCostly, pcNoMatch, pcSeedTest, pcSeedViable, pcSeedGrown, pcScanMTests, pcScanMMatches, pcScanMMatchedM, pcScanBackTracks, pcSkipbacks, pcRvsExact, pcRvsFuzzy, pcSkipbackProbs
        }
        public static Stopwatch[][] perfTimers;
        public const int noOfPerfTimers = 9;
        public const string perfTimersHeaders = "Scan\tTVMScan\tTVM1\tTVM2\tTVM3\tRVSE\tRVSF\tGVMV\tGrowSeed";
        //                                0     1        2     3    4     5     6     7     8
        enum PT
        {
            ptScan, ptTVMScan, ptTVM1, ptTVM2, ptTVM3, ptRVSE, ptRVSF, ptGVMV, ptGrowSeed
        }
#endif

        public static int FindNextDegenerateBase(string sequence, int startIdx)
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

        public static long GetSequences(string genomeFN, int kmerSize, List<string> headers, List<string> sequences, int minLength, bool quiet, bool verbose)
        {
            StreamReader genome = SeqFiles.OpenSeqStream(genomeFN);
            bool EOF = false;
            long totalLength = 0;

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
                                totalLength += contig.Length;
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
                        totalLength += sequence.Length;
                        if (headers != null)
                            headers.Add(header);
                        //trace.WriteLine(sequences.Count + "\t" + header);
                        break;
                    }
                }
            }
            if (verbose)
                Console.WriteLine("loaded " + sequences.Count + " sequences (" + totalLength + "bp) from " + genomeFN);

            return totalLength;
        }

        public static long MatchSourceToTarget(bool exactMatch, VariantsAllowed allowedVariants, int kmerSize, string sourceSeq, TiledGenome target, int sourceIdx, int targetIdx, string sourceName, string targetName,
                                       List<VariantMerSet> kmerSetPool, VariantMerSet outerMerSet, int threadNo, out int timeTaken)
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

            int mersMatchedForSeq = MatchSequenceToTarget(exactMatch, allowedVariants, ScanType.forwardScan, kmerSize, sequence, 0, target, Strand.either, sourceIdx, targetIdx, 0, kmerSetPool, outerMerSet, threadNo);
            matchingMers += mersMatchedForSeq;

            sw.Stop();
            timeTaken = (int)sw.ElapsedMilliseconds;
            return matchingMers;
        }

        public static int MatchSequenceToTarget(bool exactMatchesOnly, VariantsAllowed allowedVariants, ScanType scanDirection, int kmerSize, Sequence sequence, int start, TiledGenome target, Strand strand, int sourceIdx, int targetIdx, int sequenceIdx, List<VariantMerSet> kmerSetPool, VariantMerSet variantMerSet, int threadNo)
        {
            bool self = sourceIdx == targetIdx;
            int matches = 0;
            ulong kMer = 0;                 // valid/matched kmer 
            int kMerNo = 0;                 // kmer number for tracing etc
            int lastMatchingMerNo = 0;      // last time we had a match of any form (used to control all/last variants)
            int lastMatchingFBI = 0;        // FBI value for the last match (E or V). This is the index of the first unmatched base of a possible unmatched run. 
            bool inMismatchRegion = false;  // in a mismatch region (1 or more consecutive 'no match' k-mers)
            MatchMode matchingMode = scanDirection == ScanType.forwardScan ? MatchMode.scanning : MatchMode.continuing;
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

            MerHashSet? currentTargetMers = null;   // current set of target kMers (plus, minus or all)
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
                    Console.WriteLine("kmerNo off end of sequence:" + "kmerNo=" + kMerNo + " " + scanDirection +
                                        " (" + sourceIdx + "," + sequenceIdx + ")." + sequence.Length + " to " + targetIdx + ".");
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
                // breakpoint
                //if (kMerNo == 297 && sourceIdx != targetIdx)
                //    Debugger.Break();

                // first try for an exact match
                // ----------------------------
                if (currentTargetMers!.Contains(kMer))
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
                    int droppedCost = (int)Costs.costMatch;
                    previousVariantType = VariantType.none;

                    if (trailingCosts.Count == trailingCostsLength)
                    {
                        droppedCost = trailingCosts[0];
                        trailingCosts.RemoveAt(0);
                    }
                    cumulativeCost = cumulativeCost - (int)droppedCost + (int)Costs.costMatch;
                    trailingCosts.Add((int)Costs.costMatch);

                    if (trailingMers.Count == trailingCostsLength)
                        trailingMers.RemoveAt(0);
                    trailingMers.Add(kMer);
#if ( LLTRACE || HLTRACE )
                    if (!self)
                    {
                        //if (kmer == 0xD07F904816000000)
                        //    Debugger.Break();
                        trace.WriteLine(kMerNo + ": E\t" + kMers.ExpandMer(kMer, kmerSize) + " " + fbi +
                            " +" + (fbi == sequence.Length ? "end" : sequence.Bases[fbi].ToString()) + " " + currentStrand.ToString() + " " +
                            (scanDirection == ScanType.reverseScan ? (" [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                            " cc=" + cumulativeCost + " cm=" + consecutiveExactMatches);
                    }
#endif
                    // set current strand for the matching kMer if we're scanning, rather than continuing
                    if (currentStrand == Strand.either)
                        SetCurrentStrand(kMer, target, out currentStrand, out currentTargetMers);

                    kMerNo++;       // deferred until after trace

                    // found an exact match while scanning so it's possible that inexact matches preceding this kMer could also have succeeded (but were too expensive to try)
                    // so we'll scan/match backwards from here, allowing mismatches and continue only until the first mismatch is encountered
                    if (inMismatchRegion && scanDirection == ScanType.forwardScan && !exactMatchesOnly)
                    {
#if ( LLTRACE || HLTRACE )
                        trace.WriteLine(kMerNo + ": rescan after exact match");
#endif
#if PERF
                        perfCounts[threadNo][(int)PC.pcRvsExact]++;
                        perfTimers[threadNo][(int)PT.ptRVSE].Start();
#endif
                        variantMerSet.Reset();
                        matches += RescanMismatchRegion(exactMatchesOnly, allowedVariants, kmerSize, sequence, target, currentStrand, sourceIdx, targetIdx, sequenceIdx, kMer, lastMatchingFBI, fbi, kmerSetPool, variantMerSet, threadNo);
#if PERF
                        perfTimers[threadNo][(int)PT.ptRVSE].Stop();
#endif
                    }

                    //if (kMer == 0x4597950186400000)
                    //    Debugger.Break();

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
#if (LLTRACE)
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
                            // trailingMers holds the last 'few' exact matching kMers. These must be consecutive as trailingMers is cleared on a mismatch.
                            ulong trailingMer = trailingMers[i];
                            variantMerSet.Reset();
                            variantMerSet.AddViable(VariantType.none, trailingMer, 0, 0, (int)Costs.costMatch);
                            skippedExactMatches++;

                            // generate all viable alternatives to this exact-matching kMer (including the kMer itself)
                            if (CheckForViableVariants(trailingMer, currentTargetMers, kmerSize, variantMerSet))
                            {
                                // found an exact match with alternatives, so generate the metrics on all of them and choose the best
                                int bestAlternativeIdx = TryVariantMers(sequence, fbi - skippedExactMatches, kmerSize, MatchType.full, allowedVariants, currentTargetMers, lengthLeft + skippedExactMatches, 0, 1, kMerNo - skippedExactMatches, maxMismatches,
                                                         smallGap, 0, variantMerSet, kmerSetPool, previousVariantType, previousIndelMerNo, threadNo);

                                // if one of the alternatives is better than the initial match ([0]), switch over to this new kMer and fix up the scan state to match
                                if (bestAlternativeIdx > 0 && variantMerSet.merSet[bestAlternativeIdx].exactMatchesFollowing > variantMerSet.merSet[0].exactMatchesFollowing)
                                {
#if (LLTRACE || HLTRACE)
                                    trace.WriteLine(kMerNo + ": found alternative to exact match. " + kMers.ExpandMer(trailingMer, kmerSize) + "-->" + kMers.ExpandMer(variantMerSet.merSet[bestAlternativeIdx].mer, kmerSize) +
                                        " em=" + variantMerSet.merSet[0].exactMatchesFollowing + "-->" + variantMerSet.merSet[bestAlternativeIdx].exactMatchesFollowing); ;
#endif
                                    for (int t = i; t < trailingMers.Count; t++)
                                    {
                                        trailingMers.RemoveAt(t);
                                        cumulativeCost -= (int)trailingCosts[t];
                                        trailingCosts.RemoveAt(t);
                                    }

                                    fbi = fbi - skippedExactMatches;
                                    matches = matches - skippedExactMatches;
                                    consecutiveExactMatches = consecutiveExactMatches - skippedExactMatches;
                                    kMerNo -= skippedExactMatches;
                                    kMer = variantMerSet.merSet[bestAlternativeIdx].mer;
                                    lengthLeft += skippedExactMatches;
                                    previousVariantType = VariantType.none;
                                    previousIndelMerNo = -1;
                                    indelBalance = 0;
                                    retryCurrentMer = true;
#if (LLTRACE || HLTRACE)
                                    trace.WriteLine(kMerNo + ": skipping back over " + skippedExactMatches + " matched variants at " + kMerNo + " " + kMers.ExpandMer(trailingMer, kmerSize) + " --> " + kMers.ExpandMer(variantMerSet.merSet[bestAlternativeIdx].mer, kmerSize));
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
#if (LLTRACE )
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

#if (DEBUG || LLTRACE || HLTRACE)
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

                        int noOfVariantMers = GenerateViableMerVariants(matchType, allowedVariants, previousVariantType, checkFor2xSubs, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, sequence, fbi, lastMatchingFBI,
                                                                        (kMerNo - previousIndelMerNo), 0, currentTargetMers, target.seedsIndex, target.seedsContexts, smallGap, lastScanbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out lastScanbackFBI);
#if PERF
                        perfTimers[threadNo][(int)PT.ptGVMV].Stop();
#endif
                        // backtrack after a successful (but occasional scanning match) and fall into the checking-variants code
                        if (fbiAdjustment > 0)
                        {
#if LLTRACE || HLTRACE
                            trace.WriteLine(kMerNo + ": skipping back to " + (kMerNo - fbiAdjustment) + " (" + kMers.ExpandMer(overlookedMer, kmerSize) + ")");
#endif
                            consecutiveNoMatches -= fbiAdjustment;
                            fbi -= fbiAdjustment;
                            kMerNo -= fbiAdjustment;
                        }

                        if (noOfVariantMers == 1 && matchingMode == MatchMode.continuing)
                        {
                            foundReplacement = true;
                            foundStrongReplacement = true;
                            bestReplacement = variantMerSet.merSet[0];
                            costOfVariant = (int)Costs.costSub;
#if PERF
                            perfCounts[threadNo][(int)PC.pcScanSingleMatches]++;
#endif
#if LLTRACE
                            trace.WriteLine(kMerNo + ": " + kMers.ExpandMer(kMer, kmerSize) + " (1)--> " + kMers.ExpandMer(bestReplacement.mer, kmerSize));
#endif
                        }

                        //if (matchingMode == MatchMode.scanning && noOfVariantMers == 1)
                        //    Debugger.Break();

                        // if there's more than one viable variant or we're being cautious, consider their followers
                        if (noOfVariantMers > 1 || (matchingMode == MatchMode.scanning && noOfVariantMers == 1) || (allowedVariants == VariantsAllowed.subOnly && noOfVariantMers == 1))
                        {
#if PERF
                            perfTimers[threadNo][matchingMode == MatchMode.scanning ? (int)PT.ptTVMScan : (int)PT.ptTVM1].Start();
#endif
#if LLTRACE
                            string variants = "";
                            for (int i = 0; i < variantMerSet.setSize; i++)
                                variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                            trace.WriteLine(kMerNo + ": TVM1 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif
                            // first try for an fairly close match (small number of mismatches, small gaps)
                            replacementIdx = TryVariantMers(sequence, fbi, kmerSize, (matchingMode == MatchMode.scanning ? MatchType.scanning : MatchType.full),
                                                           allowedVariants, currentTargetMers, lengthLeft, 0, 1, kMerNo, maxMismatches,
                                                           smallGap, 0, variantMerSet, kmerSetPool, previousVariantType, previousIndelMerNo, threadNo);
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
                                int lengthCovered = bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing;
                                foundReplacement = (bestReplacement.exactMatchesFollowing >= bestReplacement.mismatchesFollowing) && (lengthCovered >= lengthLeft / 2);
                                foundStrongReplacement = foundReplacement && lengthCovered >= lengthLeft && bestReplacement.mismatchesFollowing < maxMismatches;

                                costOfVariant = bestReplacement.costOfVariant;
#if DEBUG
                                if (foundStrongReplacement)
                                    whenMatched = "x1";
#endif
                            } // first attempt at choosing between viable variants

                            // if we didn't find a strong replacement - and we're continuing a forward matching region - we'll try a more expensive exploration with more mismatches allowed
                            if (matchingMode == MatchMode.continuing && foundReplacement && !foundStrongReplacement && trailingCosts.Count > trailingCostsLength / 2 && cumulativeCost < trailingCosts.Count / 2)
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcExt2Tests]++;
                                perfTimers[threadNo][(int)PT.ptTVM2].Start();
#endif

                                //if (kMerNo == 1622)
                                //    Debugger.Break();

                                // generate the variants again, with 'full' and a larger allowable gap
                                variantMerSet.Reset();
                                noOfVariantMers = GenerateViableMerVariants(MatchType.full, allowedVariants, previousVariantType, false, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, sequence, fbi, lastMatchingFBI,
                                                                            (kMerNo - previousIndelMerNo), 0, currentTargetMers, target.seedsIndex, target.seedsContexts, biggerGap, skipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out skipbackFBI);
#if LLTRACE
                                variants = "";
                                for (int i = 0; i < variantMerSet.setSize; i++)
                                    variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                                trace.WriteLine(kMerNo + ": TVM2 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif

                                replacementIdx = TryVariantMers(sequence, fbi, kmerSize, MatchType.full, allowedVariants,
                                                                currentTargetMers, lengthLeft, 0, 1, kMerNo, maxMismatches * 2,
                                                                biggerGap, 0, variantMerSet, kmerSetPool, previousVariantType, previousIndelMerNo, threadNo);
#if PERF
                                perfTimers[threadNo][(int)PT.ptTVM2].Stop();
                                if (replacementIdx > 0)
                                    perfCounts[threadNo][(int)PC.pcExt2Matches]++;
#endif  
                                if (replacementIdx >= 0)
                                {
                                    bestReplacement = variantMerSet.merSet[replacementIdx];
                                    int lengthCovered = bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing;
                                    foundReplacement = (bestReplacement.exactMatchesFollowing > bestReplacement.mismatchesFollowing) && (lengthCovered > lengthLeft / 2);
                                    foundStrongReplacement = foundReplacement && lengthCovered >= lengthLeft && bestReplacement.mismatchesFollowing <= maxMismatches * 2;
#if DEBUG
                                    if (foundStrongReplacement)
                                        whenMatched = "x2";
#endif
                                }
                            } // second attempt

                            // if this is the first mismatch after a reasonable looking run of matches, try for a long gap or even reversal
                            if (matchingMode == MatchMode.continuing && scanDirection == ScanType.forwardScan && !foundReplacement && trailingCosts.Count >= (trailingCostsLength * 2 / 3) && cumulativeCost < trailingCostsLength / 2)
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcExt3Tests]++;
                                perfTimers[threadNo][(int)PT.ptTVM3].Start();
#endif
                                variantMerSet.Reset();
                                noOfVariantMers = GenerateViableMerVariants(MatchType.spanGap, allowedVariants, previousVariantType, false, kMer, kmerSize, retryCurrentMer, variantMerSet, kmerSetPool, sequence, fbi, lastMatchingFBI,
                                                                            (kMerNo - previousIndelMerNo), 0, currentTargetMers, target.seedsIndex, target.seedsContexts, longGap, skipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out skipbackFBI);
#if LLTRACE
                                variants = "";
                                for (int i = 0; i < variantMerSet.setSize; i++)
                                    variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + ") ";
                                trace.WriteLine(kMerNo + ": TVM3 " + kMers.ExpandMer(kMer, kmerSize) + " --> " + variants);
#endif
                                replacementIdx = TryVariantMers(sequence, fbi, kmerSize, MatchType.spanGap, allowedVariants,
                                                                target.mersFromGenomeEither, lengthLeft, 0, 1, kMerNo, maxMismatches,
                                                                longGap, 0, variantMerSet, kmerSetPool, previousVariantType, previousIndelMerNo, threadNo);
#if PERF
                                perfTimers[threadNo][(int)PT.ptTVM3].Stop();
                                if (replacementIdx > 0)
                                    perfCounts[threadNo][(int)PC.pcExt3Matches]++;
#endif
                                if (replacementIdx >= 0)
                                {
                                    bestReplacement = variantMerSet.merSet[replacementIdx];
                                    int lengthCovered = bestReplacement.exactMatchesFollowing + bestReplacement.mismatchesFollowing;
                                    foundReplacement = (bestReplacement.exactMatchesFollowing > bestReplacement.mismatchesFollowing) && (lengthCovered > lengthLeft / 2);
                                    foundStrongReplacement = foundReplacement && lengthCovered >= lengthLeft && bestReplacement.mismatchesFollowing <= maxMismatches;

                                    bool reversal = false;
                                    if (foundStrongReplacement)
                                    {
                                        bestReplacement = variantMerSet.merSet[replacementIdx];
                                        MerHashSet priorTargetMers = currentTargetMers;
                                        SetCurrentStrand(bestReplacement.mer, target, out currentStrand, out currentTargetMers);
                                        reversal = priorTargetMers != currentTargetMers;
                                    }

#if DEBUG
                                    if (foundStrongReplacement)
                                    {
                                        if (reversal)
                                            whenMatched = "reversed";
                                        else
                                            whenMatched = "gap";
                                    }
#endif
                                }
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
                                costOfVariant = (int)Costs.costTwoSubs;
                            int droppedCost = 0;
                            if (trailingCosts.Count == trailingCostsLength)
                            {
                                droppedCost = (int)trailingCosts[0];
                                trailingCosts.RemoveAt(0);
                            }
                            cumulativeCost = cumulativeCost - (int)droppedCost + (int)costOfVariant;
                            trailingCosts.Add(costOfVariant);

                            if (cumulativeCost <= maxCost || indelBalance > maxDownstream / 2)
                            {
                                matchingMode = MatchMode.continuing;
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
#if (LLTRACE || HLTRACE)
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
#if (LLTRACE || HLTRACE)
                                        trace.WriteLine(kMerNo + ": rescan after fuzzy match");
                                        //if (kMerNo == 759)
                                        //    Debugger.Break();
#endif
#if PERF
                                        perfCounts[threadNo][(int)PC.pcRvsFuzzy]++;
                                        perfTimers[threadNo][(int)PT.ptRVSF].Start();
#endif
                                        matches += RescanMismatchRegion(exactMatchesOnly, allowedVariants, kmerSize, sequence, target, currentStrand, sourceIdx, targetIdx, sequenceIdx, kMer, lastMatchingFBI, fbi, kmerSetPool, variantMerSet, threadNo);
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
#if (LLTRACE || HLTRACE)
                                if (!self)
                                    trace.WriteLine(kMerNo + ": V\t" + kMers.ExpandMer(kMer, kmerSize) +
                                        " +" + (fbi == sequence.Length ? "end" : sequence.Bases[fbi].ToString()) + " " + currentStrand.ToString() + " " +
                                    traceMer.variantType + " (" + traceMer.gapLength + ")" + " bem=" + traceMer.exactMatchesFollowing + " bmm=" + traceMer.mismatchesFollowing + " bidb=" + traceMer.indelBalance +
                                    " " + matchingMode.ToString() + " " + whenMatched + (scanDirection == ScanType.reverseScan ? ("  [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                                    " cost=" + traceMer.costOfVariant + "/" + traceMer.costOfFollowing + " cc=" + cumulativeCost + " cnm=" + consecutiveNoMatches + " idb=" + indelBalance);
#endif
                            }
                            else
                            // cumulative cost getting too high, so revert to source sequence again
                            {
#if PERF
                                perfCounts[threadNo][(int)PC.pcCostly]++;
#endif
#if (LLTRACE || HLTRACE)
                                trace.WriteLine(kMerNo + ": A\t" + kMers.ExpandMer(kMer, kmerSize) + " " + currentStrand.ToString() + " " +
                                    bestReplacement.variantType + " (" + bestReplacement.gapLength + ")" + " bem=" + bestReplacement.exactMatchesFollowing + " bmm=" + bestReplacement.mismatchesFollowing +
                                    " " + matchingMode.ToString() + " " + whenMatched + (scanDirection == ScanType.reverseScan ? ("  [" + kMers.ExpandMer(kMers.ReverseComplement(kMer, kmerSize), kmerSize) + "]") : "") +
                                    " cost=" + bestReplacement.costOfVariant + "/" + bestReplacement.costOfFollowing + " cc=" + cumulativeCost + " cnm=" + consecutiveNoMatches);
#endif
                                foundReplacement = false;
                            } // cumulative cost too high - abandon replacement
                        } // found a suitable replacement
                    } // need to check for variant

                    // either couldn't get either seed match or an acceptable variant
                    if (!foundReplacement)
                    {
                        //if (kMer == 0x81236186D8220000)
                        //    Debugger.Break();
#if (LLTRACE || HLTRACE)
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
#if (LLTRACE || HLTRACE)
                            trace.WriteLine(kMerNo + ": stopping reverse scan");
#endif
                            break;
                        }

                        // restore the kmer from the sequence if necessary
                        if (consecutiveFuzzyMatches > 0)
                        {
                            ulong kMerPrior = kMer;
                            fbi = fbi - consecutiveFuzzyMatches;
                            kMerNo = kMerNo - consecutiveFuzzyMatches;
                            Sequence.CondenseMer(sequence, fbi - kmerSize, kmerSize, out kMer);
#if (LLTRACE || HLTRACE)
                            trace.WriteLine(kMerNo + ": no match reset " + kMers.ExpandMer(kMerPrior, kmerSize) + " --> " + kMers.ExpandMer(kMer, kmerSize));
#endif
                        }

                        consecutiveNoMatches++;
                        consecutiveFuzzyMatches = 0;
                        if (checkForVariant)
                        {
                            inMismatchRegion = true;
                            matchingMode = MatchMode.scanning;
                            ulong kMerPrior = kMer;
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

        private static string ShowVariants(VariantMerSet variantMerSet, int kmerSize)
        {
            string variants = "";
            for (int i = 0; i < variantMerSet.setSize; i++)
                variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + "(" + variantMerSet.merSet[i].variantType + "@" + variantMerSet.merSet[i].costOfVariant + "+" + variantMerSet.merSet[i].exactMatchesFollowing + "+" + variantMerSet.merSet[i].mismatchesFollowing + ") ";
            return variants;
        }

        public static long PrepTargetMerSet(TiledGenome[] targets, int targetNo, string targetFN, string targetName, int kmerSize, List<string> genomeSequences, int minLength, bool verbose, bool quiet)
        {
            long kmersInTarget = 0;

            // first matcher to want this target genome? then load it up its hashtables
            if (targets[targetNo].mersFromGenomePlus == null)
            {
                lock (targets)
                {
                    // some other thread may have loaded this table while we were waiting
                    if (targets[targetNo].mersFromGenomePlus != null)
                    {
                        return 0;
                    }

                    //if (targetNo == 10)
                    //    Debugger.Break();

                    Stopwatch prepSW = new Stopwatch();
                    prepSW.Start();

                    if (genomeSequences.Count == 0)
                    {
                        List<string> headers = new List<string>(10000);
                        GetSequences(targetFN, kmerSize, headers, genomeSequences, minLength, quiet, verbose);
                    }

                    //Console.WriteLine("loading targets[" + targetNo + "]");
                    int totalLength = 0;
                    foreach (string s in genomeSequences)
                        totalLength += s.Length;

                    MerHashSet plusMers = new MerHashSet(totalLength, kmerSize, false);
                    MerHashSet minusMers = new MerHashSet(totalLength, kmerSize, false); 
                    MerHashSet allMers = new MerHashSet(totalLength*2, kmerSize, false); 
                    Dictionary<uint, ulong>? seedsIndex = null;
                    List<ulong>? seedsContexts = null;
                    if (seedLength != 0)
                    {
                        int targetSeedsLength = Math.Min((int)plusMers.Count, 4 ^ seedLength);
                        seedsIndex = new Dictionary<uint, ulong>(targetSeedsLength);
                        seedsContexts = new List<ulong>(totalLength);
                    }

                    // generate the distinct k-mers for the next target genome (not canonical) and build the 'plus' kMer set
                    Stopwatch tilingSW = new Stopwatch();
                    tilingSW.Start();
                    kmersInTarget = GetMersFromSequences(genomeSequences, kmerSize, totalLength, plusMers, minusMers, allMers, seedsContexts);
                    tilingSW.Stop();

                    if (verbose)
                        Console.WriteLine("tiled target (" + targetName + ") in " + tilingSW.Elapsed.TotalSeconds.ToString("F1") + "s" + " " + plusMers.Count + "+" + minusMers.Count + "+" + allMers.Count + " kMers");

                    if (seedsIndex != null)
                    {
                        Stopwatch seedsSW = new Stopwatch();
                        seedsSW.Start();
                        int countSeedContexts = 0;
                        for (int i = 0; i < seedsContexts.Count; i++)
                        {
                            ulong kMer = seedsContexts[i];
                            uint seed = (uint)((kMer & seedMask) >> 32);
                            uint startForSeed = (uint)i;
                            uint countForSeed = 1;
                            while (i < seedsContexts.Count - 1 && seed == (uint)((seedsContexts[i + 1] & seedMask) >> 32))
                            {
                                countForSeed++;
                                i++;
                            }
                            ulong idxForSeed = ((ulong)startForSeed << 32) | (ulong)countForSeed;
                            seedsIndex.Add(seed, idxForSeed);
                            countSeedContexts += (int)countForSeed;
                        }
                        seedsSW.Stop();
                        if (verbose)
                            Console.WriteLine(seedsIndex.Count + " seeds with " + countSeedContexts + " contexts (" + seedsContexts.Count + ") in" + seedsSW.Elapsed.TotalSeconds.ToString("F1") + "s");
                    }

                    string seedMsg = useSeeds ? (seedsIndex.Count + " " + seedLength + "-mer seeds") : "no seeds";
                    if (useSeeds)
                    {
                        targets[targetNo].seedsIndex = FrozenDictionary.ToFrozenDictionary<uint, ulong>(seedsIndex);
                        seedsIndex = null;
                        targets[targetNo].seedsContexts = seedsContexts;
                    }
                    targets[targetNo].mersFromGenomeMinus = minusMers;
                    targets[targetNo].mersFromGenomeEither = allMers;
                    targets[targetNo].mersFromGenomePlus = plusMers;        // last act before lock is released
                    prepSW.Stop();

                    if (!quiet)
                        Console.WriteLine("loaded target #" + targetNo + " (" + targetName + ")" + " plus=" + plusMers.Count + " both=" + allMers.Count + " (" + seedMsg + ")" + " in " + prepSW.Elapsed.TotalSeconds.ToString("F1") + "s");
                }
            }

            return kmersInTarget;
        }

        public static void ReleaseTargetMers(int targetNo, string targetName, TiledGenome[] targets)
        {
            lock (targets)
            {
                targets[targetNo].refCount--;
                if (targets[targetNo].refCount == 0)
                {
                    targets[targetNo].mersFromGenomePlus = null;
                    targets[targetNo].mersFromGenomeMinus = null;
                    targets[targetNo].mersFromGenomeEither = null;
                    targets[targetNo].seedsIndex = null;
                    targets[targetNo].seedsContexts = null;
                    GC.Collect(2);
                    Console.WriteLine("Finished matching to " + targetName);
                }
            }
        }

        public static long GetMersFromSequences(List<string> sequences, int kmerSize, int totalLength, MerHashSet plusMers, MerHashSet minusMers, MerHashSet allMers, List<ulong> seedsContexts)
        {
            List<ulong> allPlusMers = new List<ulong>((int)totalLength);
            List<ulong> allMinusMers = new List<ulong>((int)totalLength);

            char[] degenerateBases;
            Dictionary<char, List<char>> degenerateBaseExpansions;
            kMers.InitialiseDegenerateBaseTables(out degenerateBases, out degenerateBaseExpansions);

            //for (int c = 0; c < sequences.Count; c++)
            Parallel.For(0, sequences.Count, c =>
            {
                List<ulong> kMersFromSeqPlus = new List<ulong>(sequences[c].Length);
                int mersFound = kMers.GenerateExpandedMersFromRead(sequences[c], kmerSize, kMersFromSeqPlus, degenerateBases, degenerateBaseExpansions);
                lock (allPlusMers)
                {
                    allPlusMers.AddRange(kMersFromSeqPlus);
                    for (int i = 0; i < kMersFromSeqPlus.Count; i++)
                        allMinusMers.Add(kMers.ReverseComplement(kMersFromSeqPlus[i], kmerSize));
                }
            });
            //}

            Parallel.For(0, 4, s =>
            {
                if (s == 0)
                {
                    for (int m = 0; m < allPlusMers.Count; m++)
                        plusMers.AddIfNotPresent(allPlusMers[m]);
                }

                if (s == 1)
                {
                    for (int m = 0; m < allMinusMers.Count; m++)
                        minusMers.AddIfNotPresent(allMinusMers[m]);
                }

                if (s == 2)
                {
                    for (int m = 0; m < allPlusMers.Count; m++)
                        allMers.AddIfNotPresent(allPlusMers[m]);
                    for (int m = 0; m < allMinusMers.Count; m++)
                        allMers.AddIfNotPresent(allMinusMers[m]);
                }

                if (s == 3 && seedsContexts != null)
                {
                    seedsContexts.AddRange(allPlusMers);
                    seedsContexts.AddRange(allMinusMers);
                    seedsContexts.Sort();  
                    DeduplicateSortedList(seedsContexts);
                }
            });

            return (long)allPlusMers.Count;
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

        static int DeduplicateSortedList(List<ulong> kMersList)
        {
            int duplicates = 0;
            int nextDistinctIdx = 1;
            for (int i = 1; i < kMersList.Count; i++)
            {
                if (kMersList[i - 1] == kMersList[i])
                    duplicates++;
                else
                {
                    kMersList[nextDistinctIdx] = kMersList[i];
                    nextDistinctIdx++;
                }
            }

            kMersList.RemoveRange(nextDistinctIdx, duplicates);
            return duplicates;
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
#if LLTRACE
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

        private static int RescanMismatchRegion(bool exactMatching, VariantsAllowed allowedVariants, int kmerSize, Sequence sequence, TiledGenome target, Strand strand, int sourceIdx, int targetIdx, int sequenceIdx, ulong kmer, int lastMatchingFBI, int fbi,
                                                List<VariantMerSet> kmerSetPool, VariantMerSet outerMerSet, int threadNo)
        {
            int matches = 0;

            int rescanStart = lastMatchingFBI - kmerSize + 1;
            if (rescanStart < 0)
                rescanStart = 0;

            int rescanLength = fbi - lastMatchingFBI - 1;
            if (rescanLength < kmerSize + 1)
                return 0;

            Sequence mismatchRegion = sequence.SubSeq(rescanStart, rescanLength);
#if (LLTRACE)
            trace.WriteLine(mismatchRegion.ToString());
#endif
            mismatchRegion.Replace(mismatchRegion.Length - kmerSize + 1, kmer, kmerSize - 1);
            mismatchRegion.ReverseComplement();

#if (LLTRACE || HLTRACE)
            trace.WriteLine("reverse scan of mismatch region." + " lfbi=" + lastMatchingFBI + " fbi=" + fbi + " length=" + mismatchRegion.Length);
            trace.WriteLine(mismatchRegion.ToString());
#endif
            Strand currentStrand;
            if (strand == Strand.plus)
                currentStrand = Strand.minus;
            else
                currentStrand = Strand.plus;

            matches = MatchSequenceToTarget(exactMatching, allowedVariants, ScanType.reverseScan, kmerSize, mismatchRegion, 0, target, currentStrand, sourceIdx, targetIdx, sequenceIdx, kmerSetPool, outerMerSet, threadNo);

#if (LLTRACE || HLTRACE)
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

        private static int TryVariantMers(Sequence sequence, int fbi, int kmerSize, MatchType matchType, VariantsAllowed allowedVariants, MerHashSet targetMers,
                                          int lengthLeft, int consecutiveFixes, int startingIndent, int kmerNo, int allowedMismatches, int allowedGapLength, int basesAdded, VariantMerSet variantMerSet,
                                          List<VariantMerSet> kmerSetPool, VariantType previousVariantType, int previousIndelMerNo, int threadNo)
        {
#if LLTRACE
            string variants = "";
            for (int i = 0; i < variantMerSet.setSize; i++)
                variants += kMers.ExpandMer(variantMerSet.merSet[i].mer, kmerSize) + " " + variantMerSet.merSet[i].variantType + " ";
            trace.WriteLine(Indent(startingIndent) + "tvm @ " + kmerNo + " ll=" + lengthLeft + " fbi=" + fbi + " amm=" + allowedMismatches +
                                    " +" + ((fbi == sequence.Length) ? "end" : sequence.Bases[fbi].ToString()) +
                                    " " + previousVariantType.ToString() + " pid=" + previousIndelMerNo + " cf=" + consecutiveFixes + " lg=" + allowedGapLength + " (" + variants + ")");
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
                int mismatches = 1;

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
#if LLTRACE
                    trace.WriteLine(Indent(startingIndent + 1) + "close @ " + kmerNo + " " + kMers.ExpandMer(currentVariant.mer, kmerSize) + " " + previousVariantType +
                                    " @ " + previousIndelMerNo + " (" + currentVariant.variantType + " @" + kmerNo + ")");
#endif
                    currentVariant.variantType = VariantType.invalid;
                    continue;
                }

                // don't allow indels to get out of control. Too many of the same type of indels in a row is probably a sign of going astray.
                if (currentVariant.variantType == VariantType.ins)
                {
                    currentVariant.indelBalance = -currentVariant.gapLength;
                    basesAdded -= currentVariant.gapLength;
                }
                if (currentVariant.variantType == VariantType.del)
                {
                    currentVariant.indelBalance = +currentVariant.gapLength;
                    basesAdded += currentVariant.gapLength;
                }

                // we've got a viable k-mer variant - one that matches a k-mer in the target set, so find out how well it matches downstream
                CountFollowers(currentVariant, sequence, adjustedFBI, kmerSize, targetMers, allowedVariants, currentVariant.variantType, lengthLeft, consecutiveFixes + 1, currentIndelMerNo, indent + 1, kmerNo,
                               allowedMismatches, allowedGapLength, basesAdded, kmerSetPool, threadNo);

                exactMatches = currentVariant.exactMatchesFollowing;
                mismatches = currentVariant.mismatchesFollowing;
                int followers = exactMatches + mismatches;
                int costOfFollowing = currentVariant.costOfFollowing + currentVariant.costOfVariant;

                if (followers > 0)
                {
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

                        // indel-rich repairs can overstate exact matches, so boost sub slightly to force a cost comparison
                        if ((currentVariant.variantType == VariantType.ins || currentVariant.variantType == VariantType.del) && currentVariant.followingIndel != VariantType.none)
                            exactMatches--;

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

                        // if we got to the end of the match region with this variant, any other variants will have to do at least as well to be viable
                        if (allowedMismatches > mismatches)
                            allowedMismatches = mismatches + 1;
                    }
                }
            }

            // if no viable variants, continue the main scanning loop
            if (viableMerCount == 0)
            {
#if LLTRACE
                trace.WriteLine(Indent(startingIndent) + "tvr @ " + kmerNo + " no viable variants ");
#endif
                return noViableVariants;
            }
#if LLTRACE
            trace.WriteLine(Indent(startingIndent + 1) + "tvc @ " + kmerNo + " " + "ll=" + lengthLeft + " lc=" + OKLowestCost + " bem=" + OKHighesttExactMatches + " best=" + bestIdx);
            for (int v = 0; v < noOfVariantMers; v++)
                trace.WriteLine(Indent(startingIndent + 1) + "tvc @ " + kmerNo + " " +
                            kMers.ExpandMer(variantMers[v].mer, kmerSize) + " @ " + (fbi + variantMers[v].followingBaseIncrement) + " "
                            + variantMers[v].variantType + "(" + variantMers[v].gapLength + ")" +
                            " bem=" + variantMers[v].exactMatchesFollowing + " bmm=" + variantMers[v].mismatchesFollowing + " fid=" + variantMers[v].followingIndel +
                            " cost=" + variantMers[v].costOfVariant + "/" + variantMers[v].costOfFollowing + " cf=" + consecutiveFixes + " idb=" + variantMers[v].indelBalance);
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
                    int costOfThisVariant = variantMers[v].costOfFollowing + variantMers[v].costOfVariant;
                    if (followers > highestFollowers)
                    {
                        highestFollowers = followers;
                        lowestCost = costOfThisVariant;
                    }
                    if (followers == highestFollowers && variantMers[v].costOfFollowing < lowestCost)
                        lowestCost = costOfThisVariant;
                }

                // now find the first variant that gives these values
                for (int v = 0; v < noOfVariantMers; v++)
                {
                    // ignore any non-viable variants
                    if (variantMers[v].variantType == VariantType.invalid)
                        continue;

                    int followers = variantMers[v].exactMatchesFollowing + variantMers[v].mismatchesFollowing;
                    if (followers == highestFollowers && (variantMers[v].costOfFollowing + variantMers[v].costOfVariant) == lowestCost)
                    {
                        bestIdx = v;
                        break;
                    }
                }
            }
#if LLTRACE
            if (bestIdx >= 0)
                trace.WriteLine(Indent(startingIndent) + "tvr @ " + kmerNo + " " +
                                kMers.ExpandMer(variantMers[bestIdx].mer, kmerSize) + " @ " + fbi + " "
                                + variantMers[bestIdx].variantType + "(" + variantMers[bestIdx].gapLength + ")" +
                                " bem=" + variantMers[bestIdx].exactMatchesFollowing + " bmm=" + variantMers[bestIdx].mismatchesFollowing +
                                " cost=" + variantMers[bestIdx].costOfVariant + "/" + variantMers[bestIdx].costOfFollowing + " cfix=" + consecutiveFixes);
            else
                trace.WriteLine(Indent(startingIndent) + "tvr @ " + kmerNo + "no best choice");
#endif

            return bestIdx;
        }

        // count how many exact matches there are downstream from a given k-mer in a given sequence context. 
        // The scan is length-limited to avoid going to the end of the sequence every time, and mismatch-limited so we can give up early
        private static void CountFollowers(Mer kMer, Sequence sequence, int nbi, int kmerSize, MerHashSet targetMers, VariantsAllowed allowedVariants, VariantType previousVariantType,
                                           int lengthLeft, int consecutiveFixes, int previousIndelMerNo, int indent, int kMerNo, int allowedMismatches, int allowedGap, int basesAdded, List<VariantMerSet> kmerSetPool, int threadNo)
        {
#if LLTRACE
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
            int costOfFollowing = (int)kMer.costOfVariant;
            if (kMer.gapLength > smallGap)
                costOfFollowing++;
            if (kMer.gapLength > biggerGap)
                costOfFollowing++;
            int basesToFirstChange = 0;
            int usedMismatches = 0;
            VariantType followingIndel = VariantType.none;

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
                    previousVariantType = VariantType.none;
#if LLTRACE
                    trace.WriteLine(Indent(indent + 1) + "cf= @ " + kMerNo + " nbi=" + nbi + " em=" + exactMatches + " ll=" + lengthLeft + " " + kMers.ExpandMer(nextMer, kmerSize) + " (" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ")");
#endif
                }
                else
                {
                    allowedMismatches--;
                    usedMismatches++;

                    if (allowedMismatches > 0 && lengthLeft > 0 && consecutiveFixes < maxConsecutiveFixes)
                    {
                        lengthLeft--;
                        nbi++;
                        kMerNo++;
#if LLTRACE
                        //trace.WriteLine(Indent(indent + 1) + "tvm? @ " + kmerNo + " nbi=" + nbi + " ll=" + lengthLeft + " amm=" + allowedMismatches + " tf=" + targetExactMatches + " f=" + exactMatches + " ?=" + (lengthLeft >= targetExactMatches) + " " + kMers.ExpandMer(nextMer, kmerSize) + " (" + ((nbi == sequence.Length) ? "end" : sequence.Bases[nbi].ToString()) + ")");
#endif
                        VariantMerSet variantMerSet = GetVariantMerSet(kmerSetPool);
                        int fbiAdjustment;
                        int lastSkipbackFBI = 0;
                        ulong overlookedMer;
                        // this is an 'extension' match, not a 'scan'
                        int noOfMerVariants = GenerateViableMerVariants(MatchType.full, allowedVariants, previousVariantType, false, nextMer, kmerSize, false, variantMerSet, null, sequence, nbi, nbi,
                                                                        (kMerNo - previousIndelMerNo), 0, targetMers, null, null, allowedGap, lastSkipbackFBI, threadNo, out fbiAdjustment, out overlookedMer, out lastSkipbackFBI);
                        int bestIdx = TryVariantMers(sequence, nbi, kmerSize, MatchType.full, allowedVariants, targetMers, lengthLeft, consecutiveFixes, indent + 1, kMerNo, allowedMismatches,
                                                     allowedGap, basesAdded, variantMerSet, kmerSetPool, kMer.variantType, previousIndelMerNo, threadNo);
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

            // copy follower results back to variant kMer
            kMer.exactMatchesFollowing = exactMatches;
            kMer.mismatchesFollowing = usedMismatches;
            kMer.costOfFollowing = costOfFollowing;
            kMer.followingIndel = followingIndel;
#if LLTRACE
            trace.WriteLine(Indent(indent) + "cfr @ " + startingKMerNo + " " + kMer.variantType + "(" + kMer.gapLength + ")"
                + " em=" + exactMatches + " umm=" + usedMismatches + " sll=" + startingLengthLeft +
                " ll=" + lengthLeft + (lengthLeft == 0 ? " end" : "") + " cost=" + kMer.costOfVariant + "/" + kMer.costOfFollowing + " cf=" + consecutiveFixes + " fid=" + followingIndel + " cm=" + consecutiveMatches);
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

        private static void ReturnVariantMerSet(List<VariantMerSet> kmerSetPool, VariantMerSet kmerSet)
        {
            kmerSet.Reset();
            //Console.WriteLine("return VMS# " + variantMerSet.merSetNo);
            kmerSetPool.Add(kmerSet); 
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

    private static int GenerateViableMerVariants(MatchType matchType, VariantsAllowed allowedVariants, VariantType previousVariantType, bool doubleSub, ulong kMer, int kmerSize, bool tryCurrentMer, VariantMerSet variantMerSet, List<VariantMerSet> kmerSetPool, 
                                                 Sequence sequence, int fbi, int lastMatchingFBI, int distanceToLastIndel, int basesAdded, MerHashSet targetMers, FrozenDictionary<uint, ulong>? seedsIndex, List<ulong> seedsContexts, int allowedGapLength, int previousScanbackFBI, int threadNo,
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
                // don't allow consecutive indels, or a sub immediately followed by an indel (the sub has just deferred the indel)
                if (allowedVariants == VariantsAllowed.subInsDel && distanceToLastIndel > 1 && previousVariantType != VariantType.sub)
                {
                    delVariantsAdded = GenerateViableMerDelVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength, basesAdded);
                    insVariantsAdded = GenerateViableMerInsVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength, basesAdded);
                }
                return variantMerSet.setSize;

            case MatchType.scanning:
                    // if the genomes/genes are small enough, use seed matching/extension during scanning. Seeds have to be both fairly small and sparse
                    if (useSeeds)
                    {
#if PERF
                        perfCounts[threadNo][(int)PC.pcSeedTest]++;
#endif
                        uint seed = (uint)((kMer & seedMask) >> 32);
                        bool seedViable = false;
                        ulong seedContext; 
                        if (seedsIndex.TryGetValue(seed, out seedContext))
                            seedViable = true;

                        if (seedViable)
                        {
#if PERF
                            perfCounts[threadNo][(int)PC.pcSeedViable]++;
                            perfTimers[threadNo][(int)PT.ptGrowSeed].Start();
#endif
                            int idx = (int)(seedContext >> 32);
                            int count = (int)(seedContext & 0xFFFFFFFF);
                            subVariantsAdded = GenerateViableSubVariantsFromSeeds(kMer, kmerSize, sequence.Bases[fbi], targetMers, seedsContexts, idx, count, variantMerSet);
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
#if LLTRACE
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
                if (allowedVariants == VariantsAllowed.subInsDel && distanceToLastIndel > 1 && previousVariantType != VariantType.sub)
                {
                    delVariantsAdded = GenerateViableMerDelVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength, basesAdded);
                    insVariantsAdded = GenerateViableMerInsVariants(kMer, kmerSize, sequence, fbi, variantMerSet, targetMers, allowedGapLength, basesAdded);
                }
                return variantMerSet.setSize;
        }
        return 0;
    }

    private static int GenerateViableSubVariantsFromSeeds(ulong kmer, int kmerSize, char followingBaseInSequence, MerHashSet targetMers, List<ulong> seedsContexts, int seedFirstIdx, int seedCount, VariantMerSet variantMerSet)
    {
        int viableVariants = 0;
        ulong followingBaseShifted = (ulong)kMers.BaseCharToInt(followingBaseInSequence) << (32 - kmerSize) * 2;

        for (int i = 0; i < seedCount; i++)
        {
            ulong context = seedsContexts[seedFirstIdx + i];
            ulong followingMer = (context << 2) | followingBaseShifted;
            int targetMismatches = CountMismatches(kmer, context, kmerSize, maxSeedMismatches);
            if (targetMismatches <= maxSeedMismatches && targetMers.Contains(followingMer))
            {
                viableVariants++;
                variantMerSet.AddViable(VariantType.sub, context, 0, 0, targetMismatches);
            }
        }

        //if (viableVariants > 0)
        //    Interlocked.Increment(ref extendedSeedCounts);
        return viableVariants;
    }

    private static int CountMismatches(ulong kmer, ulong variant, int kmerSize, int maxMismatches)
    {
        int mismatches = 0;
        ulong kmerShifted = kmer;
        ulong variantShifted = variant;
        for (int i = 0; i < kmerSize; i++)
        {
            if ((kmerShifted & 0xc000000000000000) != (variantShifted & 0xc000000000000000))
                mismatches++;
            if (mismatches > maxMismatches)
                break;
            kmerShifted = kmerShifted << 2;
            variantShifted = variantShifted << 2;
        }

        return mismatches;
    }

    private static string DumpVMS(VariantMerSet vms, int kmerSize)
    {
        string result = "";
        for (int i = 0; i < vms.setSize; i++)
            result += i + "\t" + kMers.ExpandMer(vms.merSet[i].mer, kmerSize) + "(" + vms.merSet[i].variantType + ",em=" + vms.merSet[i].exactMatchesFollowing + ",mm=" + vms.merSet[i].mismatchesFollowing + ")\n";
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
    private static int GenerateViableMerDelVariants(ulong mer, int kmerSize, Sequence sequence, int fbi, VariantMerSet variantMerSet, MerHashSet targetMers, int allowedGapLength, int basesAdded)
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

        // add viable single-insertion variants (unless we've previously added a base)
        if (basesAdded != 1)
            AddViableDelVariants(mer & upperMerMask, 1, kmerSize, targetMers, variantMerSet);

        // add extended variants of the del variants already added
        for (int g = 1; g < allowedGapLength; g++)
        {
            // remember last one added
            int lastDelVariantAdded = variantMerSet.setSize;

            if (basesAdded != g)
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
    private static int GenerateViableMerInsVariants(ulong mer, int kmerSize, Sequence sequence, int fbi, VariantMerSet variantMerSet, MerHashSet targetMers, int allowedGapLength, int basesAdded)
    {
        int variantsAdded = 0;

        int consecutiveIns = 0;
        for (int g = 0; g < allowedGapLength; g++)
        {
            // don't delete bases if we've just added the same number
            if (basesAdded == -g)
                continue;

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

    public static void GetPathFN(string readsFN, out string readsPath, out string readsFNP)
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
            this.indelBalance = mer.indelBalance;
            this.costOfVariant = mer.costOfVariant;
            this.costOfFollowing = mer.costOfFollowing;
            this.followingIndel = mer.followingIndel;
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

            merSet[setSize].mer = mer;
            merSet[setSize].variantType = variantType;
            merSet[setSize].exactMatchesFollowing = 0;
            merSet[setSize].mismatchesFollowing = 0;
            merSet[setSize].followingBaseIncrement = followingBaseIncrement;
            merSet[setSize].gapLength = gapLength;
            merSet[setSize].indelBalance = 0;
            merSet[setSize].costOfVariant = cost;
            merSet[setSize].costOfFollowing = 0;
            merSet[setSize].followingIndel = VariantType.none;
            setSize++;
        }

        public void CopyfromMer(int fromIdx, int i)
        {
            merSet[i].mer = merSet[fromIdx].mer;
            merSet[i].variantType = merSet[fromIdx].variantType;
            merSet[i].exactMatchesFollowing = merSet[fromIdx].exactMatchesFollowing;
            merSet[i].mismatchesFollowing = merSet[fromIdx].mismatchesFollowing;
            merSet[i].followingBaseIncrement = merSet[fromIdx].followingBaseIncrement;
            merSet[i].gapLength = merSet[fromIdx].gapLength;
            merSet[i].indelBalance = merSet[fromIdx].indelBalance;
            merSet[i].costOfVariant = merSet[fromIdx].costOfVariant;
            merSet[i].costOfFollowing = merSet[fromIdx].costOfFollowing;
            merSet[i].followingIndel = merSet[fromIdx].followingIndel;
        }

        public void CopyfromMerSet(VariantMerSet fromSet)
        {
            if (merSet.Count < fromSet.merSet.Count)
                for (int i = merSet.Count; i < fromSet.merSet.Count; i++)
                    merSet.Add(new Mer());

            for (int i = 0; i < fromSet.setSize; i++)
            {
                merSet[i].mer = fromSet.merSet[i].mer;
                merSet[i].variantType = fromSet.merSet[i].variantType;
                merSet[i].exactMatchesFollowing = fromSet.merSet[i].exactMatchesFollowing;
                merSet[i].mismatchesFollowing = fromSet.merSet[i].mismatchesFollowing;
                merSet[i].followingBaseIncrement = fromSet.merSet[i].followingBaseIncrement;
                merSet[i].gapLength = fromSet.merSet[i].gapLength;
                merSet[i].indelBalance = fromSet.merSet[i].indelBalance;
                merSet[i].costOfVariant = fromSet.merSet[i].costOfVariant;
                merSet[i].costOfFollowing = fromSet.merSet[i].costOfFollowing;
                merSet[i].followingIndel = fromSet.merSet[i].followingIndel;
            }
            setSize = fromSet.setSize;
        }
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
    public class TiledGenome
    {
        public int refCount;
        public MerHashSet? mersFromGenomePlus;
        public MerHashSet? mersFromGenomeMinus;
        public MerHashSet? mersFromGenomeEither;
        public FrozenDictionary<uint, ulong>? seedsIndex;
        public List<ulong> seedsContexts;

        public TiledGenome(int genomesCount)
        {
            this.refCount = genomesCount;
            this.mersFromGenomePlus = null;
            this.mersFromGenomeMinus = null;
            this.mersFromGenomeEither = null;
            this.seedsIndex = null;
            this.seedsContexts = null;
        }
    }
}
