using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using WorkingDogsCore;

namespace MerCollections
{
    // A read-only (once loaded) kMer pair table built on top of Tessel's MerDictionary. These tables are lock-free - and there can be no deletions or resizes while
    // they are being loaded. If any of the partitions get full, any subsequent entries will be put into a per-thread overflow table (these can be resized). 
    // The primary partitions and any overflow tables will be merged by the calling code. 

    public class PairTable
    {
        public int noOfPartitions = 0;                              // no. of mer partitions (hash-distributed mers)
        int maxTableSize = 40000000;                                // ensure table partitions are not too big (used in initial partition sizing)
        int minTableSize = 5000000;                                 // nor too small
        int merSize = 0;

        public const int int31Mask = 0x7fffffff;                    // all but the top bit of an int 32
        public const ulong fullMerMask = 0xffffffffffffffff;        // all bits in the mer are used in comparisons

        public MerDictionary[] repeatedMers = null;                 // hash-partitioned (closed) hash tables holding repeat mers
        // these tables hold canonical (lowest of as-read/RC) mers
        // shared amongst all threads - insertion should be lock free. no deletions or resizes done
        public bool[] repeatedMersFull = null;                      // repeatedMers partition is full, new mers go in per-thread overflow tables
        public MerDictionary[] overflowMers = null;                 // per-thread overflow repeat mer tables

        // compatibility only
        public LowRepMerBuffer culledBuffer = null;                 // temporary shared structures used during multi-threaded low-rep mer flush
        public object culledLock = new object();                    // (only ever done if there are more than 2 seq files being tiled - e.g. multi-lane datasets)
        public List<LowRepMerBuffer> filledCulledBuffers = null;    // only used if culledBuffer becomes full

        // constructor
        public PairTable(long dictionarySize, int noThreads, int merSize)
        {
            this.merSize = merSize;

            // scale genome size to compensate for the number of error-tainted high depth pairs
            dictionarySize = dictionarySize * 4;

            // how many shared mer partitions are needed to safely hold this many distinct k-mers?
            this.noOfPartitions = (int)(dictionarySize / maxTableSize + 1);
            if (this.noOfPartitions < 1)
                this.noOfPartitions = 1;
            // and how big should the partitions be? 
            int partitionSize = (int)(dictionarySize / noOfPartitions);
            if (partitionSize < minTableSize)
                partitionSize = minTableSize;

            repeatedMers = new MerDictionary[noOfPartitions];               // create partitioned dictionaries
            repeatedMersFull = new bool[noOfPartitions];                    // create full flags array (default is false)
            overflowMers = new MerDictionary[noThreads];                    // create per-thread overflow tables

            // initialise per-partition structures
            for (int i = 0; i < noOfPartitions; i++)
            {
                repeatedMers[i] = new MerDictionary(partitionSize, merSize, noThreads);
            }

            // initialise per-thread structures
            for (int i = 0; i < noThreads; i++)
            {
                // overflowMers[i] = new MerDictionary(singletonSize); // allocated when first used to save on memory space
            }
        }

        public bool AddOrIncrement(ulong mer, int threadNo)
        {
            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;

            // look in the main table first - and increment the value if the pair is there
            int idx = repeatedMers[partitionNo].FindEntry(mer);
            if (idx >= 0)
            {
                // minor race here - could result in counts being slightly low sometimes but doesn't matter for the purposes for which the count is used
                // could replace with InterlockedIncrement is this ever needs to be an accurate count
                repeatedMers[partitionNo].values[idx]++;
                //Interlocked.Increment(ref repeatedMers[partitionNo].values[idx]);

                return true;
            }

            // perhaps the pair is in a per-thread overflow table
            if (repeatedMersFull[partitionNo] && overflowMers[threadNo] != null)
            {
                // is the mer in the overflow table?
                int overflowIdx = overflowMers[threadNo].FindEntry(mer);
                if (overflowIdx >= 0)
                {
                    // already there so just add to its count
                    overflowMers[threadNo].values[overflowIdx]++;
                    return true;
                }
            }

            // pair not present so add an entry for it
            if (repeatedMersFull[partitionNo])
            {
                // main table is full, so add the entry to the overflow table for this thread 
                if (overflowMers[threadNo] == null)
                {
                    // no such overflow table yet for this thread, so create one 
                    overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, merSize, 1);
                    //Console.WriteLine("added overflow for thread " + threadNo + " for [" + partitionNo + "]");
                }

                bool full = overflowMers[threadNo].Add(mer, 1);
                // add will always work but could return 'no more please' status so we'll resize in this (unlikely) case
                if (full)
                {
                    //Console.WriteLine("resized overflow for thread " + threadNo + " for [" + partitionNo + "]");
                    overflowMers[threadNo].Resize();
                }
            }
            else
            {
                // space in main table - so add it. Could be a race here where the mer is added twice. Resolved during writing phase.
                bool full = repeatedMers[partitionNo].Add(mer, 1);
                if (full)
                    repeatedMersFull[partitionNo] = true;
            }
            return false;
        }

        public bool Contains(ulong mer, int threadNo)
        {
            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;

            // try the main table
            int idx = repeatedMers[partitionNo].FindEntry(mer);
            if (idx >= 0)
                return true;

            // if not there, try the overflow table if it is in use
            if (repeatedMersFull[partitionNo])
            {
                int overflowIdx = overflowMers[threadNo].FindEntry(mer);
                if (overflowIdx >= 0)
                    return true;
            }

            return false;
        }

        public bool IncrementIfPresent(ulong mer, int threadNo)
        {
            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;

            // first look in the main table
            int idx = repeatedMers[partitionNo].FindEntry(mer);
            if (idx >= 0)
            {
                // minor race here - could result in counts being slightly low sometimes but doesn't matter for the purposes for which the count is used
                // could replace with InterlockedIncrement is this ever needs to be an accurate count
                repeatedMers[partitionNo].values[idx]++;
                //Interlocked.Increment(ref repeatedMers[partitionNo].values[idx]);

                return true;
            }

            // perhaps in an overflow table
            if (repeatedMersFull[partitionNo] && overflowMers[threadNo] != null)
            {
                int overflowIdx = overflowMers[threadNo].FindEntry(mer);
                if (overflowIdx >= 0)
                {
                    // there so just add to its count
                    overflowMers[threadNo].values[overflowIdx]++;
                    return true;
                }
            }

            return false;
        }

        public long Count
        {
            get
            {
                long totalCount = 0;
                for (int p = 0; p < noOfPartitions; p++)
                    totalCount += repeatedMers[p].Count;
                for (int t = 0; t < overflowMers.Length; t++)
                    if (overflowMers[t] != null)
                        totalCount += overflowMers[t].Count;

                return totalCount;
            }
        }

        public long Capacity
        {
            get
            {
                long totalCapacity = 0;
                for (int p = 0; p < noOfPartitions; p++)
                    totalCapacity += repeatedMers[p].Capacity;
                for (int t = 0; t < overflowMers.Length; t++)
                    if (overflowMers[t] != null)
                        totalCapacity += overflowMers[t].Capacity;

                return totalCapacity;
            }
        }

    }

}