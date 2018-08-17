using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Threading;
using System.Runtime.InteropServices;
using WorkingDogsCore;

namespace MerCollections
{
    public class MerTables
    {
        // allow 2 status bits (1 base) at the end of each singleton mer
        private const ulong singletonRCFlagMask = 0x0000000000000001;        // rc flag is lowest order bit
        private const ulong singletonRCFlag = 0x0000000000000001;            // set to 1 if kMer in table is RC form
        private const ulong singletonPromotedFlagMask = 0x0000000000000001;  // same bit re-used to resolve races during promotion
        private const ulong singletonActiveFlagMask = 0x0000000000000002;    // active/promoted flag is 2nd bit (1==active)
        private const ulong singletonMerMask = 0xfffffffffffffffc;           // mask out the last 2 bits (status bits in singleton filter) for comparisons
        private const ulong resetActiveSingletonMask = 0xfffffffffffffffd;   // just clear the active bit 
        private const int int31Mask = 0x7fffffff;                            // all but the top bit of an int 32

        private int noThreads = 0;                                  // no of parallel calling threads
        const int noSingletonFlushingThreads = 2;                   // no. of parallel flushing threads
        private int noOfPartitions = 0;                             // no. of mer partitions (hash-distributed mers)
        private int noSingletonPartitions = 0;                      // base-prefix ordered entries
        private int singletonPrefixBits = 0;                        // how many LHS bits are used for partitioning singletons
        int[] maxSingletonCapacity;                                 // max no. of entries allowed in each singleton table before flushing
        bool flushSingletons;                                       // flush full singleton tables or leave them in the deferred list?
        bool inexactCount;                                          // discard singleton tables rather than keeping them for later merging
        int maxTableSize = 40000000;                                // ensure table partitions are not too big (used in initial partition sizing)
        int minTableSize = 5000000;                                 // nor too small
        int[] possibleSingletonSizes = new int[] { 10000000, 20000000, 30000000, 40000000, 50000000 }; // quantised singleton table sizes
        const long minKeepDepth = 3;                                // only keep k-mers of at least this depth in the repeated kmers table at the end-of-file flush
        int merSize = 0;                                            // set during initialisation and passed to the MerCollection constructors
        bool canonical;                                             // save canonical or as-read form of kMer

        public MerDictionary[] repeatedMers = null;                 // hash-partitioned (closed) hash tables holding repeat mers
                                                                    // these tables hold canonical (lowest of as-read/RC) mers
                                                                    // shared amongst all threads - insertion should be lock free. no deletions or resizes done
        public bool[] repeatedMersFull = null;                      // repeatedMers partition is full, new mers go in per-thread overflow tables
        public MerDictionary[] overflowMers = null;                 // per-thread overflow repeat mer tables

        public LowRepMerBuffer culledBuffer = null;                 // temporary shared structures used during multi-threaded low-rep mer flush
        public object culledLock = new object();                    // (only ever done if there are more than 2 seq files being tiled - e.g. multi-lane datasets)
        public List<LowRepMerBuffer> filledCulledBuffers = null;    // only used if culledBuffer becomes full

        public SingletonCollection[] singletonFilters = null;       // partitioned hash sets holding possible singleton mers
                                                                    // these tables contain partitioned mers with an RC status bit at the RHS of the ulong
                                                                    // [0] says whether the mer is the RC form (true) or the as-is (false) form (used to calculate initial counts)
        public List<SingletonCollection>[] singletonsWaitingFlush = null; // (probably largely empty) list of singleton buffers waiting to be flushed. Used to consolidate flush files
        public List<string>[] flushedSingletonFNs = null;           // per-partition list of file names of flushed singleton partitions
        public List<ulong>[] firstFlushedSingletonMer = null;       // lowest mer in each flushed singleton file
        public List<ulong>[] lastFlushedSingletonMer = null;        // highest mer in each flushed singleton file
        public List<int>[] flushedSingletonsCount = null;           // number of mers in each flushed singleton file
        public int[] flushSingletonNumber = null;                   // next flush number for each partition (used to construct unique file names)
        public List<string> flushedLowRepsFNs = null;               // list of file names of flushed LowReps partitions
        public List<ulong> firstFlushedLowRepsMer = null;           // lowest mer in each flush LowReps file
        public List<ulong> lastFlushedLowRepsMer = null;            // highest mer in each flush LowReps file
        public List<int> flushedLowRepsCount = null;                // number of mers in each flushed LowReps file
        public int[] flushLowRepsNumber = null;                     // next flush number for each partition (used to construct unique file names)
        string tempDirectory = "";                                  // where the temporary flush files are to be written

        private List<SingletonCollection> singletonPool;            // pool of recyclable singleton tables. Singleton tables are placed here after they've been flushed. 

        // event for signalling singleton writer that a buffer is queued 
        EventWaitHandle singletonWriterEWH = new EventWaitHandle(false, EventResetMode.AutoReset);
        // and the queue of singleton buffers waiting to be written
        Queue<SingletonCollection> singletonsAwaitingFlush = new Queue<SingletonCollection>();
        // and the flushing threads
        Thread[] singletonFlushingThreads;

        //const bool tracing = false;
        //public Queue<TraceEntry> traceUpdates = new Queue<TraceEntry>(1000000);
        //const int maxTrace = 1000000;

        // constructor
        public MerTables(long dictionarySize, int merSize, string tempDir, int noThreads, bool flushSingletons, bool inexactCount, bool canonical)
        {
            this.noThreads = noThreads;

            // scale genome size to compensate for the number of repeated error k-mers
            dictionarySize = dictionarySize * 2;
            // and double this if we're counting non-canonical kMers as each kMer in the genome will be counted twice
            if (!canonical)
                dictionarySize = dictionarySize * 2;

            // how many shared mer partitions are needed to safely hold this many distinct k-mers?
            this.noOfPartitions = (int)(dictionarySize / maxTableSize + 1);
            if (this.noOfPartitions < 1)
                this.noOfPartitions = 1;
            // and how big should the partitions be? 
            int partitionSize = (int)(dictionarySize / noOfPartitions);
            if (partitionSize < minTableSize)
                partitionSize = minTableSize;

            // how many singleton partitions are desirable?
            // per-partition singleton tables can use of a considerable amount of memory, but too few will increase the number of concurrently open files during merge
            noSingletonPartitions = 32;
            singletonPrefixBits = 5;

            if (noOfPartitions <= 32)
            {
                noSingletonPartitions = 16;
                singletonPrefixBits = 4;
            }
            if (noOfPartitions <= 16)
            {
                noSingletonPartitions = 4;
                singletonPrefixBits = 2;
            }
            //if (noOfPartitions <= 4)
            //{
            //    noSingletonPartitions = 1;
            //    singletonPrefixBits = 0;
            //}

            this.tempDirectory = tempDir;
            this.flushSingletons = flushSingletons;
            this.inexactCount = inexactCount;
            this.merSize = merSize;
            this.canonical = canonical;

            repeatedMers = new MerDictionary[noOfPartitions];               // create partitioned dictionaries
            repeatedMersFull = new bool[noOfPartitions];                    // create full flags array (default is false)
            overflowMers = new MerDictionary[noThreads];                    // create per-thread overflow tables
            singletonFilters = new SingletonCollection[noSingletonPartitions];    // create partitioned singleton filters
            singletonsWaitingFlush = new List<SingletonCollection>[noSingletonPartitions];
            flushedSingletonFNs = new List<string>[noSingletonPartitions];
            firstFlushedSingletonMer = new List<ulong>[noSingletonPartitions];
            lastFlushedSingletonMer = new List<ulong>[noSingletonPartitions];
            flushedSingletonsCount = new List<int>[noSingletonPartitions];
            flushSingletonNumber = new int[noSingletonPartitions];
            maxSingletonCapacity = new int[noSingletonPartitions];

            flushedLowRepsFNs = new List<string>();
            firstFlushedLowRepsMer = new List<ulong>();
            lastFlushedLowRepsMer = new List<ulong>();
            flushedLowRepsCount = new List<int>();

            singletonPool = new List<SingletonCollection>(noThreads);

            // initialise per-partition structures
            for (int i = 0; i < noOfPartitions; i++)
            {
                repeatedMers[i] = new MerDictionary(partitionSize, merSize, noThreads);
            }

            // initialise per-singleton-partition structures
            for (int i = 0; i < noSingletonPartitions; i++)
            {
                int scaledSingletonSize = /*2 **/ partitionSize / noSingletonPartitions * (noSingletonPartitions - i);
                scaledSingletonSize = RoundUpSingletonSize(scaledSingletonSize, possibleSingletonSizes);
                singletonFilters[i] = new SingletonCollection(scaledSingletonSize, merSize, singletonMerMask, singletonActiveFlagMask, singletonRCFlagMask, i, singletonPrefixBits);
                //long safeMaxSingletonCapacity = (long)scaledSingletonSize * 95 / 100;   // careful of int overflow
                //maxSingletonCapacity[i] = (int)(safeMaxSingletonCapacity);
                maxSingletonCapacity[i] = scaledSingletonSize - noThreads;                         
                flushedSingletonFNs[i] = new List<string>();
                flushSingletonNumber[i] = 1;
                firstFlushedSingletonMer[i] = new List<ulong>();
                lastFlushedSingletonMer[i] = new List<ulong>();
                flushedSingletonsCount[i] = new List<int>();
            }

            // initialise per-thread structures
            for (int i = 0; i < noThreads; i++)
            {
                // overflowMers[i] = new MerDictionary(singletonSize); // allocated when first used to save on memory space
            }

            // and start the flushing threads
            singletonFlushingThreads = new Thread[noSingletonFlushingThreads];
            FlushingSingletonsThreadParams flushingParams = new FlushingSingletonsThreadParams();
            flushingParams.buffersToBeFlushed = singletonsAwaitingFlush;
            flushingParams.flushingEWH = singletonWriterEWH;
            for (int t = 0; t < noSingletonFlushingThreads; t++)
            {
                singletonFlushingThreads[t] = new Thread(new ParameterizedThreadStart(SingletonFlusherThread));
                singletonFlushingThreads[t].Priority = ThreadPriority.Normal;
                singletonFlushingThreads[t].Name = "SingletonFlusher";
                singletonFlushingThreads[t].Start(flushingParams);
            }
        }

        // ensure MerTable is in a final, stable state (called at the end of the counting phase)
        public void Stabilise()
        {
            //if (singletonsAwaitingFlush.Count > 0)
            //    Console.WriteLine("waiting for " + singletonsAwaitingFlush.Count + " singleton buffers to be flushed");
            SingletonCollection lastBuffer = GetSingletonTable(singletonFilters[0]);
            lastBuffer.lastBuffer = true;
            lock (singletonsAwaitingFlush)
            {
                for (int t = 0; t < noSingletonFlushingThreads; t++)
                    singletonsAwaitingFlush.Enqueue(lastBuffer);
            }
            singletonWriterEWH.Set();

            // wait for the flushing threads to finish
            for (int t = 0; t < noSingletonFlushingThreads; t++)
            {
                singletonFlushingThreads[t].Join();
                singletonFlushingThreads[t] = null;
            }
            //Console.WriteLine("singleton flushing threads finished");
        }

        // wait for all pending singleton buffers to be written
        public void WaitForSingletonBufferFlush()
        {
            //if (singletonsAwaitingFlush.Count > 0)
            //    Console.WriteLine("waiting for " + singletonsAwaitingFlush.Count + " singleton buffers to be flushed");

            while (singletonsAwaitingFlush.Count > 0)
                Thread.Sleep(1000);
        }

        public void AddOrIncrement(ulong mer, int merSize, int threadNo)
        {
            long addingIncrement = 0x0000000100000000;                  // assume we've got the as-read form, not the canonical form
            ulong rcFlagToBeSet = 0x0;                                  // and that we don't want to set the RC flag
            //if (mer == 0xABF70D5BD9C54E7C)
            //    Debugger.Break();
            // generate canonical k-mer first
            ulong rcMer = kMers.ReverseComplement(mer, merSize);
            if (canonical && rcMer < mer)
            {
                mer = rcMer;
                addingIncrement = 0x0000000000000001;                   // increment the low part of the count pair
                rcFlagToBeSet = singletonRCFlag;                        // remember if the canonical k-mer was the RC form
            }

            int absMerHashCode = mer.GetHashCode() & int31Mask;
            int partitionNo = absMerHashCode % noOfPartitions;
            int singletonPartitionNo = (int)(mer >> (64 - singletonPrefixBits));

            // this mer may have been seen before, so first try updating it in one of the repeated mer tables
            bool updatedRepeat = UpdateRepeatedMer(partitionNo, mer, threadNo, mer, addingIncrement);

            if (updatedRepeat)
                return;

            // handling a k-mer for the first time - try adding it to the singletons table
            // ----------------------------------------------------------------------------

            // get a stable pointer to the current singletons table (in case someone else fills it and initiates a flush while we're still busy with it)
            SingletonCollection thisSingletonPartition = singletonFilters[singletonPartitionNo];
            //Interlocked.Increment(ref thisSingletonPartition.activeCount);

            // try to add this mer to this partition's singletons collection (and fetch the existing singleton+flag if it's already there)
            int filterIdx;
            ulong fMer = (mer << singletonPrefixBits) | rcFlagToBeSet | singletonActiveFlagMask;
            bool added = thisSingletonPartition.TryInsertKey(fMer, out filterIdx);

            if (added)
            {
                // successfully added this mer so we must be seeing it for the first time 

                // if singleton table is now full, flush it out and empty the table
                if (thisSingletonPartition.Count >= maxSingletonCapacity[singletonPartitionNo])
                {
                    bool flushNeeded = true;
                    int flushNumberToUse = 0;

                    // lock this section to avoid two threads trying to flush/replace the same singleton buffer concurrently
                    lock (thisSingletonPartition)
                    {
                        // test entry condition now that we have the lock (filter may have been reset while we were waiting)
                        if (!thisSingletonPartition.flushed)
                        {
                            // allocate a replacement table for the other threads to use while we're flushing this one
                            SingletonCollection emptySingletonFilter = GetSingletonTable(thisSingletonPartition);  // allocate new local filter for the partition
                            maxSingletonCapacity[singletonPartitionNo] = emptySingletonFilter.length - noThreads;               // use whatever size it happened to be
                            singletonFilters[singletonPartitionNo] = emptySingletonFilter;                                      // make it visible to the concurrent threads (atomic assignment)
                            //Console.WriteLine("replacement singletons: " + maxSingletonCapacity[singletonPartitionNo] + "/" + singletonFilters[singletonPartitionNo].Capacity);
                            thisSingletonPartition.flushed = true;
                            flushNumberToUse = flushSingletonNumber[singletonPartitionNo];
                            flushSingletonNumber[singletonPartitionNo]++;
                        }
                        else
                            flushNeeded = false;
                    }

                    if (flushNeeded)
                    {
                        thisSingletonPartition.flushNo = flushNumberToUse;

                        if (inexactCount)
                        {
                            // just jellyfish-like inexact counts that ignore singletons
                            ReturnSingletonTable(thisSingletonPartition);
                        }
                        else
                        {
                            thisSingletonPartition.timeWhenQueued = DateTime.Now;
                            lock (singletonsAwaitingFlush)
                            {
                                singletonsAwaitingFlush.Enqueue(thisSingletonPartition);
                            }
                            singletonWriterEWH.Set();
                            //Console.WriteLine("queued " + singletonPartitionNo + "/" + flushNumberToUse + " for flushing");
                        }
                    }
                }
            }
            else
            {
                // Insert failed, so must be seeing this k-mer for second (or rarely more) time. Mark as inactive in singletons and add to a repeats table with appropriate counts. 
                // There can be a race here with two threads trying to concurrently promote the same singleton. This is resolved by atomically clearing the singleton
                // active flag - and only one of the threads will get the 'active' flag returned from the Exchange. This thread does the promotion - and then sets the 
                // promotion-complete bit for the singleton. The other threads will spin until they find this bit has been set.

                //if (tracing)
                //    lock (traceUpdates)
                //    {
                //        traceUpdates.Enqueue(new TraceEntry(threadNo, 1, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.entries[filterIdx].key));
                //        //traceUpdates.Enqueue(new TraceEntry(threadNo, 1, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.keys[filterIdx]));            
                //        if (traceUpdates.Count > maxTrace)
                //        {
                //            traceUpdates.Dequeue();
                //        }
                //    }

                // get the current value of this singleton entry (safe because the promotion changes are progressive)
                //ulong merFromFilter = (ulong)thisSingletonPartition.entries[filterIdx].key;
                ulong merFromFilter = (ulong)thisSingletonPartition.keys[filterIdx];         
                // and see if this singleton may have already been promoted
                bool activeSingleton = (merFromFilter & singletonActiveFlagMask) != 0;

                // if this singleton may be 'active', try to promote it
                if (activeSingleton)
                {
                    ulong inactiveMer = mer & singletonMerMask;                      // build what the inactive-but-being-promoted entry should look like
                    // if no-one else has altered the singleton entry, then set it to inactive-but-being-promoted ('long' for CEX requirements)
                    //long currentMerFromFilter = Interlocked.CompareExchange(ref thisSingletonPartition.entries[filterIdx].key, (long)inactiveMer, (long)merFromFilter);
                    long currentMerFromFilter = Interlocked.CompareExchange(ref thisSingletonPartition.keys[filterIdx], (long)inactiveMer, (long)merFromFilter);

                    //if (tracing)
                    //    lock (traceUpdates)
                    //    {
                    //        traceUpdates.Enqueue(new TraceEntry(threadNo, 2, singletonPartitionNo, filterIdx, (ulong)currentMerFromFilter));
                    //        if (traceUpdates.Count > maxTrace)
                    //        {
                    //            traceUpdates.Dequeue();
                    //        }
                    //    }

                    // if this thread successfully set the singleton to 'inactive', it will take care of the promotion
                    if (currentMerFromFilter == (long)merFromFilter)
                    {
                        ulong rcFlag = merFromFilter & singletonRCFlagMask;          // non-zero --> RC found in singletons 

                        long initialCount = 0;
                        if (rcFlag != 0)	                                // singleton was seen in RC form
                            initialCount = 0x0000000000000001;
                        else	                                            // singleton was seen in as-is form
                            initialCount = 0x0000000100000000;

                        if (repeatedMersFull[partitionNo])
                        {
                            if (overflowMers[threadNo] == null)
                            {
                                overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, merSize, 1);
                                //Console.WriteLine("added overflow for thread " + threadNo + " for [" + partitionNo + "]");
                            }

                            bool full = overflowMers[threadNo].Add(mer, initialCount);
                            if (full)
                                overflowMers[threadNo].Resize();
                        }
                        else
                        {
                            bool full = repeatedMers[partitionNo].Add(mer, initialCount);
                            if (full)
                                repeatedMersFull[partitionNo] = true;
                        }

                        // now that the mer has been promoted, set the 'promoted' flag
                        inactiveMer = inactiveMer | (long)singletonPromotedFlagMask;
                        //thisSingletonPartition.entries[filterIdx].key = (long)inactiveMer;
                        thisSingletonPartition.keys[filterIdx] = (long)inactiveMer;

                        //if (tracing)
                        //    lock (traceUpdates)
                        //    {
                        //        traceUpdates.Enqueue(new TraceEntry(threadNo, 3, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.entries[filterIdx].key));
                        //        //traceUpdates.Enqueue(new TraceEntry(threadNo, 3, singletonPartitionNo, filterIdx, (ulong)thisSingletonPartition.keys[filterIdx]));        
                        //        if (traceUpdates.Count > maxTrace)
                        //        {
                        //            traceUpdates.Dequeue();
                        //        }
                        //    }
                    }
                }

                // singleton is now known to be no longer active, so wait (if necessary) for the 'promoted' flag to be set and increment the repeat counter

                //merFromFilter = (ulong)thisSingletonPartition.entries[filterIdx].key;
                merFromFilter = (ulong)thisSingletonPartition.keys[filterIdx];

                //if (tracing)
                //    lock (traceUpdates)
                //    {
                //        traceUpdates.Enqueue(new TraceEntry(threadNo, 4, singletonPartitionNo, filterIdx, merFromFilter));
                //        if (traceUpdates.Count > maxTrace)
                //        {
                //            traceUpdates.Dequeue();
                //        }
                //    }

                bool promotionComplete = (merFromFilter & singletonPromotedFlagMask) != 0;
                bool alreadySlept = false;
                while (!promotionComplete)
                {
                    //promotionComplete = (((ulong)thisSingletonPartition.entries[filterIdx].key & singletonPromotedFlagMask) != 0);
                    promotionComplete = (((ulong)thisSingletonPartition.keys[filterIdx] & singletonPromotedFlagMask) != 0);           
                    if (alreadySlept && !promotionComplete)
                    {
                        //if (tracing)
                        //{
                        //    lock (traceUpdates)
                        //    {
                        //        StreamWriter trace = new StreamWriter("trace.txt");
                        //        foreach (TraceEntry t in traceUpdates)
                        //            trace.WriteLine(t.place + "\t" + t.thread + "\t" + t.partition + "\t" + t.index + "\t" + t.value.ToString("x16"));
                        //        trace.Close();
                        //    }
                        //    Console.WriteLine("promotion still not complete after sleep");
                        //}
                    }
                    if (!promotionComplete)
                        Thread.Sleep(100);
                    alreadySlept = true;
                }

                UpdateRepeatedMerAfterPromotion(partitionNo, mer, threadNo, mer, addingIncrement);
                //if (!updateSucceeded)
                //{
                //    lock (traceUpdates)
                //    {
                //        StreamWriter trace = new StreamWriter("trace.txt");
                //        foreach (TraceEntry t in traceUpdates)
                //            trace.WriteLine(t.thread + "\t" + t.place + "\t" + t.partition + "\t" + t.index + "\t" + t.value.ToString("x16"));
                //        trace.Close();
                //    }
                //    Console.WriteLine("UpdateRepeatedMerRetry failed after waiting for promotion to complete");
                //}

            }

            //Interlocked.Decrement(ref thisSingletonPartition.activeCount);
        }

        private bool UpdateRepeatedMer(int partitionNo, ulong pMer, int threadNo, ulong mer, long addingIncrement)
        {
            bool updated = false;

            // update the count for this mer if it's in the shared repeated set 
            updated = repeatedMers[partitionNo].UpdateIfPresent(pMer, addingIncrement);
            if (updated)
            {
                //repeats++;
                return true;
            }

            // mer could be in the thread-local overflow table
            updated = overflowMers[threadNo] != null && overflowMers[threadNo].UpdateIfPresent(mer, addingIncrement);
            if (updated)
            {
                //overflow++;
                return true;
            }

            return false;
        }

        private void UpdateRepeatedMerAfterPromotion(int partitionNo, ulong pMer, int threadNo, ulong mer, long addingIncrement)
        {
            bool updated = false;

            // update the count for this mer if it's in the shared repeated set 
            updated = repeatedMers[partitionNo].UpdateIfPresent(pMer, addingIncrement);
            if (updated)
            {
                //repeats++;
                return;
            }

            // mer could be in the thread-local overflow table (or it may have been promoted by another thread)
            updated = overflowMers[threadNo] != null && overflowMers[threadNo].UpdateIfPresent(mer, addingIncrement);
            if (updated)
            {
                //overflow++;
                return;
            }

            // if we get here, the promoted singleton must have gone to another thread's overflow table (and the partition must be full)
            // so create a new overflow table if we don't already have one
            if (repeatedMersFull[partitionNo] && overflowMers[threadNo] == null)
            {
                overflowMers[threadNo] = new MerDictionary(repeatedMers[partitionNo].lengthEntries / 10, merSize, 1); 
                //Console.WriteLine("added overflow for thread " + threadNo);
            }

            // finally add this promoted singleton to the overflow table for this thread 
            if (overflowMers[threadNo] != null)
            {
                bool overflowFull = overflowMers[threadNo].Add(mer, addingIncrement);
                if (overflowFull)
                    overflowMers[threadNo].Resize();
                return;
            }

            // the shared repeated table isn't full but we didn't find the singleton there, could have been a race on a bucket that resulted in an orphan
            repeatedMers[partitionNo].Add(pMer, addingIncrement);

        }

        private void SingletonFlusherThread(object threadParams)
        {
            FlushingSingletonsThreadParams theseParams = (FlushingSingletonsThreadParams)threadParams;
            EventWaitHandle flushingEWH = theseParams.flushingEWH;
            Queue<SingletonCollection> buffersToBeFlushed = theseParams.buffersToBeFlushed;           
            
            bool stopFlusher = false;

            // loop until we've flushed the final buffer
            while (!stopFlusher)
            {
                SingletonCollection singletons = null;

                // get the next buffer to be flushed 
                lock (buffersToBeFlushed)
                {
                    if (buffersToBeFlushed.Count > 0)
                        singletons = buffersToBeFlushed.Dequeue();
                }

                // wait for something to be enqueued on the flushing queue if it was empty
                if (singletons == null)
                {
                    //Console.WriteLine("flusher: sleeping");
                    flushingEWH.WaitOne(30000);
                    //Console.WriteLine("flusher: woken");
                    continue;
                }

                // check if this is the specially marked 'last buffer' - and exit if it is
                if (singletons.lastBuffer)
                {
                    stopFlusher = true;
                    singletons = null;
                    // propagate 'something waiting' so the next thread can pick up its terminating message as well
                    flushingEWH.Set();
                    break;
                }

                byte[] byteBuffer = singletons.byteBuffer;

                int partitionNo = singletons.partitionNo;
                int flushNo = singletons.flushNo;

                // pause briefly here if necessary to allow any in-flight updates to this table to settle
                if ((DateTime.Now - singletons.timeWhenQueued).TotalMilliseconds < 1000.0)
                    Thread.Sleep(1000);

                //Console.WriteLine("flushing singletons[" + partitionNo + ", " + flushNo + "]" + " " + buffersToBeFlushed.Count + " buffers still queued");

                // condense and in-place sort the singleton mers (left as un-extended)
                singletons.Condense();
                singletons.Sort();

                // if we're only writing out few active singletons, save these for the next flush by adding the singletons set to the waiting list 
                // and we'll also take this path if we've been asked not to flush singletons at all

                int pendingSingletonsCount = singletons.Count;
                if (singletonsWaitingFlush[partitionNo] != null)
                    foreach (SingletonCollection pending in singletonsWaitingFlush[partitionNo])
                        pendingSingletonsCount += pending.Count;
                if (!this.flushSingletons || pendingSingletonsCount < singletons.Capacity / 2)
                {
                    if (singletonsWaitingFlush[partitionNo] == null)
                        singletonsWaitingFlush[partitionNo] = new List<SingletonCollection>();

                    singletonsWaitingFlush[partitionNo].Add(singletons);

                    //Console.WriteLine("deferring singleton write[" + partitionNo + "] " + merIdx + " mers");
                    continue;
                }

                // get a pointer to the singletons themselves (from Entries)
                //MerCollection.Entry[] sortedSingletons = singletons.entries;;
                long[] sortedSingletons = singletons.keys;
                int merIdx = singletons.Count;

                // got a singleton set worth flushing so add any pending singleton sets to it 
                if (singletonsWaitingFlush[partitionNo] != null && singletonsWaitingFlush[partitionNo].Count != 0)
                {
                    // find out how many singletons are waiting (across all waiting sets)
                    int sumWaiting = 0;
                    foreach (SingletonCollection s in singletonsWaitingFlush[partitionNo])
                        sumWaiting += s.Count;

                    // ensure the buffer is big enough
                    if (merIdx + sumWaiting > singletons.Capacity)
                        //Array.Resize<MerCollection.Entry>(ref sortedSingletons, merIdx + sumWaiting);
                        Array.Resize<long>(ref sortedSingletons, merIdx + sumWaiting);

                    // copy the waiting singletons into this buffer
                    for (int i = 0; i < singletonsWaitingFlush[partitionNo].Count; i++)
                    {
                        int singlesInArray = singletonsWaitingFlush[partitionNo][i].Count;
                        //Console.WriteLine("writing pending singletons[" + partitionNo + ", " + i + "] with " + singlesInArray + " kMers");
                        //MerCollection.Entry[] pendingSingles = singletonsWaitingFlush[partitionNo][i].entries;
                        long[] pendingSingles = singletonsWaitingFlush[partitionNo][i].keys;
                        for (int s = 0; s < singlesInArray; s++)
                        {
                            sortedSingletons[merIdx] = pendingSingles[s];
                            merIdx++;
                        }
                    }

                    singletonsWaitingFlush[partitionNo].Clear();

                    // and re-sort the combined singleton entries
                    singletons.Sort();
                }

                //ulong previousMer = 0;

                // and write out the (possibly augmented) set of singletons     
                if (merIdx > 0)
                {
                    //Console.WriteLine("wrote  " + merIdx + " singletons for flush " + partitionNo + "-" + flushNo + " " + buffersToBeFlushed.Count + " still queued");

                    // Add process ID at the beginning of the file name
                    int processID = Process.GetCurrentProcess().Id;
                    string binaryfileName = tempDirectory + processID + "_" + partitionNo + "_" + flushNo + ".bst";

                    bool writeSucceeded = WriteSingletonFile(binaryfileName, singletons);
                    if (!writeSucceeded)
                    {
                        Console.WriteLine("retrying WriteSingletonFile after exception");
                        Thread.Sleep(1000);
                        writeSucceeded = WriteSingletonFile(binaryfileName, singletons);
                    }

                    //Console.WriteLine("wrote " + binaryfileName);
                    flushedSingletonFNs[partitionNo].Add(binaryfileName);                                       // Add name of file to list of flush files for the partition
                    //firstFlushedSingletonMer[partitionNo].Add((ulong)sortedSingletons[0].key & singletonMerMask);           // remember lowest and highest mer in each flush file 
                    //lastFlushedSingletonMer[partitionNo].Add((ulong)sortedSingletons[merIdx - 1].key & singletonMerMask);   // (without the RC bit)
                    firstFlushedSingletonMer[partitionNo].Add((ulong)sortedSingletons[0] & singletonMerMask);           // remember lowest and highest mer in each flush file 
                    lastFlushedSingletonMer[partitionNo].Add((ulong)sortedSingletons[merIdx - 1] & singletonMerMask);   // (without the RC bit)   
                    flushedSingletonsCount[partitionNo].Add(merIdx);                                            // and how many mers we just wrote
                }

                ReturnSingletonTable(singletons);

            } // end of thread loop
        }

        private bool WriteSingletonFile(string binaryfileName, SingletonCollection singletons)
        {
            bool writeSucceeded = true;
            BinaryWriter flushWriter = null;

            long[] sortedSingletons = singletons.keys;
            byte[] byteBuffer = singletons.byteBuffer;

            try
            {
                int byteIdx = 0;
                flushWriter = new BinaryWriter(new FileStream(binaryfileName, FileMode.Create, FileAccess.Write, FileShare.None, 100000));

                // Write the number of k-mers in the flush file at the start
                flushWriter.Write(singletons.Count);
                // Write out the in-use singletons
                for (int i = 0; i < singletons.Count; i++)
                {
                    //flushWriter.Write((ulong)sortedSingletons[i].key);
                    //flushWriter.Write((ulong)sortedSingletons[i]);

                    bool rcMer;
                    ulong singleton = singletons.Extend((ulong)sortedSingletons[i], out rcMer);

                    byteBuffer[byteIdx++] = (byte)singleton;
                    byteBuffer[byteIdx++] = (byte)(singleton >> 8);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 16);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 24);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 32);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 40);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 48);
                    byteBuffer[byteIdx++] = (byte)(singleton >> 56);
                    byteBuffer[byteIdx++] = Convert.ToByte(rcMer);

                    if (byteIdx == SingletonCollection.byteBufferLen)
                    {
                        flushWriter.Write(byteBuffer, 0, byteIdx);
                        byteIdx = 0;
                    }

                    //if (singleton < previousMer)
                    //    Console.WriteLine("sorted singleton failure");
                    //previousMer = singleton;
                }
                flushWriter.Write(byteBuffer, 0, byteIdx);

                flushWriter.Close();
            }
            catch (Exception e)
            {
                Exception baseException = e.GetBaseException();
                Console.WriteLine("WriteSingletonFile exception: " + e.Message + " " + e.HResult.ToString("X16") + ", base exception " + baseException.HResult.ToString("X16"));
                writeSucceeded = false;
                if (flushWriter != null)
                    flushWriter.Close();
            }


            return writeSucceeded;
        }

        // flush the low-rep mers from the repeat tables, condense the remaining repeated mers and fold in the per-thread repeats. Can only be called after all the
        // threads have finished for a seq data file. This code is *not* thread-safe.
        public void FlushLowRepMers(MerTables merTable, int fileNo)
        {
            // allocate a buffer to hold the flushed low-rep mers
            //int initialBufferLength = 500000;
            int initialBufferLength = this.repeatedMers[0].Capacity;
            culledBuffer = new LowRepMerBuffer();
            culledBuffer.keys = new ulong[initialBufferLength + noOfPartitions];
            culledBuffer.values = new long[initialBufferLength + noOfPartitions];
            culledBuffer.idx = 0;
            culledBuffer.bufferActive = true;
            culledBuffer.bufferNo = 1;
            culledBuffer.limit = initialBufferLength;
            culledLock = new object();

            FlushingMerThreadParams[] flushingParams = new FlushingMerThreadParams[noOfPartitions];
            Thread[] flushingThreads = new Thread[noOfPartitions];

            for (int p = 0; p < noOfPartitions; p++)
            {
                flushingParams[p] = new FlushingMerThreadParams();
                flushingParams[p].merTable = merTable;
                flushingParams[p].partitionNo = p;
                flushingThreads[p] = new Thread(new ParameterizedThreadStart(MerTables.FlushLowRepMersInPartition));
                flushingThreads[p].Priority = ThreadPriority.BelowNormal;
                flushingThreads[p].Start(flushingParams[p]);
            }

            for (int p = 0; p < noOfPartitions; p++)
            {
                flushingThreads[p].Join();
                flushingThreads[p] = null;
            }

            // write out any filled culled buffers
            int bufferNo = 0;
            if (filledCulledBuffers != null)
            {
                for (int i = 0; i < filledCulledBuffers.Count; i++)
                {
                    WriteLowRepMers(fileNo, bufferNo, filledCulledBuffers[i], filledCulledBuffers[i].keys.Length);
                    bufferNo++;
                    filledCulledBuffers[i] = null;
                }
                filledCulledBuffers = null;
            }
            // finally write out any remaining culled low-rep mers
            if (culledBuffer.idx > 0)
                WriteLowRepMers(fileNo, bufferNo, culledBuffer, culledBuffer.idx);

            // return the temporary buffers
            culledBuffer = null;

            // finally push the per-thread dictionaries to the shared dictionary
            for (int t = 0; t < overflowMers.Length; t++)
            {
                if (overflowMers[t] == null)
                    continue;

                MerDictionary currentOverflow = overflowMers[t];
                MerDictionary replacementOverflow = new MerDictionary(currentOverflow.Capacity, merSize, 1);

                foreach (KeyValuePair<ulong, long> kvp in currentOverflow)
                {
                    int absMerHashCode = kvp.Key.GetHashCode() & int31Mask;
                    int partitionNo = absMerHashCode % noOfPartitions;

                    if (repeatedMersFull[partitionNo])
                        replacementOverflow.Add(kvp.Key, kvp.Value);
                    else
                    {
                        bool OK = repeatedMers[partitionNo].Add(kvp.Key, kvp.Value);
                        if (!OK)
                            repeatedMersFull[partitionNo] = true;
                    }
                }

                overflowMers[t] = replacementOverflow;
            }

        }

        private static void FlushLowRepMersInPartition(object threadParams)
        {
            FlushingMerThreadParams flushingParams = (FlushingMerThreadParams)threadParams;
            MerTables merTable = flushingParams.merTable;
            int partitionNo = flushingParams.partitionNo;
            int culledMers = merTable.repeatedMers[partitionNo].Reduce(minKeepDepth, merTable);
            // and it can't possibly be full if we've just culled some entries for it
            if (culledMers > 0)
                merTable.repeatedMersFull[partitionNo] = false;
            //Console.WriteLine("flushed " + culledMers + " from partition " + partitionNo + ". " + merTable.repeatedMers[partitionNo].Count + " left");
        }

        private void WriteLowRepMers(int fileNo, int bufferNo, LowRepMerBuffer culledBuffer, int noMers)
        {
            ulong[] culledBufferKeys = culledBuffer.keys;
            long[] culledBufferValues = culledBuffer.values;

            //for (int i = 0; i < noMers; i++)
            //    if (culledBufferValues[i] == 0)
            //        Debugger.Break();

            //Array.Sort<ulong, long>(culledBufferKeys, culledBufferValues, 0, noMers);
            MerArraySorter<long> merSorter = new MerArraySorter<long>(culledBufferKeys, culledBufferValues);
            merSorter.Sort(0, noMers - 1);
            //for (int i = 0; i < noMers - 2; i++)
            //    if (culledBufferKeys[i] > culledBufferKeys[i + 1])
            //    {
            //        Console.WriteLine("sort failed");
            //        break;
            //    }

            // Add process ID at the beginning of the file name
            int processID = Process.GetCurrentProcess().Id;
            string binaryfileName = tempDirectory + processID + "_" + fileNo + "_" + bufferNo + ".bfm";

            bool writeSucceeded = WriteLowRepMersFile(binaryfileName, noMers, culledBufferKeys, culledBufferValues);

            if (!writeSucceeded)
            {
                Console.WriteLine("retrying WriteLowRepMersFile after exception");
                Thread.Sleep(1000);
                writeSucceeded = WriteLowRepMersFile(binaryfileName, noMers, culledBufferKeys, culledBufferValues);
            }

            flushedLowRepsFNs.Add(binaryfileName);                       // Add name of file to list of flush files for the partition
            firstFlushedLowRepsMer.Add(culledBufferKeys[0]);             // remember lowest and highest mer in each flush file 
            lastFlushedLowRepsMer.Add(culledBufferKeys[noMers - 1]);
            flushedLowRepsCount.Add(noMers);                             // and how many mers we wrote

            //Console.WriteLine("flushed " + noMers + " low rep mers (" + fileNo + "_" + bufferNo + ")");
        }

        private bool WriteLowRepMersFile(string binaryfileName, int noMers, ulong[] culledBufferKeys, long[] culledBufferValues)
        {
            bool writeSucceeded = true;
            BinaryWriter flushWriter = null;

            try
            {            
                // Binary writer
                flushWriter = new BinaryWriter(new FileStream(binaryfileName, FileMode.Create, FileAccess.Write, FileShare.None, 100000));

                // Write the number of k-mers in the flush file at the start
                flushWriter.Write(noMers);
                // Write out the in-use singletons
                for (int i = 0; i < noMers; i++)
                {
                    flushWriter.Write(culledBufferKeys[i]);
                    flushWriter.Write((ulong)culledBufferValues[i]);
                }
                flushWriter.Close();
            }
            catch (Exception e)
            {
                Exception baseException = e.GetBaseException();
                Console.WriteLine("WriteLowRepMersFile exception: " + e.Message + " " + e.HResult.ToString("X16") + ", base exception " + baseException.HResult.ToString("X16"));
                writeSucceeded = false;
                if (flushWriter != null)
                    flushWriter.Close();
            }

            return writeSucceeded;
        }

        private SingletonCollection GetSingletonTable(SingletonCollection singletonTable)
        {
            SingletonCollection newSingletonTable = null;

            lock (singletonPool)
            {
                // is there a pooled buffer of this size waiting to be re-used?
                for (int i = 0; i < singletonPool.Count; i++)
                {
                    if (singletonPool[i].Capacity == singletonTable.Capacity)
                    {
                        newSingletonTable = singletonPool[i];
                        singletonPool.RemoveAt(i);
                        //Console.WriteLine("allocating used singleton table (" + singletonTable.collectionNo + "). " + singletonLength + " k-mers");
                        break;
                    }
                }
            }

            if (newSingletonTable == null)
            {
                newSingletonTable = new SingletonCollection(singletonTable);
                //Console.WriteLine("allocating new singleton table (" + singletonTable.collectionNo + "). " + singletonLength + " k-mers");
            }
            else
            {
                newSingletonTable.Recycle(singletonTable.partitionNo);
                //Console.WriteLine("recycling singleton table (" + singletonTable.collectionNo + "). " + singletonLength + "/" + singletonTable.Capacity + " k-mers");
            }

            return newSingletonTable;
        }

        private void ReturnSingletonTable(SingletonCollection singletonTable)
        {
            lock (singletonPool)
            {
                singletonPool.Add(singletonTable);
            }
            //Console.WriteLine("returned singleton table (" + singletonTable.collectionNo + "). " + singletonTable.Capacity + " k-mers");
        }

        private int RoundUpSingletonSize(int requestedSize, int[] possibleSingletonSizes)
        {
            int sizeToUse = 0;

            foreach (int possibleSize in possibleSingletonSizes)
            {
                if (possibleSize > requestedSize)
                {
                    sizeToUse = possibleSize;
                    break;
                }
            }

            if (sizeToUse == 0)
                sizeToUse = possibleSingletonSizes[possibleSingletonSizes.Length - 1];

            return sizeToUse;
        }
    }

    public class TraceEntry
    {
        public int thread;
        public int place;
        public int partition;
        public int index;
        public ulong value;

        public TraceEntry(int thread, int place, int partition, int index, ulong value)
        {
            this.thread = thread;
            this.place = place;
            this.partition = partition;
            this.index = index;
            this.value = value;
        }
    }

    class FlushingMerThreadParams
    {
        public MerTables merTable;
        public int partitionNo;
    }

    public class LowRepMerBuffer
    {
        public int idx;
        public int limit;
        public bool bufferActive;
        public int bufferNo;
        public ulong[] keys;
        public long[] values;
    }

    public class FlushingSingletonsThreadParams
    {
        public EventWaitHandle flushingEWH;
        public Queue<SingletonCollection> buffersToBeFlushed;
    }
}