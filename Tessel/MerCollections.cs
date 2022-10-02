using System;
using System.Collections.Generic;
//using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.Threading;
using System.Diagnostics;

namespace MerCollections
{
    // A specialised Dictionary for holding k-mers for Tessel
    // These dictionaries never have entries deleted and are never automatically resized. These restrictions allow them to be lock-free.
    // the shared k-mer dictionaries are never resized, but the per-thread overflow dictionaries can be manually resized (safely)

    public class MerDictionary
    {
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        //public Entry[] entries;                             // hash table contents. An array of entries linked into lists
        public int[] entries;                               // next pointer for entries
        public ulong[] keys;                                // keys from entries in parallel array to support in-situ sorting (reduced memory and parallelism)
        public long[] values;                               // and the corresponding values
        private int lengthBuckets = 0;                      // current length of buckets array
        public int lengthEntries = 0;                       // and the length of the entries array
        private int count;                                  // no. of allocated (in-use) entries
        private int fullTrigger;                            // table is full when this many entries have been added 
        public bool sorted;                                 // set when table flattened into ordered sets of k-mers+counts at the end of the counting phase

        //public struct Entry : IComparable<Entry>
        //{
        //    public int next;                                // next list pointer for this entry list
        //    public ulong key;                               // the mer being stored - moved into parallel array to support in-place sorting
        //    public long value;                              // the value (ulong count pair) associated with the key. Should never hit sign bit. Interlocked.Add will only accept signed parameters
        //    // no need for hashcode as comparing ulongs is fast
        //    public int CompareTo(Entry y)                   // CompareTo used in Sort
        //    {
        //        if (this.key < y.key)
        //            return -1;
        //        if (this.key == y.key)
        //            return 0;
        //        return +1;
        //    }
        //}

        // constructor. Must always pass in the capacity. No () constructor
        public MerDictionary(int capacity, int merSize, int noThreads)
        {
            if (capacity > 0)
                Initialize(capacity, merSize, noThreads);
            else
                throw (new ApplicationException("mer dictionary size <= 0"));
        }

        private void Initialize(int capacity, int merSize, int noThreads)
        {
            this.buckets = new int[capacity + capacity / 3];
            for (int i = 0; i < this.buckets.Length; i++)
            {
                this.buckets[i] = -1;
            }
            //this.entries = new Entry[capacity];
            this.entries = new int[capacity];
            this.keys = new ulong[capacity];
            this.values = new long[capacity];
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            this.fullTrigger = capacity - noThreads;                        // trigger a resize of 'full' condition when all the entries are in-use - allowing for concurrent filling threads
            if (merSize == 32)
                HashHelpers.rightShiftBases = 1;
            else
                HashHelpers.rightShiftBases = 64 - merSize * 2 + 1;         // how much to shift to right-adjust a kmer 
        }

        // A lock-free (unchecked) add to a mer dictionary (which is guaranteed not to have deletions). The 'lock free' nature of this code will rarely allow
        // the same entry to be present multiple times in the table. Post-race increments will always update the first entry added, and the
        // merge phase will ensure correct counts.
        // 
        // Return value says whether this Add pushed the table past its 'full' point.
        //
        public bool Add(ulong key, long value)
        {
            //if (key == 0x1bb9df5864f8a0bd)
            //    Debugger.Break();
            int hashCode = HashHelpers.HashMer(key);
            int targetBucket = hashCode % this.buckets.Length;

            //if (targetBucket == 0)
            //    Debugger.Break();

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = Interlocked.Increment(ref this.count) - 1;

            //this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            //this.entries[index].key = key;
            //this.entries[index].value = value;
            this.entries[index] = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.keys[index] = key;
            this.values[index] = value;   
            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            return this.count >= this.fullTrigger;
        }

        public bool UpdateIfPresent(ulong key, long value)
        {
            int hashCode = HashHelpers.HashMer(key);

            //for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i])
            {
                //if (this.entries[i].key == key)
                if (this.keys[i] == key)    
                {
                    //Interlocked.Add(ref this.entries[i].value, value);
                    Interlocked.Add(ref this.values[i], value);   
                    return true;
                }
            }

            return false;
        }

        public int FindEntry(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key);
            //for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i])
            {
                //if (this.entries[i].key == key)
                if (this.keys[i] == key)    
                {
                    return i;
                }
            }

            return -1;
        }

        // can be called manually but is not thread-safe. Never called for shared repeatedMers but is called for per-thread overflow tables
        public void Resize()
        {
            int lengthNewBuckets = this.lengthBuckets * 5 / 4;
            int[] newBuckets = new int[lengthNewBuckets];
            for (int i = 0; i < newBuckets.Length; i++)
            {
                newBuckets[i] = -1;
            }
            int newEntriesLength = this.lengthEntries * 5 / 4;
            //Entry[] newEntries = new Entry[newEntriesLength];
            int[] newEntries = new int[newEntriesLength];
            ulong[] newKeys = new ulong[newEntriesLength];
            long[] newValues = new long[newEntriesLength];
            //Array.Copy(this.entries, 0, newEntries, 0, this.count);
            Array.Copy(this.keys, 0, newKeys, 0, this.count);
            Array.Copy(this.values, 0, newValues, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                //int bucket = HashHelpers.HashMer(newEntries[j].key) % lengthNewBuckets;
                int bucket = HashHelpers.HashMer(newKeys[j]) % lengthNewBuckets;
                //newEntries[j].next = newBuckets[bucket];
                newEntries[j] = newBuckets[bucket];
                newBuckets[bucket] = j;
            }
            this.buckets = newBuckets;
            this.entries = newEntries;
            this.keys = newKeys;
            this.values = newValues;
            this.lengthBuckets = newBuckets.Length;
            this.lengthEntries = newEntries.Length;
            this.fullTrigger = this.lengthEntries;
        }

        public void Clear()
        {
            if (count > 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, count);
                Array.Clear(keys, 0, count);
                Array.Clear(values, 0, count);
                count = 0;
            }
        }

        // An in-place reduction of a MerDictionary. All entries with a depth < minReps are copied to the culled arrays 
        // and entries to be kept are moved towards the beginning of the table. Assumes that the entries are inserted 
        // sequentially and that there are no empty slots. Culled entries are placed in the shared 'culled' arrays. If these arrays
        // are filled, they are moved to the 'filled' lists.
        //
        public int Reduce(long minReps, MerTables merTable)
        {
            int culledCount = 0;
            LowRepMerBuffer culledBuffer = merTable.culledBuffer;

            // clear out all the buckets
            for (int i = 0; i < buckets.Length; i++)
                buckets[i] = -1; 

            // and either re-add or cull each entry
            int originalCount = count;
            count = 0;
            for (int i = 0; i < originalCount; i++)
            {
                //ulong key = this.entries[i].key;
                //long value = this.entries[i].value;
                ulong key = this.keys[i];
                long value = this.values[i];      
                long summedValue = (value >> 32) + (value & 0xffffffff);

                if (summedValue > minReps)
                    // keeping this entry as it has enough counts
                    Add(key, value);                            // must always be room in the table - worst case is rewriting all entries in-place
                else
                {
                    // flush this entry to the (shared) culled buffer
                    int ci = Interlocked.Increment(ref culledBuffer.idx) - 1;

                    if (ci >= culledBuffer.limit)
                    {
                        lock (merTable.culledLock)
                        {
                            // re-test now that we have the lock (some other thread could have flushed the array while we were waiting for the lock)
                            if (culledBuffer.bufferActive)
                            {
                                // this buffer is being saved so it's no longer active
                                culledBuffer.bufferActive = false;
                                // add the full buffer to the list of filled buffers
                                if (merTable.filledCulledBuffers == null)
                                    merTable.filledCulledBuffers = new List<LowRepMerBuffer>();
                                merTable.filledCulledBuffers.Add(merTable.culledBuffer);

                                // prepare a new culledBuffer
                                LowRepMerBuffer newCulledBuffer = new LowRepMerBuffer();
                                newCulledBuffer.keys = new ulong[merTable.culledBuffer.keys.Length];
                                newCulledBuffer.values = new long[merTable.culledBuffer.values.Length];
                                newCulledBuffer.idx = 0;
                                newCulledBuffer.limit = culledBuffer.limit;
                                newCulledBuffer.bufferActive = true;
                                newCulledBuffer.bufferNo = culledBuffer.bufferNo + 1;

                                // and (atomically) make it available to any concurrent threads
                                merTable.culledBuffer = newCulledBuffer;
                            }
                        }

                        // remember to use the new buffer for this thread    
                        culledBuffer = merTable.culledBuffer;   
                        // and get a new index for this new buffer
                        ci = Interlocked.Increment(ref culledBuffer.idx) - 1;
                    }

                    culledBuffer.keys[ci] = key;
                    culledBuffer.values[ci] = value;
                    culledCount++;
                }
            } // for each entry in the original entries table

            return culledCount;
        }

        public int Sort()
        {
            // just return if we've already sorted this table
            if (sorted)
                return count;

            // return the bucket data structure for GC as these hash tables are never re-used
            buckets = null;

            // sort the entries in-place 
            //Array.Sort<Entry>(entries, 0, count);
            //Array.Sort<ulong, long>(keys, values, 0, count);
            MerArraySorter<long> merTableSorter = new MerArraySorter<long>(keys, values);
            merTableSorter.Sort(0, count - 1);

            //for (int i = 0; i < count - 2; i++)
            //    if (keys[i] > keys[i + 1])
            //    {
            //        Console.WriteLine("sort failed");
            //        break;
            //    }

            return count;
        }

        public IEnumerator<KeyValuePair<ulong, long>> GetEnumerator()
        {
            for (int i = 0; i < count; i++)
            {
                //KeyValuePair<ulong, long> kvp = new KeyValuePair<ulong, long>(entries[i].key, entries[i].value);
                KeyValuePair<ulong, long> kvp = new KeyValuePair<ulong, long>(keys[i], values[i]);       
                yield return kvp;
            }
        }

        public int Count
        {
            get
            {
                return (this.count);
            }

            set
            {
                count = value;
            }
        }

        public int Capacity
        {
            get
            {
                if (this.entries == null)
                    return this.count;
                else
                    return this.entries.Length;
            }
        }

        public long this[ulong key]
        {
            get
            {
                int index = this.FindEntry(key);
                if (index >= 0)
                {
                    //return this.entries[index].value;
                    return this.values[index];    
                }
                throw (new ApplicationException("not found"));
            }
        }

        public void CheckHashTable(out Dictionary<int, int> chains, out List<ulong> longestChain)
        {
            chains = new Dictionary<int, int>();
            longestChain = new List<ulong>();
            int maxChainLength = 0;
            int maxChainIdx = 0;

            int inUseBuckets = 0;
            for (int bi = 0; bi < lengthBuckets; bi++)
                if (buckets[bi] >= 0)
                    inUseBuckets++;
            chains.Add(-1, inUseBuckets);
            chains.Add(-2, lengthBuckets);
            chains.Add(-3, lengthEntries);

            for (int bi = 0; bi < lengthBuckets; bi++)
            {
                int bucket = buckets[bi];
                if (bucket == -1)
                    continue;

                int chainLength = 1;

                //for (int i = this.entries[bucket].next; i >= 0; i = this.entries[i].next)
                for (int i = this.entries[bucket]; i >= 0; i = this.entries[i])
                    chainLength++;

                if (chains.ContainsKey(chainLength))
                    chains[chainLength]++;
                else
                    chains.Add(chainLength, 1);

                if (chainLength > maxChainLength)
                {
                    maxChainLength = chainLength;
                    maxChainIdx = bi;
                    longestChain.Clear();
                    for (int i = this.entries[bucket]; i >= 0; i = this.entries[i])
                        longestChain.Add(this.keys[i]);
                }
            }

            chains.Add(-4, maxChainIdx);
        }

        internal void QuickSort(int left, int right)
        {
            int i = left;
            int j = right;

            int middle = (i + j) >> 1;
            ulong pivot = keys[middle];

            while (i <= j)
            {
                while (keys[i] < pivot)
                    i++;
                while (keys[j] > pivot)
                    j--;

                if (i <= j)
                {
                    ulong mer = keys[i];
                    keys[i] = keys[j];
                    keys[j] = mer;

                    long value = values[i];
                    values[i] = values[j];
                    values[j] = value;

                    i++;
                    j--;
                }
            }

            if (left < j)
                QuickSort(left, j);

            if (i < right)
                QuickSort(i, right);
        }
    }

    // A hash-based collection of unique mers. Effectively a Dictionary without the values, or a HashSet without the 'set' operations.
    // Adding entries to this collection is thread-safe as there are no deletions. Used to hold singletons only.
    // These structures are expected to be recycled. Entries are added to this collection until it is full, then it is flushed and emptied for re-use.
    // There are never any deletions, and so no need for a free list.
    // 
    public class SingletonCollection
    {
        private static int nextCollectionNo = 0;            // ordinal for allocated collections - used to track recycling only
        public int collectionNo;                            // unique no for this collection instance
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        //public Entry[] entries;                             // hash table contents. An array of entries linked into lists
        public int[] entries;                               // hash table contents. An array of entries linked into lists
        public long[] keys;                                 // actual keys from entries - taken out of Entry to allow for easy in-place sorting (and reduced memory usage, and parallelism)
        public bool sorted;                                 // has this collection been sorted already?
        public bool condensed;                              // has this collection been reduced already?
        public int length = 0;                              // current length of buckets and entries arrays
        private int count = 0;                              // count of allocated entries - and index to next unused entry slot
        //public IEqualityComparer<ulong> equalityComparer;   // mer comparer (equality and hash - masking off the status bits)
        //private IComparer<long> comparer;                   // entry comparer (longs compared as ulongs)
        private int merSize;                                // remember for re-initialisation
        private ulong statusBitsMask;                       // used to remove status bits before searches
        private ulong activeBitMask;                        // used to determine whether an entry is active or not before condensing/writing
        private ulong rcBitMask;                            // used to remember if the kMer is as-read or in RC form
        private int upperBitsLength = 0;                    // how many 'upper' bits there are (how far to >> the in-table kMer to make room)
        private ulong upperBits = 0;                        // upper bits of every kMer in this partition (shifted partitionNo)
        public byte[] byteBuffer;                           // byte buffer used when writing out to file  
        public bool[] rcFlags;                              // captures the RC bits when the active entries in an in-memory singleton table are expanded in preparation for merging

        public bool flushed = false;                        // this collection has already been flushed
        //public int merSize = 0;                             // length of the kMers in this collection
        public int partitionNo = 0;                         // which partition of the singleton table does this collection belong to (active and queued for consolidation/writing)
        public int flushNo = 0;                             // ... <XXXX>_partitionNo_flushNo.bst
        public DateTime timeWhenQueued;                     // ... when was this collection queued for flushing?
        public bool lastBuffer = false;                     // signal flushing thread that this is the end and it should terminate now

        public const int byteBufferLen = 8192*9;            // size of byte buffer used when flushing tables to storage (has to be a multiple of 8+1 bytes)

        //public struct Entry : IComparable<Entry>
        //{
        //    public int next;                                // next list pointer for this entry list
        //    public long key;                                // the mer being stored. Mers are always partitioned over multiple tables, based on their top bases, so the RHS of each mer will be empty and can be used for status flags. These are masked off for comparisons.
        //                                                    // this is really a ulong but the Interlocked primitive only accepts longs... so much casting needed

        //    public int CompareTo(Entry y)                 // needed only with Array.Sort - built-in to hand-crafted sort
        //    {
        //        if ((ulong)this.key < (ulong)y.key)
        //            return -1;
        //        if (this.key == y.key)
        //            return 0;
        //        return +1;
        //    }
        //}

        // Constructor. Must always pass in the capacity. No () constructor
        public SingletonCollection(int capacity, int merSize, ulong statusBitsMask, ulong activeBitMask, ulong rcBitMask, int partitionNo, int upperBitsLength)
        {
            if (capacity > 0)
                Initialize(capacity, merSize, statusBitsMask, activeBitMask, rcBitMask, partitionNo, upperBitsLength);
            else
                throw (new ApplicationException("mer collection size <= 0"));
        }

        public SingletonCollection(SingletonCollection previousSC)
        {
            Initialize(previousSC.Capacity, previousSC.merSize, previousSC.statusBitsMask, previousSC.activeBitMask, previousSC.rcBitMask, previousSC.partitionNo, previousSC.upperBitsLength);
        }

        private void Initialize(int capacity, int merSize, ulong statusBitsMask, ulong activeBitMask, ulong rcBitMask, int partitionNo, int upperBitsLength)
        {
            this.buckets = new int[capacity + capacity / 3];
            for (int i = 0; i < this.buckets.Length; i++)
                this.buckets[i] = -1;
            //this.entries = new Entry[capacity];
            this.entries = new int[capacity];
            this.keys = new long[entries.Length];
            this.length = entries.Length;
            this.merSize = merSize;
            this.count = 0;
            this.sorted = false;
            this.flushed = false;
            //this.equalityComparer = new MerEqualityComparer(mask);
            this.collectionNo = nextCollectionNo;
            nextCollectionNo++;
            //this.comparer = new SingletonComparer();
            this.statusBitsMask = statusBitsMask;
            this.activeBitMask = activeBitMask;
            this.rcBitMask = rcBitMask;
            this.partitionNo = partitionNo;
            this.upperBits = (ulong)partitionNo << (64 - upperBitsLength);
            this.upperBitsLength = upperBitsLength;
            this.byteBuffer = new byte[byteBufferLen];
            HashHelpers.rightShiftBases = 64 - merSize * 2 + 1;
            //Console.WriteLine("initialised singletons # " + this.collectionNo + " length= " + capacity);
        }

        public void Recycle(int partitionNo)
        {
            if (this.count >= 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, this.count);
                Array.Clear(keys, 0, this.count);
                this.count = 0;
                this.sorted = false;
                this.condensed = false;
                this.flushed = false;
                this.partitionNo = partitionNo;
                this.upperBits = (ulong)partitionNo << (64 - this.upperBitsLength);
            }
        }

        // either inserts a key into the hash table (and returns true) or, if the key is already present, returns the existing key and false.
        // the mask is used to delete the status bits at the end of the mer key.
        public bool TryInsertKey(ulong key, out int currentIdx)
        {
            ulong maskedKey = key & statusBitsMask;
            int hashCode = HashHelpers.HashMer(maskedKey);
            int targetBucket = hashCode % this.buckets.Length;

            currentIdx = -1;

            //for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i])
            {
                //if (((ulong)this.entries[i].key & mask) == (key & mask))
                if (((ulong)this.keys[i] & statusBitsMask) == maskedKey)
                {
                    currentIdx = i;
                    return false;                               // insert attempt failed - return index of entry
                }
            }

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = Interlocked.Increment(ref this.count) - 1;


            //this.entries[index].next = this.buckets[targetBucket];     // point this new entry's next to the last entry pointed by buckets[targetBucket]
            //this.entries[index].key = (long)key;
            this.entries[index] = this.buckets[targetBucket];     // point this new entry's next to the last entry pointed by buckets[targetBucket]     
            this.keys[index] = (long)key;
            this.buckets[targetBucket] = index;                   // if collision occurs, update value in buckets[index] to point to new slot in entries[]
            currentIdx = index;
            return true;                                          // return true if Insert succeeded
        }

        // adds an entry to the table, regardless of whether or not it is already there
        // used to copy entries form a sparse singleton table back to the curret one being filled
        // duplicates will be taken care of during merge
        public void ForceInsertKey(ulong key)
        {
            ulong maskedKey = key & statusBitsMask;
            int hashCode = HashHelpers.HashMer(maskedKey);
            int targetBucket = hashCode % this.buckets.Length;
            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = Interlocked.Increment(ref this.count) - 1;

            //this.entries[index].next = this.buckets[targetBucket];     // point this new entry's next to the last entry pointed by buckets[targetBucket]
            //this.entries[index].key = (long)key;
            this.entries[index] = this.buckets[targetBucket];     // point this new entry's next to the last entry pointed by buckets[targetBucket]     
            this.keys[index] = (long)key;
            this.buckets[targetBucket] = index;                   // if collision occurs, update value in buckets[index] to point to new slot in entries[]
        }

        // close a collection and sort the keys (sorting the un-expanded kMers)
        public void Sort()
        {
            // only do this once per allocated collection
            if (sorted)
                return;
            sorted = true;

            //Array.Sort<Entry>(entries, 0, this.count);
            //Array.Sort<long>(keys, 0, this.count, comparer);

            SingletonArraySorter singletonSorter = new SingletonArraySorter(keys);
            singletonSorter.Sort(0, this.count - 1);
            sorted = true;

            //for (int i = 0; i < this.count - 2; i++)
            //    if ((ulong)keys[i] > (ulong)keys[i + 1])
            //    {
            //        Console.WriteLine("sort failed @" + i + " " + keys[i].ToString("X16") + " " + keys[i + 1].ToString("X16"));
            //        break;
            //    }
        }

        // compress the singletons array in-place by skipping over inactive entries (shift active entries towards the start of the array)
        public int Condense()
        {
            if (condensed)
                return count;

            int merIdx = 0;
            for (int s = 0; s < this.count; s++)
            {
                //MerCollection.Entry mer = sortedSingletons[s];
                long mer = this.keys[s];

                //if (((ulong)mer.key & singletons.activeBitMask) != 0)
                if (((ulong)mer & this.activeBitMask) != 0)
                {
                    this.keys[merIdx] = mer;
                    merIdx++;
                }
            }

            count = merIdx;
            condensed = true;

            return count;
        }

        // extend the singleton kMers in a partition to their full length, and save the rcFlags for use during the merge
        public void ExtendAndSaveRC()
        {
            rcFlags = new bool[count];

            for (int i = 0; i < count; i++)
            {
                ulong mer = (ulong)this.keys[i];

                rcFlags[i] = (mer & rcBitMask) != 0;
                this.keys[i] = (long)((mer >> this.upperBitsLength) | this.upperBits);
                //if ((ulong)this.keys[i] == 0xa1c5aafdc356f670)
                //    Debugger.Break();
            }
        }

        public ulong Extend(ulong mer, out bool rcMer)
        {
            rcMer = (mer & rcBitMask) != 0;

            return (mer >> this.upperBitsLength) | this.upperBits;
        }

        public IEnumerator<ulong> GetEnumerator()
        {
            for (int i = 0; i < this.count; i++)
            {
                //yield return (ulong)entries[i].key;
                yield return (ulong)keys[i];   
            }
        }

        public int Count
        {
            get
            {
                return this.count;
            }

            set
            {
                this.count = value;
            }
        }

        public int Capacity
        {
            get
            {
                return (this.length);
            }
        }

    } // MerCollection

    //public class SingletonComparer : IComparer<long>
    //{
    //    public int Compare(long x, long y)
    //    {
    //        if ((ulong)x < (ulong)y)
    //            return -1;
    //        if (x == y)
    //            return 0;
    //        return 1;
    //    }
    //}

    public class MerEqualityComparer : IEqualityComparer<ulong>
    {
        // used to mask off the flag bits in the singleton tables
        ulong mask;

        public MerEqualityComparer(ulong mask)
        {
            this.mask = mask;
        }

        public bool Equals(ulong x, ulong y)
        {
            return (x & mask) == (y & mask);
        }

        public int GetHashCode(ulong key)
        {
            return HashHelpers.HashMer(key); 
        }
    }

    internal static class HashHelpers
    {
        // this class used to contain code to find prime numbers for the sizing of bucket arrays. Prime table sizes are only a defence
        // against poor hashing functions (and unwanted periodicity) and are not really needed.

        internal static int rightShiftBases = 0;

        internal static int HashMer(ulong key)
        {
            // experiments with optimised code showed that the simple folded hash performed almost as well for randomness as the Murmur hash, and was short enough
            // to be inlined. The overall result was that the folded hash was faster. Copying the body of this function in-line made no difference. 

            const ulong positiveInt32Mask = 0x7fffffff;

            //int foldedUL = (int)(key >> 32) ^ (int)(key & 0xffffffff);
            //return (int)(foldedUL) & positiveInt32Mask;
            //return (((int)key) ^ (int)(key >> 31)) & positiveInt32Mask;

            // as we're dealing with left-adjusted kmers, the hash function tries to maximise the number of bits going into the hash by taking the top 32 bits, and the low significant 32 bits

            //return (int)((key >> 32) ^ ((key >> rightShiftBases + 1) & 0xffffffff)) & positiveInt32Mask;
            return (int)(((key >> 32) ^ (key >> rightShiftBases)) & positiveInt32Mask);
        }
    } // HashHelpers

    public struct MerArraySorter<V>
    {
        private ulong[] keys;
        private V[] values;

        internal MerArraySorter(ulong[] keys, V[] values)
        {
            this.keys = keys;
            this.values = values;
        }

        // QuickSort
        internal void Sort(int left, int right)
        {
            if (right - left < 7)
            {
                InsertionSort(left, right);
                return;
            }

            int i = left;
            int j = right;

            int middle = i + ((j - i) >> 1);

            SwapIfGreater(i, middle); // swap the low with the mid point 
            SwapIfGreater(i, j);      // swap the low with the high
            SwapIfGreater(middle, j); // swap the middle with the high

            ulong pivot = keys[middle];

            while (i <= j)
            {
                while (keys[i] < pivot)
                    i++;
                while (keys[j] > pivot)
                    j--;

                if (i <= j)
                {
                    ulong mer = keys[i];
                    keys[i] = keys[j];
                    keys[j] = mer;
                    V value = values[i];
                    values[i] = values[j];
                    values[j] = value;
                    i++;
                    j--;
                }
            }

            if (left < j)
                Sort(left, j);

            if (i < right)
                Sort(i, right);
        }

        internal void SwapIfGreater(int a, int b)
        {
            if (keys[a] > keys[b])
            {
                ulong key = keys[a];
                V value = values[a];
                keys[a] = keys[b];
                values[a] = values[b];
                keys[b] = key;
                values[b] = value;
            } 
        }

        internal void InsertionSort(int start, int end)
        {
            for (int x = start + 1; x <= end; x++)
            {
                ulong key = keys[x];
                V value = values[x];
                int j = x - 1;

                while (j >= 0 && key < keys[j])
                {
                    keys[j + 1] = keys[j];
                    values[j + 1] = values[j];
                    j--;
                }

                keys[j + 1] = key;
                values[j + 1] = value;
            }
        }
    }

    public struct SingletonArraySorter
    {
        // save on parameter passing overhead
        // array is longs (because of Interlocked limitations) but really ulongs and have to be sorted as such.
        // the singletons also have status bits at the RHS but these don't matter as far as sorting goes - so no masking needed
        private long[] keys;

        internal SingletonArraySorter(long[] keys)
        {
            this.keys = keys;
        }

        // QuickSort
        internal void Sort(int left, int right)
        {
            if (right - left < 9)
            {
                InsertionSort(left, right);
                return;
            }

            int i = left;
            int j = right;

            int middle = i + ((j - i) >> 1);

            // avoids unbalanced recursion if starting array is sorted (or partially)
            SwapIfGreater(i, middle); // swap the low with the mid point 
            SwapIfGreater(i, j);      // swap the low with the high
            SwapIfGreater(middle, j); // swap the middle with the high

            ulong pivot = (ulong)keys[middle]; 

            while (i <= j)
            {
                while ((ulong)keys[i] < pivot)
                    i++;
                while ((ulong)keys[j] > pivot)
                    j--;

                if (i <= j)
                {
                    long mer = keys[i];
                    keys[i] = keys[j];
                    keys[j] = mer;

                    i++;
                    j--;
                }
            }

            if (left < j)
                Sort(left, j);

            if (i < right)
                Sort(i, right);
        }

        internal void SwapIfGreater(int a, int b)
        {
            if ((ulong)keys[a] > (ulong)keys[b])
            {
                long key = keys[a];
                keys[a] = keys[b];
                keys[b] = key;
            }
        }

        internal void InsertionSort(int start, int end)
        {
            for (int x = start + 1; x <= end; x++)
            {
                long key = keys[x];
                int j = x - 1;

                while (j >= 0 && (ulong)key < (ulong)keys[j])
                {
                    keys[j + 1] = keys[j];
                    j--;
                }

                keys[j + 1] = key;
            }
        }
    }

} // namespace