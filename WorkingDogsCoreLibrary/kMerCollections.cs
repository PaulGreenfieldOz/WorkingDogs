using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace WorkingDogsCore
{
    // ===========================================
    // ---------------- MerDictionary ------------
    // ===========================================
    //
    // A partitioned kMer dictionary class. 
    //
    // Creates an array of dictionaries, each of which should be smaller than the maximum safe size for a Dictionary (~50M)
    //
    public class MerDictionary<TV>
    {
        int keySize = 0;                                                                    // number of internal partitions is 4 ** keySize
        public int noOfPartitions = 0;
        int keyBaseShift = 0;
        bool partitioned = false;
        kMerDictionary<TV>[] dictionaryPartitions = null;
        kMerDictionary<TV> dictionary = null;
        long capacity = 0;

        const int maxTable = 25000000;

        public MerDictionary(long dictionarySize, int merSize)
        {
            if (dictionarySize > maxTable)
            {
                partitioned = true;
                keySize = (int)Math.Ceiling(Math.Log(dictionarySize / maxTable, 4));    // find minimum key size
                if (keySize < 1)
                    keySize = 1;                                                        // must partition on at least one base
                noOfPartitions = (int)Math.Pow(4, keySize);                             // giving this many partitions
                keyBaseShift = 64 - keySize * 2;                                        // shift right this many bits to extract the partition no.
                long partitionLength = (int)(dictionarySize / noOfPartitions);          // with this average length (but scaled to reflect canonical distributions)

                dictionaryPartitions = new kMerDictionary<TV>[noOfPartitions];
                for (int i = 0; i < noOfPartitions; i++)
                {
                    int scaledPartitionLength = (int)(2 * partitionLength * (noOfPartitions - i) / noOfPartitions); 
                    dictionaryPartitions[i] = new kMerDictionary<TV>(scaledPartitionLength, merSize);
                }
            }
            else
            {
                partitioned = false;
                dictionary = new kMerDictionary<TV>((int)dictionarySize, merSize);
            }

            capacity = dictionarySize;
        }

        public void Add(ulong key, TV value)
        {          
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                dictionaryPartitions[partition].Add(key, value);
            }
            else
            {
                dictionary.Add(key, value);
            }
        }

        public void AddIfNotPresent(ulong key, TV value)
        {
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                dictionaryPartitions[partition].AddIfNotPresent(key, value);
            }
            else
            {
                dictionary.AddIfNotPresent(key, value);
            }
        }

        public bool ContainsKey(ulong key)
        {
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                return dictionaryPartitions[partition].ContainsKey(key);
            }
            else
            {
                return dictionary.ContainsKey(key);
            }
        }

        public int Count
        {
            get
            {
                if (partitioned)
                {
                    int totalCount = 0;
                    for (int i = 0; i < noOfPartitions; i++)
                        totalCount += dictionaryPartitions[i].Count;
                    return totalCount;
                }
                else
                {
                    return dictionary.Count;
                }
            }
        }

        public long Capacity
        {
            get
            {
                return capacity;
            }
        }

        public int[] PartitionCounts
        {
            get
            {
                if (partitioned)
                {
                    int[] partitionCounts = new int[noOfPartitions];
                    for (int i = 0; i < noOfPartitions; i++)
                        partitionCounts[i] = dictionaryPartitions[i].Count;
                    return partitionCounts;
                }
                else
                {
                    int[] partitionCounts = new int[1];
                    partitionCounts[0] = dictionary.Count;
                    return partitionCounts;
                }
            }
        }

        public void Clear()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                    dictionaryPartitions[i].Clear();
            }
            else
            {
                dictionary.Clear();
            }
        }

        public void Optimise()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                    dictionaryPartitions[i].Optimise();
            }
            else
            {
                dictionary.Optimise();
            }
        }

        public bool TryGetValue(ulong key, out TV value)
        {
            if (partitioned)
            {
                int partition = (int)(key >> (64 - (keySize * 2)));
                return dictionaryPartitions[partition].TryGetValue(key, out value);
            }
            else
            {
                return dictionary.TryGetValue(key, out value);
            }
        }

        public TV this[ulong key]
        {
            get
            {
                if (partitioned)
                {
                    int partition = (int)(key >> (64 - (keySize * 2)));
                    return dictionaryPartitions[partition][key];
                }
                else
                {
                    return dictionary[key];
                }
            }
            set
            {
                if (partitioned)
                {
                    int partition = (int)(key >> (64 - (keySize * 2)));
                    dictionaryPartitions[partition][key] = value;
                }
                else
                {
                    dictionary[key] = value;
                }
            }
        }

        public IEnumerator<KeyValuePair<ulong, TV>> GetEnumerator()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                {
                    foreach (KeyValuePair<ulong, TV> kvp in dictionaryPartitions[i])
                    {
                        yield return kvp;
                    }
                }
            }
            else
            {
                foreach (KeyValuePair<ulong, TV> kvp in dictionary)
                {
                    yield return kvp;
                }
            }
        }

    }

    // ====================================
    // ----------- MerHashSet -------------
    // ====================================
    //
    // A (possibly) partitioned kMer hashset. The HashSet equivalent of MerDictionary)
    // 
    public class MerHashSet
    {
        int noOfPartitions = 0;
        int keyBaseShift = 0;
        const int maxTable = 25000000;                  // safe size for a single hashset (can be doubled safely in scaling)
        bool partitioned = false;                       // partitioned or not
        kMerHashSet[] hashSetPartitions = null;         // partitioned HashSets
        kMerHashSet hashSet = null;                     // or a single non-partitioned HashSet

        // new contructor that allows equi-sized partitions (better for non-canonical kMers)
        public MerHashSet(long tableSize, int merSize, bool canonical)
        {
            InitialiseMerHashSet(tableSize, merSize, canonical);
        }

        // original constructor - assumed kMers were always canonical and so the partitions should always be scaled - kept for backwards compatibility
        public MerHashSet(long tableSize, int merSize)
        {
            InitialiseMerHashSet(tableSize, merSize, true);
        }

        private void InitialiseMerHashSet(long tableSize, int merSize, bool canonical)
        {
            if (tableSize > maxTable)
            {
                partitioned = true;

                int keySizeBases = (int)Math.Ceiling(Math.Log(tableSize / maxTable, 4)); // find minimum key size
                if (keySizeBases < 1)
                    keySizeBases = 1;                                                       // must partition on at least one base
                keyBaseShift = 64 - keySizeBases * 2;
                noOfPartitions = (int)Math.Pow(4, keySizeBases);                            // giving this many partitions
                long partitionLength = (int)(tableSize / noOfPartitions);       // and each partition starts off being this big and may then be resized

                hashSetPartitions = new kMerHashSet[noOfPartitions];
                for (int i = 0; i < noOfPartitions; i++)
                {
                    int scaledPartitionLength = canonical ? (int)(2 * partitionLength * (noOfPartitions - i) / noOfPartitions) : (int)(partitionLength+partitionLength/10);
                    hashSetPartitions[i] = new kMerHashSet(scaledPartitionLength, merSize);  // scale partition sizes to allow for effects of canonical kMer distributions
                }
            }
            else
            {
                partitioned = false;
                hashSet = new kMerHashSet((int)tableSize, merSize);
            }
        }

        public void Add(ulong key)
        {
            if (partitioned)
            { 
                int partition = (int)(key >> keyBaseShift);
                hashSetPartitions[partition].Add(key);
            }
            else
            {
                hashSet.Add(key);
            }
        }

        public void AddNoCheck(ulong key)
        {
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                hashSetPartitions[partition].AddNoCheck(key);
            }
            else
            {
                hashSet.AddNoCheck(key);
            }
        }

        public void AddIfNotPresent(ulong key)
        {
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                hashSetPartitions[partition].AddIfNotPresent(key);
            }
            else
            {
                hashSet.AddIfNotPresent(key);
            }
        }

        public bool Contains(ulong key)
        {
            if (partitioned)
            {
                int partition = (int)(key >> keyBaseShift);
                return hashSetPartitions[partition].Contains(key);
            }
            else
            {
                return hashSet.Contains(key);
            }
        }

        public long Count
        {
            get
            {
                if (partitioned)
                {
                    long totalCount = 0;
                    for (int i = 0; i < noOfPartitions; i++)
                        totalCount += hashSetPartitions[i].Count;
                    return totalCount;
                }
                else
                {
                    return hashSet.Count;
                }
            }
        }

        public void Optimise()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                {
                    if ((float)hashSetPartitions[i].Count / (float)hashSetPartitions[i].Capacity < 0.8)
                    {
                        hashSetPartitions[i].Optimise();
                        GC.Collect(2);
                    }
                }
            }
            else
            {
                hashSet.Optimise();
                GC.Collect(2);
            }
        }

        public void Clear()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                    hashSetPartitions[i].Clear();
            }
            else
            {
                hashSet.Clear();
            }
        }

        public void CopyEntries(ulong[] entriesArray)
        {
            int totalEntries = (int)Count;
            if (entriesArray.Length != totalEntries)
                Array.Resize<ulong>(ref entriesArray, totalEntries);

            if (partitioned)
            {
                int nextStartingPoint = 0;
                for (int i = 0; i < noOfPartitions; i++)
                {
                    hashSetPartitions[i].CopyEntries(nextStartingPoint, entriesArray);
                    nextStartingPoint += hashSetPartitions[i].Count;
                }
            }
            else
            {
                hashSet.CopyEntries(0, entriesArray);
            }
        }

        public void AddSubVariants()
        {
            List<ulong> merVariants = new List<ulong>(200);

            if (partitioned)
            {
                for (int p = 0; p < noOfPartitions; p++)
                {
                    hashSetPartitions[p].AddSubVariants();
                }
            }
            else
            {
                hashSet.AddSubVariants();
            }
        }

        public IEnumerator<ulong> GetEnumerator()
        {
            if (partitioned)
            {
                for (int i = 0; i < noOfPartitions; i++)
                {
                    foreach (ulong k in hashSetPartitions[i])
                    {
                        yield return k;
                    }
                }
            }
            else
            {
                foreach (ulong k in hashSet)
                {
                    yield return k;
                }
            }
        }

        public Dictionary<int, int>[] CheckHashTable()
        {
            if (partitioned)
            {
                Dictionary<int, int>[] chainArray = new Dictionary<int, int>[noOfPartitions];
                for (int i = 0; i < noOfPartitions; i++)
                    chainArray[i] = hashSetPartitions[i].CheckHashTable();
                return chainArray;
            }
            else
            {
                Dictionary<int, int>[] chainArray = new Dictionary<int, int>[1];
                chainArray[0] = hashSet.CheckHashTable();
                return chainArray;
            }
        }
    }

    // ===========================
    // --------- MerTable -------- 
    // ===========================
    //
    // A potentially very large MerDictionary (normally loaded from a .cbt file)
    // 
    // Implemented using a single Dictionary for small tables (faster) or an array of length-optimised Dictionaries for larger tables (slightly slower to load).
    // These tables are assumed to be read-only once loaded (indicated by LoadFinished) being called. The mers to be loaded are assumed to be ordered. 
    //
    // Intended for holding kMers loaded from a .cbt file and the kmers always have to be loaded in order for partitioned tables.
    // Partitioned data is first loaded into arrays and these are then asynchronously loaded into the KmerDictionary partitions.
    // 
    public class MerTable<TV>
    {
        bool partitioned = false;           // single Dictionary or multiple
        int merSize = 0;                    // how big are the kmers in the table
        const int maxTable = 20000000;      // safe size for a single Dictionary (can be doubled safely in scaling)
        bool dataNeedsSorting = false;      // set if out-of-order load detected
        ulong previousMer = 0;              // loads should be monotonically increasing

        kMerDictionary<TV> mers = null;     // small sets of mers in a single dictionary (a bit faster to load)
        kMerDictionary<TV>[] pmers = null;  // large sets of mers in length-optimised Dictionaries
        ulong[][] orderedMers = null;       // arrays used (and re-used) during partition loading
        TV[][] orderedValues = null;        // values are loaded into these arrays prior to the allocation of the corresponding dictionary partition
        ulong[] partitionBoundaries;        // upper possible k-mer in each partition
        int currentPartition = 0;           // current partition being loaded
        int cpi = 0;                        // index into partition buffer currently being filled (via Add calls)
        int currentBuffer = 0;
        IAsyncResult iarCopyToTable;
        CopyToTableDelegate ctd;

        int noParts = 0;                    // how many partitions
        int keyBaseShift = 0;               // # of bits to shift to get partition-key bases (bits)

        public MerTable(long tableSize, int merSize)
        {
            this.merSize = merSize;
            if (tableSize > maxTable)
            {
                partitioned = true;

                int keySizeBases = (int)Math.Ceiling(Math.Log(tableSize / maxTable, 4));    // find minimum key size
                if (keySizeBases < 1)
                    keySizeBases = 1;                                                       // must partition on at least one base
                keyBaseShift = 64 - keySizeBases * 2;
                noParts = (int)Math.Pow(4, keySizeBases);                                   // giving this many partitions
                int partitionLength = (int)(tableSize / noParts);                           // 
                int scaledPartitionLength = 2 * partitionLength;                            // rough guess as to space needed to hold the first (largest) partition of loaded kmers 

                orderedMers = new ulong[2][];
                orderedValues = new TV[2][];

                for (int b = 0; b < 2; b++)
                {
                    orderedMers[b] = new ulong[scaledPartitionLength];
                    orderedValues[b] = new TV[scaledPartitionLength];
                }
                ctd = new CopyToTableDelegate(CopyToTable);

                pmers = new kMerDictionary<TV>[noParts];

                ulong fillBases = 0xffffffffffffffff >> (keySizeBases * 2);
                partitionBoundaries = new ulong[noParts];
                for (ulong k = 0; k < (ulong)noParts; k++)
                    partitionBoundaries[k] = k << keyBaseShift | fillBases;
            }
            else
            {
                mers = new kMerDictionary<TV>((int)tableSize, merSize);
            }
        }

        public void Add(ulong key, TV value)
        {
            if (partitioned)
            // load into buffers first, then copy these into hash tables when each partition is complete (and its size is known)
            {
                // is this the first kmer to go into the next partition?
                if (key > partitionBoundaries[currentPartition])
                {
                    //Console.WriteLine(cpi + " mers added to partition " + currentPartition);
                    // wait for previous buffer copy to complete
                    if (iarCopyToTable != null && !iarCopyToTable.IsCompleted)
                    {
                        //Console.WriteLine("waiting for copy to finish");
                        iarCopyToTable.AsyncWaitHandle.WaitOne();
                    }
                    //Console.WriteLine("calling copy for partition " + currentPartition + " for buffer " + currentBuffer);
                    iarCopyToTable = ctd.BeginInvoke(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi, merSize, null, null);
                    //CopyToTable(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi);

                    if (currentBuffer == 0)
                        currentBuffer = 1;
                    else
                        currentBuffer = 0;

                    currentPartition++;
                    cpi = 0;
                }

                orderedMers[currentBuffer][cpi] = key;
                orderedValues[currentBuffer][cpi] = value;
                cpi++;

                if (cpi == orderedMers[currentBuffer].Length)
                {
                    int newLength = cpi + cpi / 10;
                    Array.Resize<ulong>(ref orderedMers[currentBuffer], newLength);
                    Array.Resize<TV>(ref orderedValues[currentBuffer], newLength);
                }
                if (key < previousMer)
                    dataNeedsSorting = true;
                previousMer = key;
            }
            else
            // not partitioned so just load straight into the hash table
            {
                mers.AddNoCheck(key, value);
            }
        }

        private delegate void CopyToTableDelegate(int partitionNo, ulong[] orderedMers, TV[] orderedValues, int merCount, int merSize);

        private void CopyToTable(int partitionNo, ulong[] orderedMers, TV[] orderedValues, int merCount, int merSize)
        {
            //Console.WriteLine("starting copy for partition " + partitionNo);
            kMerDictionary<TV> pmp = new kMerDictionary<TV>(merCount, merSize);
            pmers[partitionNo] = pmp;
            for (int i = 0; i < merCount; i++)
            {
                pmp.AddNoCheck(orderedMers[i], orderedValues[i]);
            }
            pmp.Optimise();
            //Console.WriteLine(merCount + " mers copied for partition " + partitionNo);
        }

        public bool LoadFinished()
        {
            if (dataNeedsSorting)
            {
                Console.WriteLine("kmers being loaded into table were not in sorted order");
                return false;
            }

            if (partitioned)
            {
                // wait for previous buffer copy to complete
                if (iarCopyToTable != null && !iarCopyToTable.IsCompleted)
                {
                    //Console.WriteLine("waiting for copy to finish before final copy");
                    iarCopyToTable.AsyncWaitHandle.WaitOne();
                }

                // copy final buffer into its hash table partition
                //Console.WriteLine("calling copy from LoadFinished for partition " + currentPartition + " for buffer " + currentBuffer);
                CopyToTable(currentPartition, orderedMers[currentBuffer], orderedValues[currentBuffer], cpi, merSize);

                for (int p = 0; p < pmers.Length; p++)
                    if (pmers[p] == null)
                        pmers[p] = new kMerDictionary<TV>(1, merSize);
            }

            // finished with these buffers now
            orderedMers = null;
            orderedValues = null;

            //Console.WriteLine(mersAdded + " mers added, " + mersCopied + " mers copied");
            return true;
        }

        public bool ContainsKey(ulong key)
        {
            if (partitioned)
            {
                int partIdx = (int)(key >> keyBaseShift);
                return pmers[partIdx].ContainsKey(key);
            }
            else
                return mers.ContainsKey(key);
        }

        public long Count
        {
            get
            {
                if (partitioned)
                {
                    long sum = 0;
                    for (int partIdx = 0; partIdx < noParts; partIdx++)
                        sum += pmers[partIdx].Count;
                    return sum;
                }
                else
                    return mers.Count;
            }
        }

        public void Clear()
        {
            if (partitioned)
            {
                for (int p = 0; p < noParts; p++)
                    pmers[p].Clear();
            }
            else
                mers.Clear();
        }

        public bool TryGetValue(ulong key, out TV value)
        {
            if (partitioned)
            {
                int partIdx = (int)(key >> keyBaseShift);
                return pmers[partIdx].TryGetValue(key, out value);
            }
            else
                return mers.TryGetValue(key, out value);
        }

        public TV this[ulong key]
        {
            get
            {
                if (partitioned)
                {
                    int partIdx = (int)(key >> keyBaseShift);
                    return pmers[partIdx][key];
                }
                else
                    return mers[key];
            }
            set
            {
                if (partitioned)
                {
                    int partIdx = (int)(key >> keyBaseShift);
                    pmers[partIdx][key] = value;
                }
                else
                    mers[key] = value;
            }
        }

        public IEnumerator<KeyValuePair<ulong, TV>> GetEnumerator()
        {
            if (partitioned)
            {
                for (int p = 0; p < noParts; p++)
                    foreach (KeyValuePair<ulong, TV> kvp in pmers[p])
                        yield return kvp;
            }
            else
            {
                foreach (KeyValuePair<ulong, TV> kvp in mers)
                    yield return kvp;
            }
        }
    }

    // ========================================
    // -------------- kMerDictionary ----------
    // ========================================

    // A specialised read-only Dictionary for holding k-mers (ulongs). meant for internal use only but declared public for testing
    // Tables are loaded during initialisation and from then on are read-only. Loading is single-threaded, and code is thread-safe afterwards.
    // 

    public class kMerDictionary<TV>
    {
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        private Entry[] entries;                            // hash table contents. An array of entries linked into lists
        private int lengthBuckets = 0;                      // current length of buckets array
        private int lengthEntries = 0;                      // and the length of the entries array
        private int count;                                  // no. of allocated (in-use) entries
        private int resizeTriggerPoint;                     // force a resize when the table gets this full 
        private int rightShiftBases;                        // number of bases to shift when right-justifying kMer in hashing

        public struct Entry
        {
            public int next;                                // next list pointer for this entry list
            public ulong key;                               // the mer being stored 
            public TV value;                                // the value (ulong count pair) associated with the key. 
        }

        // constructor. Must always pass in the capacity. No () constructor
        public kMerDictionary(int capacity, int merSize)
        {
            if (capacity > 0)
                Initialize(capacity, merSize);
            else
                throw (new ApplicationException("mer dictionary size <= 0"));
        }

        private void Initialize(int capacity, int merSize)
        {
            int bucketsLength = capacity + capacity / 3;
            if (bucketsLength % 2 == 0)
                bucketsLength++;
            this.buckets = new int[bucketsLength];

            for (int i = 0; i < buckets.Length; i++)
            {
                this.buckets[i] = -1;
            }
            this.entries = new Entry[capacity];
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            this.resizeTriggerPoint = capacity;
            this.rightShiftBases = merSize > 16 ? 64 - merSize * 2 : 0;     // right-adjust the low half of the kMer if there is one 
        }

        // A lock-free add to a mer dictionary (which is guaranteed not to have deletions or concurrent callers). 
        // This version checks for the presence of the key before trying to add it - and raises an exception if it's already present
        //
        public void Add(ulong key, TV value)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                    throw (new ApplicationException("mer already present"));
            }

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = count;
            count++;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;
            this.entries[index].value = value;
   
            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        public void AddIfNotPresent(ulong key, TV value)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                    return;
            }

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = count;
            count++;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;
            this.entries[index].value = value;

            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        // An unchecked lock-free Add to a mer dictionary (which is guaranteed not to have deletions or concurrent callers). 
        // The caller is expected to have checked for the presence if the key before calling AddNoCheck, or otherwise know that there will be no duplicates.
        // The MerTable code calls this version of Add because it is loading distinct kmers from a .cbt file.
        //
        public void AddNoCheck(ulong key, TV value)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            int index = count;
            count++;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;
            this.entries[index].value = value;

            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        public bool ContainsKey(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.lengthBuckets;
            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key) 
                    return true;
            }

            return false;
        }

        public bool TryGetValue(ulong key, out TV value)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                {
                    value = this.entries[i].value;
                    return true;
                }
            }

            value = default(TV);
            return false;
        }

        // not used but could be useful in the future. All previous calls have been replaced by (faster) in-line code
        private int FindEntry(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                {
                    return i;
                }
            }

            return -1;
        }

        public void Resize()
        {
            int newBucketsLength = lengthBuckets * 5 / 4;
            int[] newBuckets = new int[newBucketsLength];
            for (int i = 0; i < newBuckets.Length; i++)
            {
                newBuckets[i] = -1;
            }
            Entry[] newEntries = new Entry[this.lengthEntries * 5 / 4];
            Array.Copy(this.entries, 0, newEntries, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                int bucket = HashHelpers.HashMer(newEntries[j].key, rightShiftBases) % newBucketsLength;
                newEntries[j].next = newBuckets[bucket];
                newBuckets[bucket] = j;
            }
            this.buckets = newBuckets;
            this.entries = newEntries;
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            this.resizeTriggerPoint = lengthEntries;
        }

        public void Optimise()
        {
            Entry[] newEntries = new Entry[count];
            int entryIdx = 0;

            for (int bi = 0; bi < lengthBuckets; bi++)
            {
                int head = buckets[bi];
                if (head == -1)
                    continue;

                int previous = -1;
                for (int i = head; i >= 0; i = this.entries[i].next)
                {
                    newEntries[entryIdx] = entries[i];
                    newEntries[entryIdx].next = previous;
                    buckets[bi] = entryIdx;
                    previous = entryIdx;
                    entryIdx++;
                }
            }

            entries = newEntries;
            lengthEntries = entries.Length;
            newEntries = null;
        }

        public void Clear()
        {
            if (count > 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, count);
                count = 0;
            }
        }

        public IEnumerator<KeyValuePair<ulong, TV>> GetEnumerator()
        {
            for (int i = 0; i < count; i++)
            {
                KeyValuePair<ulong, TV> kvp = new KeyValuePair<ulong, TV>(entries[i].key, entries[i].value);  
                yield return kvp;
            }
        }

        public int Count
        {
            get
            {
                return (this.count);
            }
        }

        public int Capacity
        {
            get
            {
                return this.entries.Length;
            }
        }

        public TV this[ulong key]
        {
            get
            {
                int hashCode = HashHelpers.HashMer(key, rightShiftBases);
                for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
                {
                    if (this.entries[i].key == key)
                    {
                        return this.entries[i].value;
                    }
                }

                throw (new ApplicationException("not found"));
            }
            set
            {
                int hashCode = HashHelpers.HashMer(key, rightShiftBases);
                int targetBucket = hashCode % this.buckets.Length;

                // search per-bucket linked list for entry (and index)
                for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
                {
                    if (this.entries[i].key == key)
                    {
                        // found the entry so just replace its value
                        this.entries[i].value = value;
                        return;
                    }
                }

                throw (new ApplicationException("not found"));
            }
        }
    }

    // =====================================
    // -------------- kMerHashSet ----------
    // =====================================

    // A specialised read-only HashSet for holding k-mers.
    // Tables are loaded during initialisation and from then on are read-only. Loading is single-threaded, and code is thread-safe afterwards.
    // 

    public class kMerHashSet
    {
        private int[] buckets;                              // hash index array. Contains indexes into starts of linked lists for each hash bucket
        private Entry[] entries;                            // hash table contents. An array of entries linked into lists
        private int lengthBuckets = 0;                      // current length of buckets array
        private int lengthEntries = 0;                      // and the length of the entries array
        private int count;                                  // no. of allocated (in-use) entries
        private int resizeTriggerPoint;                     // force a resize when the table gets this full
        private int rightShiftBases;                        // left adjust kMer by this many bits when hashing
        private int merSize;                                // length of kMers held in the table

        public struct Entry
        {
            public int next;                                // next list pointer for this entry list
            public ulong key;                               // the mer being stored
        }

        // constructor. Must always pass in the capacity. No () constructor
        public kMerHashSet(int capacity, int merSize)
        {
            if (capacity > 0)
                Initialize(capacity, merSize);
            else
                throw (new ApplicationException("mer hashset size <= 0"));
        }

        private void Initialize(int capacity, int merSize)
        {
            this.buckets = new int[capacity + capacity / 3];
            for (int i = 0; i < this.buckets.Length; i++)
            {
                this.buckets[i] = -1;
            }
            this.entries = new Entry[capacity];
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            this.resizeTriggerPoint = capacity;
            //this.rightShiftBases = merSize > 16 ? 64 - merSize * 2 : 0;     // right-adjust the low half of the kMer if there is one 
            this.rightShiftBases = 64 - merSize * 2;        // right-adjust the kMer for hashing
            this.merSize = merSize;
        }

        // A lock-free add to a hashset (which is guaranteed not to have deletions or concurrent callers). 
        // No check is made for duplicate values - it is the callers responsibility to check for their presence before calling Add
        //
        public void AddNoCheck(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            count++;
            int index = count - 1;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;

            this.buckets[targetBucket] = index;     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        public void Add(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                    throw (new ApplicationException("mer already present"));
            }

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            count++;
            int index = count - 1;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;

            this.buckets[targetBucket] = index;                     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        public void AddIfNotPresent(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            int targetBucket = hashCode % this.buckets.Length;

            for (int i = this.buckets[targetBucket]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                    return;
            }

            // no need to search the free list as there isn't one (no deletions). Just get the next unused entry.
            count++;
            int index = count - 1;

            this.entries[index].next = this.buckets[targetBucket];  // point this new entry's next to the last entry pointed by buckets[index]
            this.entries[index].key = key;

            this.buckets[targetBucket] = index;                     // if collision occurs, update value in buckets[index] to point to new slot in entries[]

            if (count == resizeTriggerPoint)
                Resize();
        }

        public bool Contains(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);

            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                    return true;
            }

            return false;
        }

        // not used but could be useful in the future. All previous calls have been replaced by (faster) in-line code
        private int FindEntry(ulong key)
        {
            int hashCode = HashHelpers.HashMer(key, rightShiftBases);
            for (int i = this.buckets[hashCode % this.lengthBuckets]; i >= 0; i = this.entries[i].next)
            {
                if (this.entries[i].key == key)
                {
                    return i;
                }
            }

            return -1;
        }

        public void Resize()
        {
            int lengthNewBuckets = lengthBuckets * 5 / 4;
            int[] newBuckets = new int[lengthNewBuckets];
            for (int i = 0; i < lengthNewBuckets; i++)
            {
                newBuckets[i] = -1;
            }
            Entry[] newEntries = new Entry[this.lengthEntries * 5 / 4];
            Array.Copy(this.entries, 0, newEntries, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                int bucket = HashHelpers.HashMer(newEntries[j].key, rightShiftBases) % lengthNewBuckets;
                newEntries[j].next = newBuckets[bucket];
                newBuckets[bucket] = j;
            }
            this.buckets = newBuckets;
            this.entries = newEntries;
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            this.resizeTriggerPoint = lengthEntries;
        }

        public void Optimise()
        { 
            int lengthNewBuckets = count + count / 2;
            int[] newBuckets = new int[lengthNewBuckets];
            for (int i = 0; i < lengthNewBuckets; i++)
            {
                newBuckets[i] = -1;
            }
            Entry[] newEntries = new Entry[count];
            Array.Copy(this.entries, 0, newEntries, 0, this.count);
            for (int j = 0; j < this.count; j++)
            {
                int bucket = HashHelpers.HashMer(newEntries[j].key, rightShiftBases) % lengthNewBuckets;
                newEntries[j].next = newBuckets[bucket];
                newBuckets[bucket] = j;
            }
            this.buckets = newBuckets;
            this.entries = newEntries;
            this.lengthBuckets = buckets.Length;
            this.lengthEntries = entries.Length;
            //Entry[] newEntries = new Entry[count];
            //int entryIdx = 0;

            //for (int bi = 0; bi < lengthBuckets; bi++)
            //{
            //    int head = buckets[bi];
            //    if (head == -1)
            //        continue;

            //    int previous = -1;
            //    for (int i = head; i >= 0; i = this.entries[i].next)
            //    {
            //        newEntries[entryIdx].next = previous;
            //        newEntries[entryIdx].key = entries[i].key;
            //        buckets[bi] = entryIdx;
            //        previous = entryIdx;
            //        entryIdx++;
            //    }
            //}

            //entries = newEntries;
            //lengthEntries = entries.Length;
        }

        public void AddSubVariants()
        {
            List<ulong> merVariants = new List<ulong>(200);

            int initialCount = this.count;
            for (int i = 0; i < initialCount; i++)
            {
                ulong mer = entries[i].key;
                merVariants.Clear();

                int vNo = kMers.GenerateMerLastSubVariants(mer, merVariants, merSize);
                for (int v = 0; v < vNo; v++)
                {
                    ulong merVariant = merVariants[v];
                    AddIfNotPresent(merVariant);
                }
            }
        }

        public void Clear()
        {
            if (count > 0)
            {
                for (int i = 0; i < buckets.Length; i++)
                    buckets[i] = -1;
                Array.Clear(entries, 0, count);
                count = 0;
            }
        }

        public IEnumerator<ulong> GetEnumerator()
        {
            for (int i = 0; i < count; i++)
            {
                yield return entries[i].key;
            }
        }

        public void CopyEntries(int starting, ulong[] entriesArray)
        {
            for (int i = 0; i < count; i++)
                entriesArray[starting + i] = entries[i].key;
        }

        public int Count
        {
            get
            {
                return (this.count);
            }
        }

        public int Capacity
        {
            get
            {
                return this.entries.Length;
            }
        }

        public Dictionary<int, int> CheckHashTable()
        {
            Dictionary<int, int> chains = new Dictionary<int, int>();

            int inUseBuckets = 0;
            for (int bi = 0; bi < lengthBuckets; bi++)
                if (buckets[bi] >= 0)
                    inUseBuckets++;
            chains.Add(-1, inUseBuckets);

            for (int bi = 0; bi < lengthBuckets; bi++)
            {
                int bucket = buckets[bi];
                if (bucket == -1)
                    continue;

                int chainLength = 1;

                for (int i = this.entries[bucket].next; i >= 0; i = this.entries[i].next)
                {
                    chainLength++;
                }

                if (chains.ContainsKey(chainLength))
                    chains[chainLength]++;
                else
                    chains.Add(chainLength, 1);
            }

            return chains;
        }
    }

    internal static class HashHelpers
    {
        // this class used to contain code to find prime numbers for the sizing of bucket arrays. Prime table sizes are only a defence
        // against poor hashing functions (and unwanted periodicity) and are not really needed.

        internal static int HashMer(ulong key, int rightShiftBases)
        {
            // experiments with optimised code showed that the simple folded hash performed almost as well for randomness as the Murmur hash, and was short enough
            // to be inlined. The overall result was that the folded hash was faster. Copying the body of this function in-line made no difference. 

            const int positiveInt32Mask = 0x7fffffff;

            key = key >> rightShiftBases;
            ulong foldedUL = (key >> 32) ^ (key & 0xffffffff);
            return (int)(foldedUL) & positiveInt32Mask;

            //return (int)((key >> 32) ^ (key & 0xffffffff)) & positiveInt32Mask;
            //return (int)((key >> 32) ^ ((key >> rightShiftBases) & 0xffffffff)) & positiveInt32Mask;
        }
     }
}
