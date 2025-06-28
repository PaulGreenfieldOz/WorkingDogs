using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace WorkingDogsCore
{
	// =======================================
    // -------------- Sequence ---------------
    // =======================================
    public class Sequence
    {
        public int Length;
        public int Capacity;
        public char[] Bases;

        public Sequence(int capacity)
        {
            Length = 0;
            Capacity = capacity;
            Bases = new char[Capacity];
        }

        public Sequence(string s)
        {
            Length = s.Length;
            Capacity = s.Length + 10;
            Bases = new char[Capacity];
            s.CopyTo(0, Bases, 0, s.Length);
        }

        public Sequence(string s, int capacity)
        {
            Length = s.Length;
            Capacity = capacity;
            Bases = new char[Capacity];
            s.CopyTo(0, Bases, 0, s.Length);
        }

        public Sequence(Sequence s)
        {
            Length = s.Length;
            Capacity = s.Capacity;
            Bases = new char[Capacity];
            Array.Copy(s.Bases, 0, Bases, 0, this.Length);
        }

        public Sequence(Sequence s, int capacity)
        {
            Length = s.Length;
            Capacity = capacity;
            Bases = new char[Capacity];
            Array.Copy(s.Bases, 0, Bases, 0, this.Length);
        }

        public void Resize(int newCapacity)
        {
            Array.Resize<char>(ref this.Bases, newCapacity);
            this.Capacity = newCapacity;
        }

        public bool Matches(string s)
        {
            if (s.Length != Length)
                return false;

            for (int i = 0; i < Length; i++)
            {
                if (Bases[i] != s[i])
                    return false;
            }

            return true;
        }

        public new string ToString()
        {
            return new string(Bases, 0, Length);
        }

        public void ToUpper()
        {
            for (int i = 0; i < Length; i++)
                Bases[i] = Char.ToUpper(Bases[i]);
        }

        public string ToString(int start, int length)
        {
            return new string(Bases, start, length);
        }

        public void Trim(int start, int length)
        { 
            for (int i = 0; i < length; i++)
                Bases[i] = Bases[start + i];
            this.Length = length;
        }

        public Sequence SubSeq(int start, int length)
        {
            if (start + length > this.Length)
                length = this.Length - start;
            Sequence newSeq = new Sequence(length);

            for (int i = 0; i < length; i++)
                newSeq.Bases[i] = Bases[start + i];
            newSeq.Length = length;

            return newSeq;
        }
        public string SubString(int start, int length)
        {
            if (start + length > this.Length)
                length = this.Length - start;

            return new string(this.Bases, start, length);
        }
        public void Reverse()
        {
            int halfWay = this.Length / 2;
            for (int i = 0; i < halfWay; i++)
            {
                int rhsIdx = this.Length - 1 - i;
                char savedChar = this.Bases[rhsIdx];
                this.Bases[rhsIdx] = Bases[i];
                this.Bases[i] = savedChar;
            }
        }

        public void ReverseComplement()
        {
            int halfWay = this.Length / 2;
            for (int i = 0; i < halfWay; i++)
            {
                int rhsIdx = this.Length - 1 - i;
                char savedChar = this.Bases[rhsIdx];
                this.Bases[rhsIdx] = this.Bases[i];
                this.Bases[i] = savedChar;
            }
            for (int i = 0; i < this.Length; i++)
            {
                int b = (int)kMers.BaseCharToInt(this.Bases[i]);
                switch (b)
                {
                    case 0:
                        this.Bases[i] = 'T';
                        break;
                    case 1:
                        this.Bases[i] = 'G';
                        break;
                    case 2:
                        this.Bases[i] = 'C';
                        break;
                    case 3:
                        this.Bases[i] = 'A';
                        break;
                    case -1:
                        this.Bases[i] = 'N';
                        break;
                }
            }
        }

        public void Append(char b)
        {
            if (Capacity - Length < 10)
            {
                Capacity += 50;
                Array.Resize(ref Bases, Capacity);
            }

            Bases[Length] = b;
            Length++;
        }
        public void Append(string s)
        {
            if (Length == Capacity)
            {
                Capacity += 50 + s.Length;
                Array.Resize(ref Bases, Capacity);
            }

            for (int b = 0; b < s.Length; b++)
            {
                Bases[Length] = s[b];
                Length++;
            }
        }

        public void CopyTo(Sequence copy)
        {
            if (this.Capacity > copy.Capacity)
            {
                Array.Resize<char>(ref copy.Bases, this.Capacity);
                copy.Capacity = this.Capacity;
            }
            copy.Length = this.Length;
            Array.Copy(this.Bases, 0, copy.Bases, 0, this.Length);
        }

        public void CopyFrom(string source)
        {
            if (source == null)
            {
                this.Length = 0;
                return;
            }

            if (source.Length > this.Capacity)
            {
                Array.Resize<char>(ref this.Bases, source.Length);
                this.Capacity = source.Length;
            }

            this.Length = source.Length;
            source.CopyTo(0, this.Bases, 0, source.Length);
        }

        public void Insert(int m, char c)
        {
            if (this.Length == this.Capacity)
            {
                Array.Resize(ref this.Bases, this.Capacity + 20);
                this.Capacity += 20;
            }

            for (int i = this.Length; i > m; i--)
                this.Bases[i] = this.Bases[i - 1];

            this.Bases[m] = c;
            this.Length++;
        }

        public void Replace(int m, ulong packedMer, int merSize)
        {
            ulong intBase = 0;                              // current base in binary form
            ulong tempMer = packedMer;                      // progressively shifted left
            char ACGT;                                      // char form of base

            for (int i = 0; i < merSize; i++)
            {
                intBase = (ulong)tempMer >> 62;
                tempMer = (ulong)(tempMer << 2);

                ACGT = kMers.baseToChar[intBase];
                this.Bases[m + i] = ACGT;
            }
        }

        public void Remove(int m)
        {
            for (int i = m; i < this.Length - 1; i++)
                this.Bases[i] = this.Bases[i + 1];
            this.Length--;
        }

        public void Remove(int m, int n)
        {
            for (int i = m; i < this.Length - n; i++)
                this.Bases[i] = this.Bases[i + n];
            this.Length -= n;
        }

        public int IndexOf(char c)
        {
            int idx = -1;
            for (int i = 0; i < this.Length; i++)
                if (this.Bases[i] == c)
                {
                    idx = i;
                    break;
                }

            return idx;
        }

        public static bool CondenseMer(Sequence seq, int start, int merSize, out ulong packedMer)
        {
            packedMer = 0;
            bool allBasesACGT = true;

            for (int m = 0; m < merSize; m++)
            {
                char nextBase = seq.Bases[m + start];
                long packedBase = kMers.BaseCharToInt(nextBase);

                if (packedBase < 0)
                {
                    packedBase = 0;
                    allBasesACGT = false;
                }

                packedMer = (packedMer << 2) | (ulong)packedBase;
            }
            packedMer = packedMer << (64 - merSize * 2);
            return allBasesACGT;
        }

        public static bool CondenseMerIncremental(int merSize, ulong previousMer, Sequence read, int m, out ulong nextMer)
        {
            char nextBase = read.Bases[m + merSize - 1];
            long packedBase = kMers.BaseCharToInt(nextBase);
            nextMer = (previousMer << 2) | (ulong)packedBase << (64 - merSize * 2);
            if (packedBase < 0)
            {
                nextMer = 0;
                return false;
            }
            return true;
        }

        public static void ReverseComplement(Sequence s)
        {
            int sl = s.Length;
            char[] rcs = new char[sl];

            for (int i = 0; i < sl; i++)
                rcs[sl - i - 1] = kMers.baseToComplement[s.Bases[i] & 0x1f];

            rcs.CopyTo(s.Bases, 0);
        }

        public static int GenerateMersFromRead(Sequence read, int merSize, ref ulong[] mers, ref bool[] merValid)
        {
            int mersInRead = read.Length - merSize + 1;
            bool merIsValid = false;
            ulong mer = 0;

            if (mersInRead < 1)
                return 0;

            if (mers.Length < mersInRead)
            {
                Array.Resize<ulong>(ref mers, mersInRead + 100);
                Array.Resize<bool>(ref merValid, mersInRead + 100);
            }

            for (int i = 0; i < mersInRead; i++)
            {
                if (merIsValid)
                {
                    merIsValid = CondenseMerIncremental(merSize, mer, read, i, out mer);
                    merValid[i] = merIsValid;
                    mers[i] = mer;
                }
                else
                {
                    merIsValid = CondenseMer(read, i, merSize, out mer);
                    merValid[i] = merIsValid;
                    mers[i] = mer;
                }
            }

            return mersInRead;
        }

        public int ComputeHash()
        {
            // not used for now. Could possibly do with a better hash algorithm as well. 
            unchecked
            {
                const int p = 16777619;
                int hash = (int)2166136261;

                for (int i = 0; i < this.Bases.Length; i++)
                    hash = (hash ^ this.Bases[i]) * p;

                hash += hash << 13;
                hash ^= hash >> 7;
                hash += hash << 3;
                hash ^= hash >> 17;
                hash += hash << 5;
                return hash;
            }
        }
    }
}