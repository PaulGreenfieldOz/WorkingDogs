using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace WorkingDogsCore
{
    public class SeqFiles
    {
        public const int formatNone = 0;            // no format specified
        public const int formatFASTA = 1;			// fasta with multiple data lines (.fna, .fa, .fasta, .fas, ...)
        public const int formatFNA = 1;             // synonym
        public const int formatFASTQ = 2;           // fastq 

        enum Formats
        {
            None = 0,
            FASTA = 1,
            FASTQ = 2
        }

        static char[] spaceDelimiter = new char[] { ' ' };

        public static string ReadRead(StreamReader reads, int readFormat)
        {
            string read = null;
            string readHeader;

            if (readFormat == formatFASTQ)
            {
                readHeader = reads.ReadLine();
                if (readHeader == null)
                    return null;
                read = reads.ReadLine();
                reads.ReadLine();
                reads.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out readHeader);
            }

            return read;
        }

        public static string ReadRead(StreamReader reads, int readFormat, out string header)
        {
            string read = null;
            header = "";

            if (readFormat == formatFASTQ)
            {
                header = reads.ReadLine();
                if (header == null)
                    return null;
                read = reads.ReadLine();
                reads.ReadLine();
                reads.ReadLine();
            }
            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out header);
            }

            return read;
        }

        public static string ReadRead(StreamReader reads, int readFormat, out string readHeader, out string qualHeader, out string quals)
        {
            string read = null;
            readHeader = null;
            qualHeader = null;
            quals = null;

            if (readFormat == formatFASTQ)
            {
                readHeader = reads.ReadLine();
                if (readHeader == null)
                    return null;
                read = reads.ReadLine();
                qualHeader = reads.ReadLine();
                quals = reads.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                read = ReadFASTA(reads, false, out readHeader);
            }

            return read;
        }

        public static string ReadRead(StreamReader readsFile, StreamReader qualsFile, int readFormat, out string readHeader, out string qualHeader, out string quals)
        {
            string read = null;
            readHeader = null;
            qualHeader = null;
            quals = null;

            if (readFormat == formatFASTQ)
            {
                readHeader = readsFile.ReadLine();
                if (readHeader == null)
                    return null;
                read = readsFile.ReadLine();
                qualHeader = readsFile.ReadLine();
                quals = readsFile.ReadLine();
            }

            if (readFormat == formatFNA)
            {
                try
                {
                    read = ReadFASTA(readsFile, false, out readHeader);
                }
                catch (Exception e)
                {
                    Console.WriteLine("exception from ReadFASTA: " + e.Message);
                    throw e;
                }
                if (qualsFile != null)
                {
                    quals = ReadFASTA(qualsFile, true, out qualHeader);
                    if (quals != null)
                    {
                        Sequence qualsSeq = new Sequence(quals);
                        ConvertQualsToCanonicalForm(qualsSeq, readFormat, 0);
                        quals = qualsSeq.ToString();
                    }
                }
            }

            return read;
        }

        public static bool ReadRead(StreamReader reads, StreamReader quals, int readFormat, Sequence readHeader, Sequence read, Sequence qualHeader, Sequence qual)
        {
            string readString = null;
            string readHeaderString = null;
            string qualHeaderString = null;
            string qualString = null;

            readString = ReadRead(reads, quals, readFormat, out readHeaderString, out qualHeaderString, out qualString);

            if (readString != null)
            {
                readHeader.CopyFrom(readHeaderString);
                read.CopyFrom(readString);
                qualHeader.CopyFrom(qualHeaderString);
                qual.CopyFrom(qualString);
            }
            else
            {
                readHeader.Length = 0;
                read.Length = 0;
                qualHeader.Length = 0;
                qual.Length = 0;
            }

            return read.Length > 0;
        }

        public static int ReadReads(int readsWanted, StreamReader readsFile, StreamReader qualsFile, int readsFormat, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals)
        {
            int readsRead = 0;

            for (int i = 0; i < readsWanted; i++)
            {
                string readHeader;
                string readSeq;
                string qualHeader;
                string qual;

                readSeq = ReadRead(readsFile, qualsFile, readsFormat, out readHeader, out qualHeader, out qual);

                if (readSeq == null)
                    break;

                readHeaders[i].CopyFrom(readHeader);
                readSeqs[i].CopyFrom(readSeq);
                if (qualHeaders != null)
                {
                    qualHeaders[i].CopyFrom(qualHeader);
                    quals[i].CopyFrom(qual);
                }

                readsRead++;
            }

            return readsRead;
        }

        // Read a FASTA read. The 'addSpace' parameter controls whether the lines of data are concatenated with a gap. For base data, no gap is wanted
        // but the same function is called for quals (for the 454 separate-qual format) and a space is needed to separate the two adjacent quals.
        private static string ReadFASTA(StreamReader reads, bool addSpace, out string readHeader)
        {
            readHeader = "";
            string currentLine = reads.ReadLine();          // read the next line from the reads file
            if (currentLine == null)                        // and return if we're at EOF
                return null;

            if (currentLine[currentLine.Length - 1] == ' ')
                currentLine = currentLine.TrimEnd();

            readHeader = currentLine;

            // multiple file lines per read
            StringBuilder currentRead = new StringBuilder(1000);
            bool keepReading = true;
            int peekChar = 0;
            currentLine = reads.ReadLine();

            while (keepReading)
            {
                if (currentLine.Length != 0)
                {
                    currentRead.Append(currentLine);        // add the just read line to the read under construction
                    if (addSpace)
                        currentRead.Append(' ');            // leave space between lines if needed
                }

                peekChar = reads.Peek();                    // and peek at the next line to see if we're at the last line for this current read
                if (peekChar < 0 || peekChar == '>')
                    keepReading = false;                    // yes... either EOF or the header for the next read
                else
                {
                    currentLine = reads.ReadLine();         // no... next line is a read line
                    if (currentLine.Length != 0 && currentLine[currentLine.Length - 1] == ' ')
                        currentLine = currentLine.TrimEnd();
                }
            }

            return currentRead.ToString();
        }

        public static void ConvertQualsToCanonicalForm(Sequence quals, int readsFormat, int qualBase)
        {
            if (readsFormat == formatFASTQ)
            {
                for (int i = 0; i < quals.Length; i++)
                    quals.Bases[i] -= (char)qualBase;
                return;
            }

            if (readsFormat == formatFNA && quals.Length != 0)
            {
                string qualString = quals.ToString();
                string[] qualStrings = qualString.Split(spaceDelimiter, StringSplitOptions.RemoveEmptyEntries);
                for (int i = 0; i < qualStrings.Length; i++)
                    quals.Bases[i] = (char)(Convert.ToInt16(qualStrings[i]));
                quals.Length = qualStrings.Length;
            }
        }

        public static void WriteRead(StreamWriter readsFile, StreamWriter qualsFile, string header, string read, string qualHeader, string quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header);
                readsFile.WriteLine(read);
                readsFile.WriteLine(qualHeader);
                readsFile.WriteLine(quals);
            }
            else
            {
                WriteRead(readsFile, header, read, readsFormat);
                WriteQual(qualsFile, qualHeader, quals, readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, string header, string read, string qualHeader, string quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header);
                readsFile.WriteLine(read);
                readsFile.WriteLine(qualHeader);
                readsFile.WriteLine(quals);
            }
            else
            {
                WriteRead(readsFile, header, read, readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, Sequence header, Sequence read, Sequence qualHeader, Sequence quals, int readsFormat)
        {
            if (readsFormat == formatFASTQ)
            {
                readsFile.WriteLine(header.Bases, 0, header.Length);
                readsFile.WriteLine(read.Bases, 0, read.Length);
                readsFile.WriteLine(qualHeader.Bases, 0, qualHeader.Length);
                readsFile.WriteLine(quals.Bases, 0, quals.Length);
            }
            else
            {
                WriteRead(readsFile, header.ToString(), read.ToString(), readsFormat);
            }
        }

        public static void WriteRead(StreamWriter readsFile, string readHeader, string read, int readFormat)
        {
            // @1:1:0:686#0/1
            // NTGGAGAATTCTGGATCCTCGGACTAAAACAATAGCAGTTGATTCGCTCACAGTTCTGGAGGCTAGAGGTATGAAA
            // +1:1:0:686#0/1
            // @hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\hhhhhhhhhhhhhhhhVg[hhhU^hhWhfgVhc^Dhh`_V
            //
            // >FUOCO2J02IPX2K rank=0000260 x=3458.0 y=3626.0 length=70
            // GAGCAGCTCCTCAAGCAACTGCAACTGGATTGAGCAAGGAATGTTCGCAGCTACCCGACT
            // GACCCGTCTT

            if (readsFile == null)
                return;

            //// adjust length for 454 reads
            //if (readFormat == formatFNA && readHeader != null && readHeader != "")
            //{
            //    int lengthIdx = readHeader.IndexOf("length=");
            //    if (lengthIdx > 0)
            //    {
            //        int numIdx = lengthIdx + "length=".Length;
            //        readHeader = readHeader.Substring(0, numIdx) + read.Length;
            //    }
            //}

            // write out the header if it exists
            if (readHeader != null)
                readsFile.WriteLine(readHeader);

            if (readFormat == formatFNA)
            {
                // if we concatenated shorter source lines (e.g. from 454 reads) when the data was read,
                // we'll now split it again as it's being written
                int m = 0;
                int hrLen = read.Length;
                while (m < hrLen)
                {
                    int wLen = Math.Min(60, hrLen - m);
                    readsFile.WriteLine(read.Substring(m, wLen));
                    m += 60;
                }
            }
            else
                readsFile.WriteLine(read);
        }

        private static void WriteQual(StreamWriter qualFile, string qualHeader, string quals, int readFormat)
        {
            if (qualFile == null)
                return;

            int qualsLength = quals.Length;

            if (readFormat == formatFNA && qualHeader != null && qualHeader != "")
            {
                int lengthIdx = qualHeader.IndexOf("length=");
                if (lengthIdx > 0)
                {
                    int numIdx = lengthIdx + "length=".Length;
                    qualHeader = qualHeader.Substring(0, numIdx) + qualsLength;
                }
            }

            if (qualHeader != null)
                qualFile.WriteLine(qualHeader);

            if (readFormat == formatFNA)
            {
                // if we concatenated shorter qual lines (e.g. from 454 reads) when they were being read,
                // we'll now reformat and split them again as they're being written
                int m = 0;
                while (m < quals.Length)
                {
                    int wLen = Math.Min(20, qualsLength - m);
                    for (int i = m; i < m + wLen; i++)
                    {
                        int qualInt = (int)quals[i];
                        qualFile.Write(qualInt);
                        qualFile.Write(' ');
                    }
                    qualFile.WriteLine();
                    m += 20;
                }
            }

            if (readFormat == formatFASTQ)
            {
                qualFile.WriteLine(quals);
            }
        }

        // The qual scores used with Fastq reads are ambiguous. Illumina data may use a base of 64 ('@') while Sanger formatted data
        // uses a base of 33 ('!'). This method reads quals from a reads file until it finds a value that can only be Sanger or until it's
        // looked at the first 100 reads, and then returns the appropriate qual offset.
        public static int ResolveFastqQualAmbiguity(string fastqReadsFN, out bool fullQualHeader)
        {
            fullQualHeader = true;
            StreamReader fastqReads = new StreamReader(fastqReadsFN);
            int qualBase = 64;                          // assume Illumina format until we find a Sanger-only value - this forces Illumina to be the default
            const char minIlluminaQual = '@';           // lowest possible Illumina qual
            const int readsToExamine = 100;             // only look at the first 100 reads
            int readsExamined = 0;
            bool foundIndicativeBase = false;           // found a base that can only appear in Sanger

            while (readsExamined < readsToExamine)
            {
                string readHeader = fastqReads.ReadLine();
                if (readHeader == null)
                {
                    Console.WriteLine("failed to resolve fastq qual ambiguity before reaching EOF on " + fastqReadsFN);
                    break;
                }
                readsExamined++;
                string read = fastqReads.ReadLine();
                string qualHeader = fastqReads.ReadLine();
                if (qualHeader == "+")
                    fullQualHeader = false;
                string quals = fastqReads.ReadLine();

                for (int i = 0; i < quals.Length; i++)
                {
                    if (quals[i] < minIlluminaQual)         // lower than any possible Illumina qual --> Sanger formatted data
                    {
                        qualBase = 33;
                        // Console.WriteLine("Sanger");
                        foundIndicativeBase = true;
                        break;
                    }
                }

                if (foundIndicativeBase)
                    break;
            }
            fastqReads.Close();

            return qualBase;
        }
   
        // finds any trailing poor qual region in a Sequence read and trims the read back to good bases
        public static int TrimTrailingPoorQuals(Sequence read, Sequence quals, int trimQual, int qualOffset)
        {
            const int windowLength = 10;
            const int passesNeededInWindow = windowLength * 3 / 4;
            int basesTrimmed = 0;
            char[] qualCharArray = quals.Bases;
            int goodBaseIdx = 0;

            for (int i = quals.Length - 1; i >= 0; i--)
                if (((int)qualCharArray[i] - qualOffset) > trimQual)
                {
                    goodBaseIdx = i;
                    break;
                }

            basesTrimmed = quals.Length - goodBaseIdx;

            for (int i = goodBaseIdx; i > windowLength; i--)
            {
                int passedInWindow = 0;
                for (int w = 0; w < windowLength; w++)
                {
                    int qualInt = (int)qualCharArray[i - w] - qualOffset;
                    if (qualInt >= trimQual)
                        passedInWindow++;
                }
                if (passedInWindow < passesNeededInWindow)
                    basesTrimmed++;
                else
                    break;
            }

            if (basesTrimmed > 0)
            {
                read.Length = read.Length - basesTrimmed;
                quals.Length = quals.Length - basesTrimmed;
            }

            return basesTrimmed;
        }

        // finds the trailing poor qual region in a read and returns the length that should be trimmed to get back to more robust data
        public static int FindTrailingPoorQuals(string read, string quals, int trimQual, int qualOffset)
        {
            const int windowLength = 10;
            const int passesNeededInWindow = windowLength * 3 / 4;
            int basesToTrim = 0;
            int goodBaseIdx = 0;

            int[] qualScores = new int[quals.Length];
            for (int i = 0; i < quals.Length; i++)
                qualScores[i] = (int)quals[i] - qualOffset;

            for (int i = quals.Length - 1; i >= 0; i--)
                if (((int)quals[i] - qualOffset) > trimQual)
                {
                    goodBaseIdx = i;
                    break;
                }

            basesToTrim = quals.Length - goodBaseIdx;

            for (int i = goodBaseIdx; i > windowLength; i--)
            {
                int passedInWindow = 0;
                for (int w = 0; w < windowLength; w++)
                {
                    int qualInt = (int)quals[i - w] - qualOffset;
                    if (qualInt >= trimQual)
                        passedInWindow++;
                    if (read[i - w] == 'N')
                        passedInWindow--;
                }
                if (passedInWindow < passesNeededInWindow)
                    basesToTrim++;
                else
                    break;
            }

            return basesToTrim;
        }

        public static string LFConvention(string readsFN)
        {
            StreamReader reads = new StreamReader(readsFN);
            bool foundCR = false;
            bool foundLF = false;

            while (!foundLF)
            {
                int c = reads.Read();
                if (c == '\r')
                    foundCR = true;
                if (c == '\n')
                    break;
            }

            reads.Close();
            return foundCR ? "\r\n" : "\n";
        }

        public static int DetermineFileFormat(string fn)
        {
            int fileFormat = formatNone;
            StreamReader file = new StreamReader(fn);
            string firstLine = file.ReadLine();
            if (firstLine != null)
            {
                if (firstLine[0] == '@')
                    fileFormat = formatFASTQ;
                if (firstLine[0] == '>')
                    fileFormat = formatFNA;
            }
            return fileFormat;
        }

        public static void SplitFileName(string mappingFNP, out string dirPart, out string fnPart)
        {
            // split the file name/pattern if necessary
            int lastFS = mappingFNP.LastIndexOf(Path.DirectorySeparatorChar);
            fnPart = mappingFNP;
            dirPart = Directory.GetCurrentDirectory();
            if (lastFS >= 0)
            {
                fnPart = mappingFNP.Substring(lastFS + 1);
                dirPart = mappingFNP.Substring(0, lastFS);
            }
        }

        public static bool ValidateOutputDir(string outputDir)
        {
            bool directoryOK = true;

            try // in case the Exists/CreateDirectory calls fail
            {
                if (outputDir[outputDir.Length - 1] == Path.DirectorySeparatorChar)
                    outputDir = outputDir.Substring(0, outputDir.Length - 1);
                if (!Directory.Exists(outputDir))
                    Directory.CreateDirectory(outputDir);
                string testOutputFN = outputDir + Path.DirectorySeparatorChar + "43EDD23F-5F68-47f0-B7B9-66AE9EE3BF0B.txt";
                StreamWriter testTemp = new StreamWriter(testOutputFN);
                testTemp.Close();
                File.Delete(testOutputFN);
            }
            catch
            {
                directoryOK = false;
            }

            return directoryOK;
        }
    }

    // =======================================
    // ----------- BufferedReader ------------
    // =======================================

    public class BufferedReader
    {
        public int fileFormat = SeqFiles.formatNone;        // fasta or fastq
        StreamReader readsFile = null;                      // fastq: seqs+quals; fasta: seqs only
        StreamReader qualsFile = null;                      // fasta only: separate quals file (rare and not buffered)

        FillBufferDelegate fbd = null;                      // async FillBuffer delegate
        IAsyncResult iarFillBuffer;                         // and the result of calling FillBuffer
        char[][] buffers = new char[2][];                   // a pair of buffers - one being filled/full and one being emptied
        const int bufferLength = 1000000;                   // length of these buffers
        int[][] lineLengths = new int[2][];                 // lengths of the lines found in the currently-being-emptied buffer (including \n)
        int[] endOfLastCompleteRead = new int[2];           // index of the start of the last complete read in each buffer
        int currentBufferIdx;                               // pointers to current and previous buffers (0 or 1)
        int previousBufferIdx;                              // swapped when buffers swap roles
        bool EOF;                                           // EOF reached on readsFile - stop calling FillBuffer
        bool readsFromBuffer;                               // buffered reads (usual) or direct calls (fasta with separate quals)
        bool bufferEmpty = true;                            // previous ReadRead call returned the last full read - time to swap buffers 
        int nextLineIdx = 0;                                // next lineLength index to be used
        int startOfNextRead = 0;                            // index into buffer of start of next line

        public BufferedReader(int fileFormat, StreamReader readsFile, StreamReader qualsFile)
        {
            this.fileFormat = fileFormat;
            this.readsFile = readsFile;
            this.qualsFile = qualsFile;

            readsFromBuffer = fileFormat == SeqFiles.formatFASTQ || (fileFormat == SeqFiles.formatFNA && qualsFile == null);

            if (readsFromBuffer)
            {
                buffers[0] = new char[bufferLength+1];          // in case we need to insert a \n at the end of the final buffer
                buffers[1] = new char[bufferLength+1];
                lineLengths[0] = new int[bufferLength / 20];    // 5% 
                lineLengths[1] = new int[bufferLength / 20];

                currentBufferIdx = 0;
                previousBufferIdx = 1;
                EOF = false;
                // start filling the initial buffer
                fbd = new FillBufferDelegate(FillBuffer);
                iarFillBuffer = fbd.BeginInvoke(readsFile, previousBufferIdx, currentBufferIdx, out EOF, null, null);
            }
        }

        public void Close()
        {
            if (iarFillBuffer != null && !EOF)
            {
                // wait for any in-process buffer fill to complete 
                if (!iarFillBuffer.IsCompleted)
                    iarFillBuffer.AsyncWaitHandle.WaitOne();
                // get out parameters from the async call
                fbd.EndInvoke(out EOF, iarFillBuffer);
                // and reset the wait handle
                iarFillBuffer.AsyncWaitHandle.Close();
            }
            readsFile.Close();
            if (qualsFile != null)
                qualsFile.Close();
        }

        public int ReadReads(int readsWanted, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals)
        {
            int readsIdx = 0;

            if (readsFromBuffer)
            {
                while (true)
                {
                    if (bufferEmpty)
                    {
                        if (EOF)
                            break;

                        // wait for previous buffer fill to complete
                        if (!iarFillBuffer.IsCompleted)
                            iarFillBuffer.AsyncWaitHandle.WaitOne();

                        // get out parameters from the async call
                        fbd.EndInvoke(out EOF, iarFillBuffer);
                        // and reset the wait handle
                        iarFillBuffer.AsyncWaitHandle.Close();

                        // and start the next buffer read
                        if (!EOF)
                        {
                            iarFillBuffer = fbd.BeginInvoke(readsFile, currentBufferIdx, previousBufferIdx, out EOF, null, null);
                        }

                        previousBufferIdx = currentBufferIdx;
                        if (currentBufferIdx == 0)
                            currentBufferIdx = 1;
                        else
                            currentBufferIdx = 0;

                        nextLineIdx = 0;
                        startOfNextRead = 0;
                    }

                    bufferEmpty = GetNextRead(currentBufferIdx, readHeaders[readsIdx], readSeqs[readsIdx], qualHeaders[readsIdx], quals[readsIdx]);

                    readsIdx++;
                    if (readsIdx == readsWanted)
                        break;
                }
            }
            else
            {
                // if we're dealing with the rare file type of fasta with a separate quals, lock the pair of files before
                // reading to ensure that we pick up pairs of reads/quals and just do unbuffered reads
                lock (readsFile)
                {
                    readsIdx = SeqFiles.ReadReads(readsWanted, readsFile, qualsFile, fileFormat, readHeaders, readSeqs, qualHeaders, quals);
                }
            }
            return readsIdx;
        }

        private bool GetNextRead(int cbi, Sequence readHeader, Sequence readSeq, Sequence qualHeader, Sequence quals)
        {
            char[] currentBuffer = buffers[cbi];
            int[] currentLineLengths = lineLengths[cbi];
            int snr = startOfNextRead;
            int nli = nextLineIdx;
            int currentEndOfLastCompleteRead = endOfLastCompleteRead[cbi];

            if (fileFormat == SeqFiles.formatFASTQ)
            {
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readHeader, 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readSeq, 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], qualHeader, 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], quals, 0);
                snr += currentLineLengths[nli];
                nli++;
            }

            if (fileFormat == SeqFiles.formatFNA)
            {
                bool continueCopying = true;

                // copy the header
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readHeader, 0);
                snr += currentLineLengths[nli];
                nli++;

                // find out out long the read is and resize the readseq appropriately
                int readLength = 0;
                int nliForLength = nli;
                int snrForLength = snr; 
                while (continueCopying)
                {
                    int lineLength = currentLineLengths[nliForLength];
                    readLength += lineLength;
                    snrForLength += lineLength;
                    nliForLength++;

                    if (snrForLength == currentEndOfLastCompleteRead || currentBuffer[snrForLength] == '>')
                        break;
                }

                if (readSeq.Capacity < readLength)
                    readSeq.Resize(readLength + 100);

                //if (readHeader.ToString() == ">NODE_543_length_150414_cov_144.508102" || readHeader.Bases[0] !='>')
                //    Debugger.Break();

                int copyStartIdx = 0;
                while (continueCopying)
                {
                    CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readSeq, copyStartIdx);
                    copyStartIdx = readSeq.Length;
                    snr += currentLineLengths[nli];
                    nli++;

                    if (snr == currentEndOfLastCompleteRead || currentBuffer[snr] == '>')
                        continueCopying = false;
                }
            }

            nextLineIdx = nli;
            startOfNextRead = snr;
            return startOfNextRead == currentEndOfLastCompleteRead || currentLineLengths[nli] == 0;
        }

        private void CopySeqFromBuffer(char[] buffer, int startOfSeq, int length, Sequence destination, int startIdx)
        {
            // copy up to but not including the \n
            int lengthAfterCopy = length + startIdx - 1;
            if (lengthAfterCopy > destination.Capacity)
                destination.Resize(lengthAfterCopy + lengthAfterCopy/2);
            Array.Copy(buffer, startOfSeq, destination.Bases, startIdx, length-1);      
            if (length > 1 && destination.Bases[lengthAfterCopy - 1] == '\r')
                lengthAfterCopy--;
            destination.Length = lengthAfterCopy;
        }

        private delegate void FillBufferDelegate(StreamReader reads, int cbi, int pbi, out bool EOF);

        private void FillBuffer(StreamReader reads, int cbi, int pbi, out bool EOF)
        {
            char[] currentBuffer = buffers[cbi];
            char[] previousBuffer = buffers[pbi];

            int bi = 0;

            // copy any remaining chars from the previous buffer
            int startOfLeftover = endOfLastCompleteRead[pbi];
            if (startOfLeftover > 1)
            {
                int leftOverLength = bufferLength - startOfLeftover;
                Array.Copy(previousBuffer, startOfLeftover, currentBuffer, bi, leftOverLength);
                bi = leftOverLength;
            }

            int charsToRead = bufferLength - bi;
            int charsRead = reads.ReadBlock(currentBuffer, bi, charsToRead);
            EOF = charsRead != charsToRead;
            int charsInBuffer = bi + charsRead;

            // find the line boundaries to save having to scan the buffer again when extracting actual reads
            // and also finds the end of the incomplete read at the end of the buffer
            endOfLastCompleteRead[cbi] = FindLinesInBuffer(cbi, charsInBuffer, EOF);
        }

        private int FindLinesInBuffer(int cbi, int charsInBuffer, bool EOF)
        {
            // finds all the start-of-lines (char after /n) in the buffer and remembers them in startsOfLines. These are the indices in the current
            // buffer of the first char of every line. The last entry in startsOfLines will be -1.
            // It then goes through and finds the start/end of all the reads (format-dependent).
            // The start of the last incomplete read is passed back so this region can be copied into the next buffer when it's filled.

            char[] buffer = buffers[cbi];
            int[] lineLengthsInBuffer = lineLengths[cbi];

            int li = 0;                     // start from the first line in the buffer
            int startOfNextLine = 0;        // start of the next line (starting point of scan for next \n)
            int startOfFinalRead = 0;       // char index of start of last read
            int endOfLastCompleteRead = 0;  // index of end of the last complete read in the buffer (return value)
            int fastqLineCount = 0;         // count 4 line sets to work out when a fastq line is starting
            bool lastLineComplete = true;  // \n at very end of buffer

            // force \n at end of buffer at EOF 
            if (EOF && buffer[charsInBuffer-1] != '\n')
            {
                buffer[charsInBuffer] = '\n';
                charsInBuffer++;
            }

            // find and remember all the \n locations (plus 1 to get to the start of the next line - may be incomplete)
            while (startOfNextLine < charsInBuffer)
            {
                int bi = Array.IndexOf<char>(buffer, '\n', startOfNextLine);        // find the next \n
                if (bi == -1)                                                       // no EOL at end of buffer
                {
                    lineLengthsInBuffer[li] = 0;                                    // incomplete line so don't record a length
                    lastLineComplete = false;
                    break;
                }

                // now have the next complete line - so track reads
                if (fileFormat == SeqFiles.formatFASTQ)                             // tracking fastq reads
                {
                    fastqLineCount++;                                               // keep track of where we are in fastq quad
                    if (fastqLineCount > 4)
                        fastqLineCount = 1;

                    if (fastqLineCount == 1 && buffer[startOfNextLine] == '@')      // remember each header line as it goes past
                        startOfFinalRead = startOfNextLine;
                }
                if (fileFormat == SeqFiles.formatFNA)                               // tracking fasta files
                {
                    if (buffer[startOfNextLine] == '>')
                        startOfFinalRead = startOfNextLine;
                }

                int lengthOfCurrentLine = bi - startOfNextLine + 1;                 // line length (including the trailing /n)
                lineLengthsInBuffer[li] = lengthOfCurrentLine;                      // remember the start of line index and move on
                startOfNextLine = bi + 1;                                           // move to the start of the next line
                li++;

                // running out of line-start slots so resize this array
                if (li == lineLengthsInBuffer.Length)
                {
                    Array.Resize<int>(ref lineLengths[cbi], lineLengthsInBuffer.Length + lineLengthsInBuffer.Length / 2);
                    lineLengthsInBuffer = lineLengths[cbi];                         // resize is tricky so need to refresh the local copy of the descriptor
                }
            }

            // find the incomplete read at the end of the buffer and work out where it ends
            if (fileFormat == SeqFiles.formatFASTQ)
            {
                if (fastqLineCount == 4 && lastLineComplete)                        // was last read complete? (4 lines)
                    endOfLastCompleteRead = charsInBuffer;                          // yes - no incomplete read to copy to start of next buffer
                else
                    endOfLastCompleteRead = startOfFinalRead;                       // the end of the line before that start of this read
            }

            if (fileFormat == SeqFiles.formatFNA)
            {
                if (EOF)
                    endOfLastCompleteRead = charsInBuffer;                          // EOF - last read is complete
                else
                    endOfLastCompleteRead = startOfFinalRead;                       // --> '>' of final read   
            }

            return endOfLastCompleteRead;
        }
    }

    // =========================================
    // ------------  BufferedWriter ------------
    // =========================================

    public class BufferedWriter
    {
        public class WriteBuffer
        {
            internal int bufferNo;
            internal int capacity;
            internal int currentIdx;
            internal char[] charBuffer;
            internal char[] newline;

            internal StreamWriter file;
            internal int setIdx;
            internal int setSize;
            internal string tag;

            internal List<WriteBuffer> availableBuffers;
            internal AutoResetEvent buffersAvailable;
        }

        StreamWriter seqFile;
        char[] newline;

        List<WriteBuffer> bufferPool;
        List<WriteBuffer> availableBuffers;
        AutoResetEvent buffersAvailable;

        static Queue<WriteBuffer> queuedBuffers;
        static AutoResetEvent buffersQueued;
        static List<WriteBuffer> buffersToWrite;

        //static List<string> bufferTrace;
        static int BWNo = 0;

        static Thread bufferWritingThread;

        // class constructor
        static BufferedWriter()
        {
            //bufferTrace = new List<string>(1000);

            // set of buffers to write atomically
            buffersToWrite = new List<WriteBuffer>();
            // signal that buffers are queued for writing
            buffersQueued = new AutoResetEvent(false);
            // buffers waiting to be written
            queuedBuffers = new Queue<WriteBuffer>();

            bufferWritingThread = new Thread(BufferWriter);
            bufferWritingThread.Priority = ThreadPriority.AboveNormal;
            bufferWritingThread.Name = "BufferedWriter";
            bufferWritingThread.Start();
        }

        // instance constructor
        public BufferedWriter(StreamWriter seqFile, int initialBufferSize, int noCallers, string tag)
        {
            BWNo++;
            this.seqFile = seqFile;
            this.newline = seqFile.NewLine.ToCharArray();

            int buffersWanted = noCallers * 2;

            // per-BufferWriter set of buffers
            bufferPool = new List<WriteBuffer>(buffersWanted);
            // free buffers available for use
            availableBuffers = new List<WriteBuffer>(buffersWanted);
            // signal that buffers are available
            buffersAvailable = new AutoResetEvent(false);

            for (int i = 0; i < buffersWanted; i++)
            {
                WriteBuffer buffer = new WriteBuffer();
                buffer.capacity = initialBufferSize;
                buffer.currentIdx = 0;
                buffer.charBuffer = new char[initialBufferSize];
                buffer.newline = seqFile.NewLine.ToCharArray();
                buffer.file = seqFile;
                buffer.tag = BWNo + "-" + tag;
                buffer.setIdx = 0;
                buffer.setSize = 0;

                buffer.bufferNo = bufferPool.Count;
                buffer.buffersAvailable = buffersAvailable;
                buffer.availableBuffers = availableBuffers;

                bufferPool.Add(buffer);
                availableBuffers.Add(buffer);
            }

            //lock (bufferTrace)
            //{
            //    bufferTrace.Add("added " + buffersWanted + " buffers to pool for " + tag);
            //}

            buffersAvailable.Set();
        }

        public void CloseBufferedWriter()
        {
            // waits for all buffers to be returned to the available pool by the writing thread
            while (true)
            {
                buffersAvailable.WaitOne();
                if (availableBuffers.Count == bufferPool.Count)
                    break;
            }
        }

        public WriteBuffer GetWriteBuffer()
        {
            WriteBuffer buffer = FindFreeBuffer();

            buffer.currentIdx = 0;
            buffer.setIdx = 0;
            buffer.setSize = 0;

            return buffer;
        }

        private WriteBuffer FindFreeBuffer()
        {
            WriteBuffer buffer = null;

            while (buffer == null)
            {
                lock (availableBuffers)
                {
                    if (availableBuffers.Count > 0)
                    {
                        buffer = availableBuffers[0];
                        availableBuffers.RemoveAt(0);
                        //lock (bufferTrace)
                        //{
                        //    bufferTrace.Add("allocated buffer " + buffer.bufferNo + " for " + buffer.tag);
                        //}
                    }
                }

                if (buffer == null)
                {
                    //lock (bufferTrace)
                    //{
                    //    bufferTrace.Add("waiting for buffer");
                    //}
                    buffersAvailable.WaitOne();
                }
            }

            return buffer;
        }

        public void WriteLine(WriteBuffer buffer, Sequence line, int start, int length)
        {
            int spaceNeeded = length + buffer.newline.Length * 2;
            int spaceAvailable = buffer.capacity - buffer.currentIdx;
            if (spaceNeeded > spaceAvailable)
            {
                buffer.capacity = buffer.capacity + (spaceNeeded * 100);
                Array.Resize<char>(ref buffer.charBuffer, buffer.capacity);
            }

            Array.Copy(line.Bases, start, buffer.charBuffer, buffer.currentIdx, length);
            buffer.currentIdx += length;
            Array.Copy(buffer.newline, 0, buffer.charBuffer, buffer.currentIdx, buffer.newline.Length);
            buffer.currentIdx += buffer.newline.Length;
        }

        public void WriteLine(WriteBuffer buffer, Sequence line)
        {
            WriteLine(buffer, line, 0, line.Length);
        }

        public void WriteLine(WriteBuffer buffer)
        {
            EnsureSpaceAvailable(ref buffer, buffer.newline.Length);
            Array.Copy(buffer.newline, 0, buffer.charBuffer, buffer.currentIdx, buffer.newline.Length);
            buffer.currentIdx += buffer.newline.Length;
        }

        public void Write(WriteBuffer buffer, string s)
        {
            EnsureSpaceAvailable(ref buffer, s.Length);
            Array.Copy(s.ToCharArray(), 0, buffer.charBuffer, buffer.currentIdx, s.Length);
            buffer.currentIdx += s.Length;
        }

        public void WriteLine(WriteBuffer buffer, string s)
        {
            EnsureSpaceAvailable(ref buffer, s.Length + buffer.newline.Length);
            Array.Copy(s.ToCharArray(), 0, buffer.charBuffer, buffer.currentIdx, s.Length);
            buffer.currentIdx += s.Length;
            Array.Copy(buffer.newline, 0, buffer.charBuffer, buffer.currentIdx, buffer.newline.Length);
            buffer.currentIdx += buffer.newline.Length;
        }

        private void EnsureSpaceAvailable(ref WriteBuffer buffer, int length)
        {
            int spaceNeeded = length + buffer.newline.Length * 2;
            int spaceAvailable = buffer.capacity - buffer.currentIdx;
            if (spaceNeeded > spaceAvailable)
            {
                buffer.capacity = buffer.capacity + (spaceNeeded * 100);
                Array.Resize<char>(ref buffer.charBuffer, buffer.capacity);
                //lock (bufferTrace)
                //{
                //    bufferTrace.Add("resized buffer " + buffer.bufferNo + " to " + buffer.capacity + " for " + buffer.tag);
                //}
            }
        }

        public static void WriteABuffer(WriteBuffer buffer)
        {
            if (buffer == null)
                return;

            // queue a single buffer for writing
            lock (queuedBuffers)
            {
                buffer.setIdx = 1;
                buffer.setSize = 1;

                queuedBuffers.Enqueue(buffer);
                //lock (bufferTrace)
                //{
                //    bufferTrace.Add("queued buffer " + buffer.bufferNo + " for " + buffer.tag);
                //}
            }

            buffersQueued.Set();
        }

        public static void WriteBufferSet(WriteBuffer[] readBuffers, WriteBuffer[] qualBuffers)
        {
            for (int i = 0; i < readBuffers.Length; i++)
            {
                int noOfBuffers = 1;
                if (qualBuffers != null && qualBuffers[i] != null)
                    noOfBuffers = 2;

                WriteBuffer readBuffer = readBuffers[i];
                WriteBuffer qualBuffer = null;

                if (readBuffer != null)
                {
                    readBuffer.setIdx = i + 1;
                    readBuffer.setSize = noOfBuffers;
                }

                if (qualBuffers != null && qualBuffers[i] != null)
                {
                    qualBuffer = qualBuffers[i];
                    qualBuffer.setIdx = i + 1;
                    qualBuffer.setSize = noOfBuffers;
                }

                lock (queuedBuffers)
                {
                    if (readBuffer != null)
                    {
                        queuedBuffers.Enqueue(readBuffer);
                        //lock (bufferTrace)
                        //{
                        //    bufferTrace.Add("queued read buffer " + readBuffer.bufferNo + " for " + readBuffer.tag);
                        //}
                    }
                    if (qualBuffer != null)
                    {
                        queuedBuffers.Enqueue(qualBuffer);
                        //lock (bufferTrace)
                        //{
                        //    bufferTrace.Add("queued qual buffer " + qualBuffer.bufferNo + " for " + qualBuffer.tag);
                        //}
                    }
                }
            }

            buffersQueued.Set();
        }

        public static void WriteBufferSet(WriteBuffer[] readBuffers)
        {
            WriteBufferSet(readBuffers, null);
        }

        public static void WritePairedBufferSet(WriteBuffer[] readBuffers)
        {
            WritePairedBufferSet(readBuffers, null);
        }

        public static void WritePairedBufferSet(WriteBuffer[] readBuffers, WriteBuffer[] qualBuffers)
        {
            // queue a set of (paired) read/quals buffers together

            int noOfBuffers = 0;
            for (int i = 0; i < readBuffers.Length; i++)
                if (readBuffers[i] != null)
                    noOfBuffers++;
            if (qualBuffers != null)
                for (int i = 0; i < qualBuffers.Length; i++)
                    if (qualBuffers[i] != null)
                        noOfBuffers++;

            int bufferIdx = 1;

            lock (queuedBuffers)
            {
                for (int i = 0; i < readBuffers.Length; i++)
                {
                    WriteBuffer buffer = readBuffers[i];
                    buffer.setIdx = bufferIdx;
                    buffer.setSize = noOfBuffers;
                    queuedBuffers.Enqueue(buffer);
                    //lock (bufferTrace)
                    //{
                    //    bufferTrace.Add("queued read buffer " + buffer.bufferNo + " for " + buffer.tag + " " + buffer.currentIdx + "/" + buffer.capacity + " bytes");
                    //}
                    bufferIdx++;
                }

                if (qualBuffers != null)
                    for (int i = 0; i < qualBuffers.Length; i++)
                    {
                        WriteBuffer qualBuffer = qualBuffers[i];
                        if (qualBuffer != null)
                        {
                            qualBuffer.setIdx = bufferIdx;
                            qualBuffer.setSize = noOfBuffers;
                            queuedBuffers.Enqueue(qualBuffer);
                            //lock (bufferTrace)
                            //{
                            //    bufferTrace.Add("queued qual buffer " + qualBuffer.bufferNo + " for " + qualBuffer.tag);
                            //}
                            bufferIdx++;
                        }
                    }

                buffersQueued.Set();
            }
        }

        public static void FinishBufferWriter()
        {
            lock (queuedBuffers)
            {
                // queue an end marker (null)
                queuedBuffers.Enqueue(null);
                //lock (bufferTrace)
                //{
                //    bufferTrace.Add("queued termination signal");
                //}
            }

            buffersQueued.Set();

            bufferWritingThread.Join();
        }

        private static void BufferWriter()
        {
            bool stopWritingThread = false;

            object fileSetLock = new object();

            while (!stopWritingThread)
            {
                // wait for something to appear in the queue
                buffersQueued.WaitOne();

                // always keep writing until the queue is empty
                while (queuedBuffers.Count > 0)
                {
                    // get the next queued buffer set (null --> finish thread)
                    lock (queuedBuffers)
                    {
                        if (queuedBuffers.Count > 0)
                        {
                            buffersToWrite.Clear();
                            WriteBuffer buffer = queuedBuffers.Dequeue();

                            if (buffer == null)
                            {
                                //lock (bufferTrace)
                                //{
                                //    bufferTrace.Add("terminating buffer writing thread");
                                //}
                                stopWritingThread = true;
                                break;
                            }

                            buffersToWrite.Add(buffer);

                            //lock (bufferTrace)
                            //{
                            //    bufferTrace.Add("removed buffer " + buffer.bufferNo + " from write queue" + " for " + buffer.tag);
                            //}

                            int currentBuffer = buffer.setIdx;
                            int bufferSetSize = buffer.setSize;
                            while (currentBuffer < bufferSetSize)
                            {
                                // buffers are always added as sets under lock protection so we can always read a complete set
                                buffer = queuedBuffers.Dequeue();
                                buffersToWrite.Add(buffer);
                                //lock (bufferTrace)
                                //{
                                //    bufferTrace.Add("removed buffer " + buffer.bufferNo + " from write queue" + " for " + buffer.tag);
                                //}
                                currentBuffer = buffer.setIdx;
                            }
                        }
                    }

                    // write out a set of buffers to their files (under lock protection to maintain pairing)
                    lock (fileSetLock)
                    {
                        for (int i = 0; i < buffersToWrite[0].setSize; i++)
                        {
                            WriteBuffer buffer = buffersToWrite[i];
                            buffer.file.Write(buffer.charBuffer, 0, buffer.currentIdx);
                            //lock (bufferTrace)
                            //{
                            //    bufferTrace.Add("wrote " + buffer.currentIdx + " bytes from buffer " + buffer.bufferNo + " for " + buffer.tag);
                            //}
                        }
                    }

                    // and return the now-written buffers
                    for (int i = 0; i < buffersToWrite.Count; i++)
                    {
                        lock (buffersToWrite[i].availableBuffers)
                        {
                            buffersToWrite[i].availableBuffers.Add(buffersToWrite[i]);
                            buffersToWrite[i].buffersAvailable.Set();
                            //lock (bufferTrace)
                            //{
                            //    bufferTrace.Add("returned buffer " + buffersToWrite[i].bufferNo + " to available list" + " for " + buffersToWrite[i].tag);
                            //}
                        }
                    }   
                                    
                } // write buffers from queue

            } // until told to stop

        } // buffer writing (thread)

    }
}
