using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.IO.Compression;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;

namespace WorkingDogsCore
{
    public class SeqFiles
    {
        public const int formatNone = 0;            // no format specified (ask Readers to determine it if possible)
        public const int formatFASTA = 1;			// fasta with multiple data lines (.fna, .fa, .fasta, .fas, ...)
        public const int formatFNA = 1;             // synonym (bacwards compatiblity with very old code)
        public const int formatSFA = 2;             // fasta with single data lines (assumed to fit inside a single buffer)
        public const int formatFASTQ = 3;           // fastq 

        enum Formats
        {
            None = 0,
            FASTA = 1,
            SFA = 2,
            FASTQ = 3
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

            if (readFormat == formatFNA || readFormat == formatSFA)
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

            if (readFormat == formatFNA || readFormat == formatSFA)
            {
                try
                {
                    read = ReadFASTA(readsFile, false, out readHeader);
                }
                catch (Exception e)
                {
                    Console.WriteLine("exception from ReadFASTA: " + e.Message);
                    return null;
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
            WriteRead(readsFile, readHeader, read, readFormat, 60);
        }

        public static void WriteRead(StreamWriter readsFile, string readHeader, string read, int readFormat, int lineLength)
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
                    int wLen = Math.Min(lineLength, hrLen - m);
                    readsFile.WriteLine(read.Substring(m, wLen));
                    m += lineLength;
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

        public static StreamReader OpenSeqStream(string seqFN)
        {
            FileStream FS = new FileStream(seqFN, FileMode.Open, FileAccess.Read, FileShare.Read, 262144);
            StreamReader RS;
            if (seqFN.EndsWith(".gz"))
            {
                GZipStream decompressor = new GZipStream(FS, CompressionMode.Decompress);
                RS = new StreamReader(decompressor, Encoding.ASCII, false);
            }
            else
                RS = new StreamReader(FS, Encoding.ASCII, false);

            return RS;
        }

        // The qual scores used with Fastq reads are ambiguous. Illumina data may use a base of 64 ('@') while Sanger formatted data
        // uses a base of 33 ('!'). This method reads quals from a reads file until it finds a value that can only be Sanger or until it's
        // looked at the first 100 reads, and then returns the appropriate qual offset.
        public static int ResolveFastqQualAmbiguity(string fastqReadsFN, out bool fullQualHeader)
        {
            fullQualHeader = true;
            StreamReader fastqReads = OpenSeqStream(fastqReadsFN);
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

        // find the best qual score in the first few reads
        public static int GetBestQual(string fastqReadsFN, int qualOffset)
        {
            StreamReader fastqReads = OpenSeqStream(fastqReadsFN);
            const int readsToExamine = 10;              // only look at the first 100 reads
            int readsExamined = 0;
            int bestQualScore = 0;

            while (readsExamined < readsToExamine)
            {
                string readHeader = fastqReads.ReadLine();
                if (readHeader == null)
                    break;
                
                readsExamined++;
                fastqReads.ReadLine();
                fastqReads.ReadLine();
                string quals = fastqReads.ReadLine();

                for (int i = 0; i < quals.Length; i++)
                {
                    if ((int)quals[i] - qualOffset > bestQualScore)
                        bestQualScore = (int)quals[i] - qualOffset;
                }
            }
            fastqReads.Close();

            return bestQualScore;
        }

        // finds any trailing poor qual region in a Sequence read and trims the read back to good bases
        public static int TrimTrailingPoorQuals(Sequence read, Sequence quals, int trimQual, int qualOffset)
        {
            const int windowLength = 10;
            const int passesNeededInWindow = windowLength * 3 / 4;
            int basesTrimmed;
            char[] qualCharArray = quals.Bases;
            int goodBaseIdx = 0;
            int lastQualIdx = quals.Length - 1;

            for (int i = lastQualIdx; i >= 0; i--)
                if (((int)qualCharArray[i] - qualOffset) > trimQual)
                {
                    goodBaseIdx = i;
                    break;
                }

            basesTrimmed = lastQualIdx - goodBaseIdx;

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

            //int[] qualScores = new int[quals.Length];
            //for (int i = 0; i < quals.Length; i++)
            //    qualScores[i] = (int)quals[i] - qualOffset;

            for (int i = quals.Length - 1; i >= 0; i--)
                if (((int)quals[i] - qualOffset) > trimQual)
                {
                    goodBaseIdx = i;
                    break;
                }

            basesToTrim = quals.Length - (goodBaseIdx + 1);

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

        // counts the total number of 'poor' qual scores in a read
        public static int CountPoorQuals(string quals, int minQual, int qualOffset)
        {
            int poorQualCount = 0;

            for (int i = 0; i < quals.Length; i++)
                if ((int)quals[i] - qualOffset < minQual)
                    poorQualCount++;

            return poorQualCount;
        }

        // sums the total qual scores for a read - higher numbers --> better quality read
        public static int SumQualScores(string quals, int qualOffset)
        {
            int summedQualCount = 0;

            for (int i = 0; i < quals.Length; i++)
                summedQualCount += (int)quals[i] - qualOffset;

            return summedQualCount;
        }

        // calculates a weighted poor qual score for a read - lower numbers --> better quality read
        // sum is weighted to reflect the cost of errors outside the overlap region in a read pair
        public static int SumWeightedPoorQualScore(string quals, int qualOffset, int bestQualScore, int overlapLength, int weighting)
        {
            int poorQualSum= 0;

            // part of read likely to be not covered by the paired read
            // errors here cannot be corrected from the paired read
            for (int i = 0; i < quals.Length - overlapLength; i++)
                poorQualSum += bestQualScore - ((int)quals[i] - qualOffset);
            poorQualSum = poorQualSum * weighting;

            // part of read likely to be covered by the paired read
            // errors here may be corrected from the paired read
            for (int i = overlapLength; i < quals.Length; i++)
                poorQualSum += bestQualScore - (int)quals[i] - qualOffset;

            return poorQualSum;
        }

        public static string LFConvention(string readsFN)
        {
            StreamReader reads = OpenSeqStream(readsFN);
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
            StreamReader file = OpenSeqStream(fn);

            string[] first4Lines = new string[4];
            for (int i = 0; i < first4Lines.Length; i++)
            {
                first4Lines[i] = file.ReadLine();
                if (first4Lines[i] == null)
                    break;
            }

            if (first4Lines[0][0] == '@' && first4Lines[2][0] == '+')   // normal fastq
                fileFormat = SeqFiles.formatFASTQ;
            if (first4Lines[0][0] == '@' && first4Lines[2][0] == '@')   // old Kelpie temp files (header+seq from fastq)
                fileFormat = SeqFiles.formatSFA;
            if (first4Lines[0][0] == '>' && first4Lines[2][0] == '>')   // appears to be single-line fasta (can always force .FASTA in caller)
                fileFormat = SeqFiles.formatSFA;
            if (first4Lines[0][0] == '>' && first4Lines[2][0] != '>')   // looks like multi-line fasta
                fileFormat = SeqFiles.formatFASTA;

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
    //
    // BufferedReader is NOT thread-safe. Most uses will be calling BR under lock protection - usually to ensure read-pairing is maintained, making the locks etc 
    // needed to make BR calls thread-safe somewhat redundant.

    public class BufferedReader
    {
        const int bufferLength = 1000000;                   // length of each buffer

        internal class ReadBuffer
        {
            internal char[] buffer;                         // the set of buffers. Each buffer will start at the start of a read.
            internal int[] lineLengths;                     // lengths of the lines found in these buffers (including \n)
            internal int endOfLastCompleteRead = -1;        // index of the end of the last complete read in each buffer
            internal int startOfBuffer = 0;                 // where reads into this buffer should start (after left overs from previous read have been copied in)
            internal int nextLineIdx = 0;                   // next lineLength index to be used 
            internal int startOfNextRead = 0;               // index into buffer of start of next read

            public ReadBuffer()
            {
                buffer = new char[bufferLength + 1];        // in case we need to insert a \n at the end of the final buffer
                lineLengths = new int[bufferLength / 20];   // 5% (bytes --> lines)
                endOfLastCompleteRead = -1;
                startOfBuffer = 0;
            }
        }

        int fileFormat = SeqFiles.formatNone;               // fasta or fastq
        StreamReader readsFile = null;                      // fastq: seqs+quals; fasta: seqs only
        StreamReader qualsFile = null;                      // fasta only: separate quals file (rare and not buffered)

        bool readsFromBuffer;                               // buffered reads (usual) or direct calls (fasta with separate quals file)
        Task<bool> bufferFilling;                           // last async buffer filling task
        const int noOfBuffers = 2;                          // initial number of buffers - both in-use and available
        List<ReadBuffer> buffers = new List<ReadBuffer>(noOfBuffers);   // the set of buffers. Each buffer will start at the start of a read.
        List<int> freeBuffers = new List<int>(noOfBuffers);             // available buffers waiting to be filled - indexes into buffers
        List<int> filledBuffers = new List<int>(noOfBuffers);           // filled buffers, waiting to be read - in file order

        int currentBufferIdx;                               // index to current buffer
        int nextBufferIdx;                                  // index to next buffer to use (holds unused end of current buffer)
        bool EOF = false;                                   // EOF reached on readsFile - stop calling FillBuffer
        bool bufferEmpty = true;                            // previous ReadRead call returned the last full read in the buffer - time to swap buffers 

        bool eolKnown = false;                              // has EOL convention been set yet?
        int eolLength;                                      // length of end-of-line (CR, LF or CR-LF). 
        char eolChar;                                       // char at end of each line (LF & CR-LF --> LF, CR --> CR)
        int linesPerRead;                                   // #lines/read for fastq (4) & single-line fasta (sfa)(2)
        char headerChar;                                    // headers start with this char

        public BufferedReader(int fileFormat, StreamReader readsFile, StreamReader qualsFile)
        {
            this.fileFormat = fileFormat;
            this.readsFile = readsFile;
            this.qualsFile = qualsFile;

            // only do buffering for FASTQ files and known single-line FASTA (such as Kelpie temp files)
            // FASTA files tend to be small and reads could cross buffer boundaries if they are multi-line.
            readsFromBuffer = fileFormat == SeqFiles.formatFASTQ || fileFormat == SeqFiles.formatSFA;
            if (fileFormat == SeqFiles.formatFASTQ)
                linesPerRead = 4;
            if (fileFormat == SeqFiles.formatSFA)
                linesPerRead = 2;
            // multi-line fasta is not buffered, so linesPerRead is irrelevant

            // unresolved file type - tell caller and set EOF
            if (fileFormat == SeqFiles.formatNone)
            {
                Console.WriteLine("file type could not be determined - returning EOF");
                EOF = true;
            }

            headerChar = (char)readsFile.Peek();

            if (readsFromBuffer)
            {
                for (int i = 0; i < noOfBuffers; i++)
                {
                    buffers.Add(new ReadBuffer());          // allocate an empty buffer           
                    freeBuffers.Add(i);                     // all buffers are free initially
                }

                EOF = false;
                currentBufferIdx = -1;
                // initiate filling the initial buffer
                bufferFilling = Task.Run(() => FillBuffer(readsFile, GetFreeBuffer(), GetFreeBuffer()));
            }
        }

        public void Close()
        {
            if (!EOF)
                bufferFilling.Wait(); ;
            readsFile.Close();
            if (qualsFile != null)
                qualsFile.Close();
        }

        public int ReadReads(int readsWanted, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals/*, int client*/)
        {
            int readsIdx = 0;
            bool EOFReached = false;

            if (readsFromBuffer)
            {
                while (true)
                {
                    if (bufferEmpty)
                    {
                        if (EOF)
                            break;

                        // return previous buffer (if one exists)
                        if (currentBufferIdx >= 0)
                            freeBuffers.Add(currentBufferIdx);

                        // wait for in-process buffer fill to complete
                        //Console.WriteLine(client + " waiting for FillBuffer");
                        bufferFilling.Wait();

                        EOF = bufferFilling.Result;
                        currentBufferIdx = filledBuffers[0];
                        //Console.WriteLine(client + " FillBuffer returned (" + EOF + "," + currentBufferIdx + ")");
                        filledBuffers.RemoveAt(0);
                        buffers[currentBufferIdx].nextLineIdx = 0;
                        buffers[currentBufferIdx].startOfNextRead = 0;

                        // and start the next buffer fill
                        if (!EOFReached)
                        {
                            int nbi = nextBufferIdx;
                            int fbi = GetFreeBuffer();
                            bufferFilling = Task.Run(() => FillBuffer(readsFile, nbi, fbi));
                        }
                    }

                    bufferEmpty = GetNextRead(currentBufferIdx, readHeaders[readsIdx], readSeqs[readsIdx], qualHeaders[readsIdx], quals[readsIdx]);

                    readsIdx++;
                    if (readsIdx == readsWanted)
                        break;
                }
            }
            else
            {
                // muliti-line FASTA files - just do unbuffered reads - add lock in rare (obsolete) case of separate fasta & quals to make reading the pairs of files thread-safe
                if (qualsFile != null)
                    lock (readsFile)
                    {
                        readsIdx = SeqFiles.ReadReads(readsWanted, readsFile, qualsFile, fileFormat, readHeaders, readSeqs, qualHeaders, quals);
                    }
                else
                    readsIdx = SeqFiles.ReadReads(readsWanted, readsFile, qualsFile, fileFormat, readHeaders, readSeqs, qualHeaders, quals);
                EOF = readsWanted != readsIdx;
            }

            return readsIdx;
        }

        public bool ReadRead(Sequence readHeader, Sequence readSeq, Sequence qualHeader, Sequence quals)
        {
            if (readsFromBuffer)
            {
                if (bufferEmpty && !EOF)
                {
                    if (currentBufferIdx >= 0)
                        freeBuffers.Add(currentBufferIdx);

                    // wait for previous buffer fill to complete
                    bufferFilling.Wait();
                    EOF = bufferFilling.Result;

                    currentBufferIdx = filledBuffers[0];
                    filledBuffers.RemoveAt(0);
                    buffers[currentBufferIdx].nextLineIdx = 0;
                    buffers[currentBufferIdx].startOfNextRead = 0;

                    // and start the next buffer read
                    if (!EOF)
                    {
                        int nbi = nextBufferIdx;
                        int fbi = GetFreeBuffer();
                        bufferFilling = Task.Run(() => FillBuffer(readsFile, nbi, fbi));
                    }
                }

                bufferEmpty = GetNextRead(currentBufferIdx, readHeader, readSeq, qualHeader, quals);
            }
            else
            {
                // general (possibly multi-line) FASTA - just use unbuffered reads - adding lock if we have two files to read (seq + quals)
                if (qualsFile != null)
                    lock (readsFile)
                    {
                        EOF = SeqFiles.ReadRead(readsFile, qualsFile, fileFormat, readHeader, readSeq, qualHeader, quals);
                    }
                else
                    EOF = SeqFiles.ReadRead(readsFile, qualsFile, fileFormat, readHeader, readSeq, qualHeader, quals);
            }

            return EOF;
        }

        private int GetFreeBuffer()
        {
            int bufferToUse;

            if (freeBuffers.Count == 0)
            {
                // always return a buffer
                buffers.Add(new ReadBuffer());
                freeBuffers.Add(buffers.Count - 1);
                //Console.WriteLine("added buffer #" + buffers.Count);
            }

            bufferToUse = freeBuffers[0];
            freeBuffers.RemoveAt(0);

            return bufferToUse;
        }

        private bool GetNextRead(int cbi, Sequence readHeader, Sequence readSeq, Sequence qualHeader, Sequence quals)
        {
            char[] currentBuffer = buffers[cbi].buffer;
            int[] currentLineLengths = buffers[cbi].lineLengths;
            int snr = buffers[cbi].startOfNextRead;
            int nli = buffers[cbi].nextLineIdx;
            int currentEndOfLastCompleteRead = buffers[cbi].endOfLastCompleteRead;

            CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readHeader, 0);
            snr += currentLineLengths[nli];
            nli++;
            CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readSeq, 0);
            snr += currentLineLengths[nli];
            nli++;
            if (fileFormat == SeqFiles.formatFASTQ)
            {
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], qualHeader, 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], quals, 0);
                snr += currentLineLengths[nli];
                nli++;
            }
            else
            {
                qualHeader.Length = 0;
                quals.Length = 0;
            }

            buffers[cbi].nextLineIdx = nli;
            buffers[cbi].startOfNextRead = snr;
            return snr == currentEndOfLastCompleteRead || currentLineLengths[nli] == 0;
        }

        private void CopySeqFromBuffer(char[] buffer, int startOfSeq, int length, Sequence destination, int startIdx)
        {
            // copy up to but not including the \n
            int lengthAfterCopy = length + startIdx - eolLength;
            if (lengthAfterCopy > destination.Capacity)
                destination.Resize(lengthAfterCopy + lengthAfterCopy / 2);
            Array.Copy(buffer, startOfSeq, destination.Bases, startIdx, lengthAfterCopy);
            destination.Length = lengthAfterCopy;
        }

        private bool FillBuffer(StreamReader reads, int cbi, int nbi)
        {
            //Console.WriteLine("filling buffer #" + cbi);

            char[] currentBuffer = buffers[cbi].buffer;
            char[] nextBuffer = buffers[nbi].buffer;

            int startIdx = buffers[cbi].startOfBuffer;
            int charsToRead = bufferLength - startIdx;
            int charsRead = reads.ReadBlock(currentBuffer, startIdx, charsToRead);
            bool EOFReached = charsRead != charsToRead;
            int charsInBuffer = startIdx + charsRead;

            // find the line boundaries to save having to scan the buffer again when extracting actual reads
            // and also finds the end of the incomplete read at the end of the buffer
            buffers[cbi].endOfLastCompleteRead = FindLinesInBuffer(cbi, charsInBuffer, EOF);

            filledBuffers.Add(cbi);
            //Console.WriteLine(filledBuffers.Count + " filled buffers");

            // copy any remaining chars into the next buffer
            int startOfLeftover = buffers[cbi].endOfLastCompleteRead;
            if (startOfLeftover > 1)
            {
                int leftOverLength = bufferLength - startOfLeftover;
                Array.Copy(currentBuffer, startOfLeftover, nextBuffer, 0, leftOverLength);
                buffers[nbi].startOfBuffer = leftOverLength;
            }
            nextBufferIdx = nbi;

            return EOFReached;
        }

        private int FindLinesInBuffer(int cbi, int charsInBuffer, bool EOF)
        {
            // finds all the start-of-lines (char after /n) in the buffer and remembers them in startsOfLines. These are the indices in the current
            // buffer of the first char of every line. The last entry in startsOfLines will be -1.
            // It then goes through and finds the start/end of all the reads (format-dependent).
            // The start of the last incomplete read is passed back so this region can be copied into the next buffer when it's filled.

            char[] buffer = buffers[cbi].buffer;
            int[] lineLengthsInBuffer = buffers[cbi].lineLengths;

            int li = 0;                     // start from the first line in the buffer
            int startOfNextLine = 0;        // start of the next line (starting point of scan for next \n)
            int startOfFinalRead = 0;       // char index of start of last read
            int endOfLastCompleteRead;      // index of end of the last complete read in the buffer (return value)
            int readLineCount = 0;          // count 2 0r 4 lines for each read
            bool lastLineComplete = true;   // \n at very end of buffer

            // find out what the EOL convention is in the file (if not already known)
            if (!eolKnown)
            {
                for (int i = 0; i < buffer.Length; i++)
                {
                    if (buffer[i] == '\n' || buffer[i] == '\r')
                    {
                        eolKnown = true;

                        // CR-LF
                        if (buffer[i] == '\r' && buffer[i + 1] == '\n')
                        {
                            eolChar = '\n';
                            eolLength = 2;
                            break;
                        }

                        // CR or LF
                        eolChar = buffer[i];
                        eolLength = 1;
                        break;
                    }
                }
            }

            // force \n at end of buffer at EOF 
            if (EOF && buffer[charsInBuffer - 1] != eolChar)
            {
                buffer[charsInBuffer] = eolChar;
                charsInBuffer++;
            }

            // find and remember all the \n locations (plus 1 to get to the start of the next line - may be incomplete)
            while (startOfNextLine < charsInBuffer)
            {
                int bi = Array.IndexOf<char>(buffer, eolChar, startOfNextLine);     // find the next \n
                if (bi == -1)                                                       // no EOL at end of buffer
                {
                    lineLengthsInBuffer[li] = 0;                                    // incomplete line so don't record a length
                    lastLineComplete = false;
                    break;
                }

                // now have the next complete line - so track reads
                readLineCount++;                                                   // keep track of where we are in fastq quad
                if (readLineCount > linesPerRead)
                    readLineCount = 1;

                if (readLineCount == 1 && buffer[startOfNextLine] == headerChar)    // remember each header line as it goes past
                    startOfFinalRead = startOfNextLine;

                int lengthOfCurrentLine = bi - startOfNextLine + 1;                 // line length (including the trailing \n)
                lineLengthsInBuffer[li] = lengthOfCurrentLine;                      // remember the start of line index and move on
                startOfNextLine = bi + 1;                                           // move to the start of the next line
                li++;

                // running out of line-start slots so resize this array
                if (li == lineLengthsInBuffer.Length)
                {
                    Array.Resize<int>(ref lineLengthsInBuffer, lineLengthsInBuffer.Length + lineLengthsInBuffer.Length / 2);
                    buffers[cbi].lineLengths = lineLengthsInBuffer;                 // resize is tricky so need to refresh the local copy of the descriptor first

                }
            }

            // find the incomplete read at the end of the buffer and work out where it ends
            if (readLineCount == linesPerRead && lastLineComplete)                  // was last read complete? (4 lines)
                endOfLastCompleteRead = charsInBuffer;                              // yes - no incomplete read to copy to start of next buffer
            else
                endOfLastCompleteRead = startOfFinalRead;                           // the end of the line before that start of this read

            return endOfLastCompleteRead;
        }
    }

#if NET6_0_OR_GREATER
    // ============================================
    // ----------- BufferedBlockReader ------------
    // ============================================

    // Getting blocks from a BufferedBlockReader is thread-safe (but need to be lock-protected for paired reads). Calls to get reads from within a buffer block are thread-safe.
    // Only fastq and single-line fatsa are supported.
    // Client threads get a buffer block for their exclusive use by calling GetBufferBlock. They can then get reads from this buffer block via Read... calls

    public class BufferedBlockReader
    {
        const int bufferLength = 1000000;                   // length of each buffer

        internal class ReadBuffer
        {
            internal char[] buffer;                         // the set of buffers. Each buffer will start at the start of a read.
            internal int[] lineLengths;                     // lengths of the lines found in these buffers (including \n)
            internal int endOfLastCompleteRead = -1;        // index of the end of the last complete read in each buffer
            internal int startOfBuffer = 0;                 // where reads into this buffer should start (after left overs from previous read have been copied in)
            internal int nextLineIdx = 0;                   // next lineLength index to be used 
            internal int startOfNextRead = 0;               // index into buffer of start of next read

            public ReadBuffer()
            {
                buffer = new char[bufferLength + 2];        // in case we need to insert a \n at the end of the final buffer
                lineLengths = new int[bufferLength / 20];   // 5% (bytes --> lines)
                endOfLastCompleteRead = -1;
                startOfBuffer = 0;
            }
        }

        int fileFormat = SeqFiles.formatNone;               // fasta or fastq (no support for fasta with quals)
        StreamReader readsFile = null;                      // fastq: seqs+quals; fasta: single-line seqs only

        Task<bool> bufferFilling;                           // last async buffer filling task
        int noOfBuffers;                                    // initial number of buffers - both in-use and available
        List<ReadBuffer> buffers = new List<ReadBuffer>();  // the set of buffers. Each buffer will start at the start of a read.
        List<int> freeBuffers = new List<int>();            // available buffers waiting to be filled - indexes into buffers
        List<int> filledBuffers = new List<int>();          // filled buffers, waiting to be read - in file order
        private AutoResetEvent buffersAvailable;            // event used to ensure single-threaded access to read-next-buffer code (file needs to be read sequentially)

        bool EOF;                                           // EOF reached on readsFile - stop calling FillBuffer
        int nextBufferIdx = -1;                             // partially filled buffer - left overs from previous fill
        bool eolKnown = false;                              // has EOL convention been set yet?
        int eolLength;                                      // length of end-of-line (CR, LF or CR-LF). 
        char eolChar;                                       // char at end of each line (LF & CR-LF --> LF, CR --> CR)
        int linesPerRead;                                   // #lines/read for fastq (4) & single-line fasta (sfa)(2)
        char headerChar;                                    // headers start with this char

        public BufferedBlockReader(int fileFormat, StreamReader readsFile, int noOfCallers)
        {
            this.fileFormat = fileFormat;
            this.readsFile = readsFile;
            noOfBuffers = noOfCallers + 2;
            EOF = false;

            headerChar = (char)readsFile.Peek();    
            
            if (fileFormat == SeqFiles.formatFASTQ)
                linesPerRead = 4;
            if (fileFormat == SeqFiles.formatSFA)
                linesPerRead = 2;

            // BufferedBlockReader does not support multi-line fasta - so write out a message and force EOF 
            if (this.fileFormat == SeqFiles.formatFASTA)
            {
                Console.WriteLine("BufferedBlockReader does not support multi-line fasta files");
                EOF = true;
            }
            // file format couldn't be determined
            if (this.fileFormat == SeqFiles.formatFASTA)
            {
                Console.WriteLine("File format could not be determined");
                EOF = true;
            }

            for (int i = 0; i < noOfBuffers; i++)
            {
                buffers.Add(new ReadBuffer());              // allocate an empty buffer           
                freeBuffers.Add(i);                         // all buffers are free initially
            }


            buffersAvailable = new AutoResetEvent(true);    // allow the first caller to wait for FillBuffer to complete
            // initiate filling the initial buffer
            bufferFilling = Task.Run(() => FillBuffer(readsFile, GetFreeBuffer(), GetFreeBuffer()));
        }

        public void Close()
        {
            if (!EOF)
                bufferFilling.Wait(); 
            readsFile.Close();
        }

        public int GetBufferBlock(/*int client*/)
        {
            int bufferIdx = -1;

            // gate access into the Wait code
            //Console.WriteLine("client " + client + " waiting on GetBufferBlock");
            buffersAvailable.WaitOne();

            if (EOF)
            {
                //Console.WriteLine("GetBufferBlock found EOF after waiting for client " + client);
                buffersAvailable.Set();         // propagate event to any other waiters
                return -1;
            }

            // wait for any previous buffer fill to complete
            //Console.WriteLine("client " + client + " waiting on FillBuffer");
            bufferFilling.Wait();
            EOF = bufferFilling.Result;

            //if (EOF)
            //    Console.WriteLine("FillBuffer returned EOF for client " + client);

            if (filledBuffers.Count > 0)
            {
                bufferIdx = filledBuffers[0];
                filledBuffers.RemoveAt(0);
            }

            buffers[bufferIdx].nextLineIdx = 0;
            buffers[bufferIdx].startOfNextRead = 0;

            // and start the next buffer read
            if (!EOF)
            {
                int nbi = nextBufferIdx;
                int fbi = GetFreeBuffer();
                bufferFilling = Task.Run(() => FillBuffer(readsFile, nbi, fbi));
            }

            buffersAvailable.Set();

            //Console.WriteLine("giving buffer #" + bufferIdx + " to client #" + client);
            return bufferIdx;
        }

        public bool ReadReadFromBlock(int bufferIdx, out Span<char> readHeader, out Span<char> readSeq, out Span<char> qualHeader, out Span<char> quals)
        {
            if (buffers[bufferIdx].startOfNextRead == buffers[bufferIdx].endOfLastCompleteRead)
            {
                //Console.WriteLine("returned buffer #" + bufferIdx + " to free buffers");
                ReturnBuffer(bufferIdx);
                readHeader = Span<char>.Empty;
                readSeq = Span<char>.Empty;
                qualHeader = Span<char>.Empty;
                quals = Span<char>.Empty;
                return false;
            }

            char[] currentBuffer = buffers[bufferIdx].buffer;
            int[] currentLineLengths = buffers[bufferIdx].lineLengths;
            int snr = buffers[bufferIdx].startOfNextRead;
            int nli = buffers[bufferIdx].nextLineIdx;

            readHeader = new Span<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
            snr += currentLineLengths[nli];
            nli++;
            readSeq = new Span<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
            snr += currentLineLengths[nli];
            nli++;
            if (fileFormat == SeqFiles.formatFASTQ)
            {
                qualHeader = new Span<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
                quals = new Span<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
            }
            else
            {
                qualHeader = Span<char>.Empty;
                quals = Span<char>.Empty;
            }

            buffers[bufferIdx].nextLineIdx = nli;
            buffers[bufferIdx].startOfNextRead = snr;

            // return true when buffer is exhausted (client will need to call GetBufferBlock to get next buffer)
            return true;
        }

        public bool ReadReadFromBlock(int bufferIdx, out Memory<char> readHeader, out Memory<char> readSeq, out Memory<char> qualHeader, out Memory<char> quals)
        {
            if (buffers[bufferIdx].startOfNextRead == buffers[bufferIdx].endOfLastCompleteRead)
            {
                //Console.WriteLine("returned buffer #" + bufferIdx + " to free buffers");
                ReturnBuffer(bufferIdx);
                readHeader = Memory<char>.Empty;
                readSeq = Memory<char>.Empty;
                qualHeader = Memory<char>.Empty;
                quals = Memory<char>.Empty;
                return false;
            }

            char[] currentBuffer = buffers[bufferIdx].buffer;
            int[] currentLineLengths = buffers[bufferIdx].lineLengths;
            int snr = buffers[bufferIdx].startOfNextRead;
            int nli = buffers[bufferIdx].nextLineIdx;

            readHeader = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
            snr += currentLineLengths[nli];
            nli++;
            readSeq = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
            snr += currentLineLengths[nli];
            nli++;
            if (fileFormat == SeqFiles.formatFASTQ)
            {
                qualHeader = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
                quals = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
            }
            else
            {
                qualHeader = Memory<char>.Empty;
                quals = Memory<char>.Empty;
            }

            buffers[bufferIdx].nextLineIdx = nli;
            buffers[bufferIdx].startOfNextRead = snr;

            // return true when buffer is exhausted (client will need to call GetBufferBlock to get next buffer)
            return true;
        }

        public bool ReadReadFromBlock(int bufferIdx, Sequence readHeader, Sequence readSeq, Sequence qualHeader, Sequence quals)
        {
            if (buffers[bufferIdx].startOfNextRead == buffers[bufferIdx].endOfLastCompleteRead)
            {
                //Console.WriteLine("returned buffer #" + bufferIdx + " to free buffers");
                ReturnBuffer(bufferIdx);
                readHeader.Length = 0;
                readSeq.Length = 0; 
                qualHeader.Length = 0;
                quals.Length = 0;
                return false;
            }

            char[] currentBuffer = buffers[bufferIdx].buffer;
            int[] currentLineLengths = buffers[bufferIdx].lineLengths;
            int snr = buffers[bufferIdx].startOfNextRead;
            int nli = buffers[bufferIdx].nextLineIdx;

            CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readHeader, 0);
            snr += currentLineLengths[nli];
            nli++;
            CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readSeq, 0);
            snr += currentLineLengths[nli];
            nli++;
            if (fileFormat == SeqFiles.formatFASTQ)
            {
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], qualHeader, 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], quals, 0);
                snr += currentLineLengths[nli];
                nli++;
            }
            else
            {
                qualHeader.Length = 0;
                quals.Length = 0;
            }

            buffers[bufferIdx].nextLineIdx = nli;
            buffers[bufferIdx].startOfNextRead = snr;

            // return true when buffer is exhausted (client will need to call GetBufferBlock to get next buffer)
            return true;
        }

        public int ReadReadsFromBlock(int bufferIdx, int readsWanted, Memory<char>[] readHeaders, Memory<char>[] readSeqs, Memory<char>[] qualHeaders, Memory<char>[] quals)
        {
            int readsRead = 0;
            char[] currentBuffer = buffers[bufferIdx].buffer;
            int[] currentLineLengths = buffers[bufferIdx].lineLengths;
            int snr = buffers[bufferIdx].startOfNextRead;
            int nli = buffers[bufferIdx].nextLineIdx;

            while (readsRead < readsWanted)
            {
                if (snr == buffers[bufferIdx].endOfLastCompleteRead)
                    break;

                readHeaders[readsRead] = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
                readSeqs[readsRead] = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                snr += currentLineLengths[nli];
                nli++;
                if (fileFormat == SeqFiles.formatFASTQ)
                {
                    qualHeaders[readsRead] = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                    snr += currentLineLengths[nli];
                    nli++;
                    quals[readsRead] = new Memory<char>(currentBuffer, snr, currentLineLengths[nli] - eolLength);
                    snr += currentLineLengths[nli];
                    nli++;
                }
                else
                {
                    qualHeaders[readsRead] = Memory<char>.Empty;
                    quals[readsRead] = Memory<char>.Empty;
                }

                readsRead++;
            }

            buffers[bufferIdx].startOfNextRead = snr;
            buffers[bufferIdx].nextLineIdx = nli;
            return readsRead;
        }

        public int ReadReadsFromBlock(int bufferIdx, int readsWanted, Sequence[] readHeaders, Sequence[] readSeqs, Sequence[] qualHeaders, Sequence[] quals)
        {
            int readsRead = 0;
            char[] currentBuffer = buffers[bufferIdx].buffer;
            int[] currentLineLengths = buffers[bufferIdx].lineLengths;
            int snr = buffers[bufferIdx].startOfNextRead;
            int nli = buffers[bufferIdx].nextLineIdx;

            while (readsRead < readsWanted)
            {
                if (snr == buffers[bufferIdx].endOfLastCompleteRead)
                    break;

                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readHeaders[readsRead], 0);
                snr += currentLineLengths[nli];
                nli++;
                CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], readSeqs[readsRead], 0);
                snr += currentLineLengths[nli];
                nli++;
                if (fileFormat == SeqFiles.formatFASTQ)
                {
                    CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], qualHeaders[readsRead], 0);
                    snr += currentLineLengths[nli];
                    nli++;
                    CopySeqFromBuffer(currentBuffer, snr, currentLineLengths[nli], quals[readsRead], 0);
                    snr += currentLineLengths[nli];
                    nli++;
                }
                else
                {
                    qualHeaders[readsRead].Length = 0;
                    quals[readsRead].Length = 0;
                }

                readsRead++;
            }

            buffers[bufferIdx].startOfNextRead = snr;
            buffers[bufferIdx].nextLineIdx = nli;
            return readsRead;
        }

        private void CopySeqFromBuffer(char[] buffer, int startOfSeq, int length, Sequence destination, int startIdx)
        {
            // copy up to but not including the \n
            int lengthAfterCopy = length + startIdx - eolLength;
            if (lengthAfterCopy > destination.Capacity)
                destination.Resize(lengthAfterCopy + lengthAfterCopy / 2);
            Array.Copy(buffer, startOfSeq, destination.Bases, startIdx, lengthAfterCopy);
            destination.Length = lengthAfterCopy;
        }

        private int GetFreeBuffer()
        {
            int bufferToUse;

            lock (buffers)
            {
                if (freeBuffers.Count == 0)
                {
                    // always return a buffer
                    buffers.Add(new ReadBuffer());
                    freeBuffers.Add(buffers.Count - 1);
                    //Console.WriteLine("added new buffer #" + buffers.Count + " to free buffers");
                }

                bufferToUse = freeBuffers[0];
                freeBuffers.RemoveAt(0);
                //Console.WriteLine("assigned buffer #" + bufferToUse + " from free buffers");
            }

            return bufferToUse;
        }

        private void ReturnBuffer(int bufferIdx)
        {
            lock (buffers)
            {
                //Console.WriteLine("returning buffer #" + bufferIdx);
                freeBuffers.Add(bufferIdx);
            }
        }

        private bool FillBuffer(StreamReader reads, int cbi, int nbi)
        {
            //Console.WriteLine("filling buffer #" + cbi);

            char[] currentBuffer = buffers[cbi].buffer;
            char[] nextBuffer = buffers[nbi].buffer;

            int startIdx = buffers[cbi].startOfBuffer;
            int charsToRead = bufferLength - startIdx;
            int charsRead = reads.ReadBlock(currentBuffer, startIdx, charsToRead);
            bool EOFReached = charsRead != charsToRead;
            int charsInBuffer = startIdx + charsRead;

            // find the line boundaries to save having to scan the buffer again when extracting actual reads
            // and also finds the end of the incomplete read at the end of the buffer
            buffers[cbi].endOfLastCompleteRead = FindLinesInBuffer(cbi, charsInBuffer, EOF);

            //Console.WriteLine("added buffer #" + cbi + " to filled");
            lock (buffers)
            {
                filledBuffers.Add(cbi);
            }

            // copy any remaining chars into the next buffer
            int startOfLeftover = buffers[cbi].endOfLastCompleteRead;
            int leftOverLength = bufferLength - startOfLeftover;
            if (startOfLeftover > 1 && !EOFReached)
            {
                Array.Copy(currentBuffer, startOfLeftover, nextBuffer, 0, leftOverLength);
                buffers[nbi].startOfBuffer = leftOverLength;
                //Console.WriteLine("copied " + leftOverLength + " leftovers to buffer #" + nbi + " chars read=" + charsRead);      
            }
            nextBufferIdx = nbi;


            return EOFReached;
        }

        private int FindLinesInBuffer(int cbi, int charsInBuffer, bool EOF)
        {
            // finds all the start-of-lines (char after /n) in the buffer and remembers them in startsOfLines. These are the indices in the current
            // buffer of the first char of every line. The last entry in startsOfLines will be -1.
            // It then goes through and finds the start/end of all the reads (format-dependent).
            // The start of the last incomplete read is passed back so this region can be copied into the next buffer when it's filled.

            char[] buffer = buffers[cbi].buffer;
            int[] lineLengthsInBuffer = buffers[cbi].lineLengths;

            int li = 0;                     // start from the first line in the buffer
            int startOfNextLine = 0;        // start of the next line (starting point of scan for next \n)
            int startOfFinalRead = 0;       // char index of start of last read
            int endOfLastCompleteRead;      // index of end of the last complete read in the buffer (return value)
            int readLineCount = 0;          // count 2 0r 4 lines for each read
            bool lastLineComplete = true;   // \n at very end of buffer

            // find out what the EOL convention is in the file (if not already known)
            if (!eolKnown)
            {
                for (int i = 0; i < buffer.Length; i++)
                {
                    if (buffer[i] == '\n' || buffer[i] == '\r')
                    {
                        eolKnown = true;

                        // CR-LF
                        if (buffer[i] == '\r' && buffer[i + 1] == '\n')
                        {
                            eolChar = '\n';
                            eolLength = 2;
                            break;
                        }

                        // CR or LF
                        eolChar = buffer[i];
                        eolLength = 1;
                        break;
                    }
                }
            }

            // force \n at end of buffer at EOF 
            if (EOF && buffer[charsInBuffer - 1] != eolChar)
            {
                buffer[charsInBuffer] = eolChar;
                charsInBuffer++;
            }

            // find and remember all the \n locations (plus 1 to get to the start of the next line - may be incomplete)
            while (startOfNextLine < charsInBuffer)
            {
                int bi = Array.IndexOf<char>(buffer, eolChar, startOfNextLine);     // find the next \n
                if (bi == -1)                                                       // no EOL at end of buffer
                {
                    lineLengthsInBuffer[li] = 0;                                    // incomplete line so don't record a length
                    lastLineComplete = false;
                    break;
                }

                // now have the next complete line - so track reads
                readLineCount++;                                                   // keep track of where we are in fastq quad
                if (readLineCount > linesPerRead)
                    readLineCount = 1;

                if (readLineCount == 1 && buffer[startOfNextLine] == headerChar)    // remember each header line as it goes past
                    startOfFinalRead = startOfNextLine;

                int lengthOfCurrentLine = bi - startOfNextLine + 1;                 // line length (including the trailing \n)
                lineLengthsInBuffer[li] = lengthOfCurrentLine;                      // remember the start of line index and move on
                startOfNextLine = bi + 1;                                           // move to the start of the next line
                li++;

                // running out of line-start slots so resize this array
                if (li == lineLengthsInBuffer.Length)
                {
                    Array.Resize<int>(ref lineLengthsInBuffer, lineLengthsInBuffer.Length + lineLengthsInBuffer.Length / 2);
                    buffers[cbi].lineLengths = lineLengthsInBuffer;                 // resize is tricky so need to refresh the local copy of the descriptor first

                }
            }

            // find the incomplete read at the end of the buffer and work out where it ends
            if (readLineCount == linesPerRead && lastLineComplete)                  // was last read complete? (4 lines)
                endOfLastCompleteRead = charsInBuffer;                              // yes - no incomplete read to copy to start of next buffer
            else
                endOfLastCompleteRead = startOfFinalRead;                           // the end of the line before that start of this read

            return endOfLastCompleteRead;
        }
    }

#endif


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

        private List<WriteBuffer> bufferPool;
        private List<WriteBuffer> availableBuffers;
        private AutoResetEvent buffersAvailable;

        private static Queue<WriteBuffer> queuedBuffers;
        private static AutoResetEvent buffersQueued;
        private static List<WriteBuffer> buffersToWrite;

        //static List<string> bufferTrace;
        private static int BWNo = 0;

        private static Thread bufferWritingThread;

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
