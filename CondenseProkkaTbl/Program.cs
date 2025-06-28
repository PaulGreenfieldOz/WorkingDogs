using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using WorkingDogsCore;

namespace CondenseProkkaTbl
{
    // Takes a .tbl (gene calls) file and turns into a table for easier comparison
    // Can also take a reference sequence for the 'same' organism and give a kMer-based identity score for each gene

    class Program
    {
        // >Feature NODE_1
        // 2030	1374	CDS
        // 			EC_number	5.3.1.6
        // 			gene	rpiA_1
        // 			inference	ab initio prediction:Prodigal:2.6
        // 			inference	similar to AA sequence:UniProtKB:Q7MHL9
        // 			locus_tag	BJKLFMJH_00001
        // 			product	Ribose-5-phosphate isomerase A
        // 2773	2174	CDS
        // 			EC_number	6.3.3.2
        // 			gene	ygfA
        // 			inference	ab initio prediction:Prodigal:2.6
        // 			inference	similar to AA sequence:UniProtKB:P0AC28
        // 			locus_tag	BJKLFMJH_00002
        // 			product	5-formyltetrahydrofolate cyclo-ligase
        // 3282	2974	CDS
        // 			gene	zapA
        // 			inference	ab initio prediction:Prodigal:2.6
        // 			inference	similar to AA sequence:UniProtKB:P0ADS2
        // 			locus_tag	BJKLFMJH_00003
        // 			product	Cell division protein ZapA
        // 3590	4165	CDS
        // 			inference	ab initio prediction:Prodigal:2.6
        // 			locus_tag	BJKLFMJH_00004
        // 			product	hypothetical protein
        // 39603 38440	CDS
        //          EC_number	3.1.1.31
        //          db_xref     COG:COG2706
        //          gene        pgl_1
        //          inference   ab initio prediction:Prodigal:002006
		//          inference   similar to AA sequence:UniProtKB:O34499
        //          locus_tag   MLKJEHBF_00117
        //          product	    6-phosphogluconolactonase 

        //>MLKJEHBF_00001 Putative sugar phosphate isomerase YwlF
        //ATGAAAATTGCGATCGGGTCCGATCACGCAGGCTTTCAATACAAGGAGAAAATAAAAGAG
        //TTTCTGCAGGAACAAGGTCATGTTGTGAACGATTTCGGCACCTTTTCGGAGGAGCCGGTT
        //GACTATCCATTGTTCATTCGCCCCGTCGCCGAGGCAGTGGCGCGTGGCGAATACGAGCGC
        //GGGATCGTCCTGGGTGGTTCAGGCAATGGTGAAGCGATGGCGGCAAACAGGATCAAGAAA
        //GTGCGTTGCGCGCTGTGTTGGAACGAAGAGCTCGCGCGTCTCGCGCGTCAGCACAATGAC
        //TCCAACGTCCTTTCGCTGGGCCAACGGGTCATCTCCGAGGAAACTGCCCTCGAGATTGTA
        //CGCGTGTGGTTGGAAACGGCATTCGAAGGTGGTCGTCACGTGCGGCGGATTGAATTGCTG
        //GACGCGGATTGA

        static void Main(string[] args)
        {
            if (args.Length != 2 && args.Length != 4)
            {
                Console.WriteLine("usage: CondenseProkkaTbl prokkaTblFN [prokkaFfnFN refFN] tableFN");
                return;
            }

            const int k = 25;
            string prokkaTblFN = args[0];       // .tbl to be summarised
            string prokkaFfnFN = null;          // file with seqs for each called gene
            string refFN = null;                // genome for comparison
            string tableFN;
            if (args.Length == 4)
            {
                prokkaFfnFN = args[1];
                refFN = args[2];
                tableFN = args[3];
            }
            else
                tableFN = args[1];

            if (!File.Exists(prokkaTblFN))
            {
                Console.WriteLine("prokka .tbl file not found: " + prokkaTblFN);
                return;
            }

            if (prokkaFfnFN != null && !File.Exists(prokkaFfnFN))
            {
                Console.WriteLine("prokka .ffn file not found: " + prokkaFfnFN);
                return;
            }

            if (refFN != null && !File.Exists(refFN))
            {
                Console.WriteLine("reference seq file not found: " + refFN);
                return;
            }

            StreamReader prokkaTbl = new StreamReader(prokkaTblFN);
            StreamWriter table = null;
            try
            {
                table = new StreamWriter(tableFN);    
            }
            catch
            {
                Console.WriteLine("could not open table file: " + tableFN);
                return;
            }

            Dictionary<string, float> geneIdentity = null;
            if (prokkaFfnFN != null)
                geneIdentity = MapGenesToReference(prokkaFfnFN, refFN, k);

            table.WriteLine("contig" + "\t" + "locus_tag" + "\t" + "start" + "\t" + "end" + "\t" + "reverse" + "\t" + "gene" + "\t" + "UniProtKB" + "\t" + "EC_number" + "\t" + "COG" + "\t" + "ref%" + "\t" + "product");

            Console.WriteLine("Summarising " + prokkaTblFN);

            bool EOF = false;
            char[] spaceSeparator = new char[] { ' ' };
            char[] tabSeparator = new char[] { '\t' };
            char[] colonSeparator = new char[] { ':' };

            string contig = null;
            int start = 0;
            int end = 0;
            bool reverse = false;
            string EC_number = "";
            string gene = "";
            string UniProtKB = "";
            string locus_tag = "";
            string product = "";
            string COG = "";

            while (!EOF)
            {
                string line = prokkaTbl.ReadLine();
                if (line == null)
                    break;
                if (line == "")
                    continue;

                if (line.StartsWith(">Feature"))
                {
                    // starting a new contig - so write the last gene from the previous contig if there is one
                    if (locus_tag != "")
                    {
                        string identity = "";
                        if (geneIdentity != null)
                        {
                            if (geneIdentity.ContainsKey(locus_tag))
                                identity = geneIdentity[locus_tag].ToString("F3");
                        }
                        table.WriteLine(contig + "\t" + locus_tag + "\t" + start + "\t" + end + "\t" + (reverse ? "reverse" : "") + "\t" + gene + "\t" + UniProtKB + "\t" + EC_number + "\t" + COG + "\t" + identity + "\t" + product);
                        locus_tag = "";
                    }

                    string[] splitNode = line.Split(spaceSeparator);
                    contig = splitNode[1];
                    continue;
                }

                string[] splitLine = line.Split(tabSeparator);

                if (splitLine.Length == 3 && (splitLine[2] == "CDS" || splitLine[2] == "tRNA" || splitLine[2] == "tmRNA" || splitLine[2] == "repeat_region" || splitLine[2] == "CRISPR" || splitLine[2] == "rRNA"))
                {
                    // starting a new gene - so write the previous one
                    if (locus_tag != "")
                    {
                        string identity = "";
                        if (geneIdentity != null)
                        {
                            if (geneIdentity.ContainsKey(locus_tag))
                                identity = geneIdentity[locus_tag].ToString("F3");
                        }
                        table.WriteLine(contig + "\t" + locus_tag + "\t" + start + "\t" + end + "\t" + (reverse ? "reverse" : "") + "\t" + gene + "\t" + UniProtKB + "\t" + EC_number + "\t" + COG + "\t" + identity + "\t" + product);
                    }
                    start = Convert.ToInt32(splitLine[0]);
                    end = Convert.ToInt32(splitLine[1]);
                    if (start > end)
                    {
                        reverse = true;
                        int swap = start;
                        start = end;
                        end = swap;
                    }
                    else
                        reverse = false;

                    EC_number = "";
                    gene = "";
                    UniProtKB = "";
                    locus_tag = "";
                    COG = "";
                    if (splitLine[2] == "CRISPR")
                        locus_tag = "CRISPR";
                    product = "";

                    continue;
                }

                if (splitLine[3] == "EC_number")
                {
                    EC_number = splitLine[4];
                    continue;
                }

                if (splitLine[3] == "gene")
                {
                    gene = splitLine[4];
                    continue;
                }

                if (splitLine[3] == "db_xref" && splitLine[4].StartsWith("COG:"))
                {
                    //          db_xref     COG:COG2706
                    COG = splitLine[4].Substring("COG:".Length);
                    continue;
                }

                if (splitLine[3] == "inference")
                {
                    // inference	similar to AA sequence:UniProtKB:P0ADS2
                    // inference	protein motif:HAMAP:MF_00208
                    string[] splitInference = splitLine[4].Split(colonSeparator);
                    if (splitInference.Length == 3 && splitInference[1] == "UniProtKB")
                        UniProtKB = "UniProt:" + splitInference[2];
                    if (splitInference.Length == 3 && splitInference[1] == "HAMAP")
                        UniProtKB = "HAMAP: " + splitInference[2];
                    continue;
                }

                if (splitLine[3] == "locus_tag")
                {
                    locus_tag = splitLine[4];
                    continue;
                }

                if (splitLine[3] == "product")
                {
                    product = splitLine[4];
                    continue;
                }
            }

            // write the last gene from the last contig if there is one
            if (locus_tag != "")
            {
                string identity = "";
                if (geneIdentity != null)
                {
                    if (geneIdentity.ContainsKey(locus_tag))
                        identity = geneIdentity[locus_tag].ToString("F3");
                }
                table.WriteLine(contig + "\t" + locus_tag + "\t" + start + "\t" + end + "\t" + (reverse ? "reverse" : "") + "\t" + gene + "\t" + UniProtKB + "\t" + EC_number + "\t" + COG + "\t" + identity + "\t" + product);
            }

            table.Close();
        }

        private static Dictionary<string, float> MapGenesToReference(string prokkaFfnFN, string refFN, int k)
        {
            Console.WriteLine("Tiling " + refFN + " for " + k + "-mers");
            Dictionary<string, float> geneIdentity = new Dictionary<string, float>();

            StreamReader refGenome = new StreamReader(refFN);
            bool EOF = false;

            HashSet<ulong> kMersFromRef = new HashSet<ulong>();
            ulong[] kMersInSeq = new ulong[10000000];
            bool[] kMersValid = new bool[10000000];

            while (!EOF)
            {
                string contig = SeqFiles.ReadRead(refGenome, SeqFiles.formatFNA);
                if (contig == null)
                    break;

                int kMerCount = kMers.GenerateMersFromRead(contig, k, ref kMersInSeq, ref kMersValid);

                for (int i = 0; i < kMerCount; i++)
                {
                    if (kMersValid[i])
                    {
                        ulong kMer = kMersInSeq[i];
                        ulong kMerRC = kMers.ReverseComplement(kMer, k);
                        if (kMerRC < kMer)
                            kMer = kMerRC;
                        if (!kMersFromRef.Contains(kMer))
                            kMersFromRef.Add(kMer);
                    }
                }
            }
            refGenome.Close();
            kMersInSeq = null;
            kMersValid = null;

            Console.WriteLine("Mapping genes in " + prokkaFfnFN);
            StreamReader calledGenes = new StreamReader(prokkaFfnFN);

            while (!EOF)
            {
                string geneName;
                string geneString = SeqFiles.ReadRead(calledGenes, SeqFiles.formatFNA, out geneName);
                if (geneString == null)
                    break;

                Sequence geneSeq = new Sequence(geneString);
                int kMerCount = geneSeq.Length - k + 1;
                int kMerMatches = 0;

                for (int i = 0; i < kMerCount; i++)
                {
                    ulong kMer;
                    bool kMerValid = Sequence.CondenseMer(geneSeq, i, k, out kMer);
                    if (kMerValid)
                    {
                        ulong kMerCanonical = kMer;
                        ulong kMerRC = kMers.ReverseComplement(kMer, k);
                        if (kMerRC < kMer)
                            kMerCanonical = kMerRC;

                        if (kMersFromRef.Contains(kMerCanonical))
                            kMerMatches++;
                        else
                        {
                            // variant mers (for 1-base substitutions)
                            List<ulong> kMerVariants = new List<ulong>(k * 4);
                            kMerVariants.Add(kMerCanonical);
                            int variantCount = kMers.GenerateMerSubVariants(kMer, kMerVariants, k);
                            foreach (ulong kMerVariant in kMerVariants)
                            {
                                ulong kMerVariantRC = kMers.ReverseComplement(kMerVariant, k);
                                kMerCanonical = kMerVariant;
                                if (kMerVariantRC < kMerVariant)
                                    kMerCanonical = kMerVariantRC;

                                if (kMersFromRef.Contains(kMerCanonical))
                                {
                                    geneSeq.Replace(i, kMerVariant, k);
                                    break;
                                }
                            }
                        }
                    }
                }

                geneName = geneName.Substring(1);                   // remove the leading >
                int spaceIdx = geneName.IndexOf(' ');               // find the end of the name
                if (spaceIdx > 0)
                    geneName = geneName.Substring(0, spaceIdx);     // and extract just the name (and discard the product)

                float identity = (float)kMerMatches / (float)kMerCount;

                geneIdentity.Add(geneName, identity);
            }

            return geneIdentity;
        }
    }
}
