using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ExtractHeadersFromFasta
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Length != 2)
            {
                Console.WriteLine("usage: ExtractHeadersFromFasta fastaFN headersFN");
                return;
            }

            string fastaFN = args[0];
            string headersFN = args[1];

            if (!File.Exists(fastaFN))
            {
                Console.WriteLine("fasta file not found: " + fastaFN);
                return;
            }

            StreamReader fasta = new StreamReader(fastaFN);
            StreamWriter headers = new StreamWriter(headersFN);
            bool EOF = false;

            while (!EOF)
            {
                string line = fasta.ReadLine();
                if (line == null)
                    break;

                if (line != "" && line[0] == '>')
                {
                    int spaceIdx = line.IndexOf(' ');
                    if (spaceIdx > 0)
                        line = line.Substring(0, spaceIdx);
                    headers.WriteLine(line);
                }
            }

            headers.Close();
        }
    }
}
