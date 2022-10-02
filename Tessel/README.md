# Tessel

Tessel is a fast, scalable kMer tiling program that creates the kMer sets (and kMer pairs) used by Blue. 
It takes any number of sequence data files (fasta or fastq), tiles them for kMers, and writes a set of the 
distinct kMers present in the reads files, and how many times each kMer was found.

Tessel usually writes these kMers+counts to a '.cbt' file that contains binary (2-bits/base) canonical
kMers. Canonical here means that only the lexicographically lowest of a kMer and its reverse complement 
is kept. This is useful for real sequence data where reads can come from either the forward or 
reverse strand, and saves programs that want to look for a kMer in a kMer set from having to check for 
both its as-read and reverse-complement forms. It is also possible to get Tessel to only track 'as-read' kMers,  
rather than converting them to canonical form, and this is useful when tiling genomes and assembled contigs.
Tessel also supports a number of other text-based output formats, both for debugging and 
for compatibility with JellyFish.

Tessel keeps two counts for each kMers - a separate count for both the kMer as it appears in the hash table and its 
reverse-complement. Keeping these two counts separate lets Blue check for unbalanced kMers/reads, often a sign of sequencing errors 
or remnant adapter sequences. 

As of v2.1, Tessel has subsumed the GenerateMerPairs program and, by default, will also tile the sequence data reads for kMer 'pairs'. These
are pairs of short kMers, separated by a gap, and are effectively longer kMers. They are used by Blue to
partially resolve the ambiguity that comes from using shorter kMers.

v2.2 is the result of porting Tessel to the .NET platform. The only significant change in this release is that Tessel can now process gzipped sequence data files 
directly. 

Tessel is a command-line program with the following cryptic usage hint:  
```
Tessel -k kMerLength -g genomeLength [-t #threads] [-tmp tempDir] [-min minCount] [-trim nn] [-trimQual nn] [-s] [-pairs] [-nopairs] [-canonical|asread] [-text textFN] [-textFormat pairs|sum|faPairs|faSum] cbtName readsFN... or readsPattern
```  

A typical run of Tessel is something like… 

`
Tessel -k 25 -g 10M -t 16 -tmp /tmp/DH10B -min 2 -trimqual 30 DH10B MiSeq_Ecoli_DH10B_110721_PF_R?_paired.fastq  
`  

…which will tile whatever files match MiSeq_Ecoli_DH10B_110721_PF_R?_paired.fastq for 25-mers and write out a file of all the distinct
canonical kMers found and their counts. The kMers will be written to DH10B_25.cbt. Tessel produces a number of temporary files and 
these will be written to the /tmp/DH10B directory. Reads will be tail-trimmed before tiling, using a qual score of 30 as the acceptable-quality 
level. Tessel will also tile the reads again for kMer pairs, and these will be written to DH10B_25.prs.

The full set of Tessel parameters is:

| Option/Parameter        | Description       |
| ------------- |-------------| 
| -k *kMerSize* |Length of the kMers to be generated. 'k' should typically be >=20 and must be <=32. The default is k=25. | 
| -h or -help | Writes out a more verbose help text.
| -g or -genome *genomeSize* | A guess at the total size of the final (assembled) genome or genomes in the sample. This is used for the initial sizing of some arrays and any rough guess will probably be OK. Default is 10,000,000 (10M) which will handle bacterial genomes quite well. K,M&G suffixes can be used. |  
| -t or -threads *number_of_threads* | Number of parallel tiling threads. Tessel scales well with the number of threads, up to the point where it is limited by reading and writing the files. Default is 1. |  
| -tmp tmpDir | Temporary files are written to this directory. This can be used to improve performance by pushing these temporary files to the fastest available storage. This is especially helpful with remote HDD storage and less important with local SSDs. Default is to write these temp files to the current directory. |
| -min *minCount* | kMers found fewer that minCount times will be dropped. Useful for reducing the size of the kMer tables, especially when Blue will probably discard low abundance kMers while loading the kMer tables. The default is to discard singletons (-min 2). | 
| -trim *trimLength* | The tails of all reads will be trimmed by this length before tiling. -trimQual is a better tail-trimming technique in most cases. Default is no trimming. |
| -tq or -trimQual *minimum_qual_value* | Tail-trims each read prior to tiling and reduces the number of kMers coming from the error-prone read tails. Uses a sliding window to cut back the tails of reads until the specified quality level is reached. Default is no trimming. |
| -s | Recursively search through sub-directories for matching file names. Default is to only look in the current directory if no explicit directory name is given as part of the readsFN parameter.  |
| -pairs | Tile the reads for paired short kMers. This option replaces the use of GenerateMerPairs the Blue v1 pipeline. This is the default setting. |
| -nopairs | Do not tile the reads for kMer pairs. |
| -canonical | Regard kMers and their reverse complements as equivalent. Separate counts are maintained for both forms of each kMer. Default option. Causes a .cbt file to be used for the binary kMers+counts. |
| -asread | kMers are saved in the form that they are found, and not converted into a canonical form. This option is intended only for use on contigs/genomes. Causes a .abt file to be used for the binary kMers+counts.|
| -text *textFN* | The following formats tell Tessel how to write out the kMers and their counts to a text file, as well as in binary format (a .cbt or .abt file). The format of the text file is specified by the -textFormat option |
| -textFormat ...| Specifies the format used when writing kMers+counts to the text file. The ‘sum’ and ‘faSum’ options are compatible with Jellyfish text file formats. | 
| -textFormat *pairs* | One line per kMer: kMer (tab) count (tab) rc-count. This is the default text format. | 
| -textFormat *sum* | One line per kMer: kMer (tab) summed-count. |
| -textFormat *faPairs* | FASTA file. Header is >count rcCount. Sequence line is kMer. |
| -textFormat *faSum* |	FASTA file. Header is >summed-count. Sequence line is kMer. |
| *cbtname* | This parameter is used to construct the output file names. The binary-formatted tiled kMers+counts will be written to cbtName_kmerLength.cbt (e.g. DH10B_25.cbt) and the kMer repetition depth histogram will be saved as cbtName_kmerLength_histo.txt (e.g. DH10B_25_histo.txt). If this parameter is simply a name, such as ‘DH10B’, then the output files will be placed in the current directory. If it includes a directory, such as ‘Healed/DH10B’, then the file will be written to this directory. If the ‘-asread’ option is used, a ‘.abt’ file suffix is used instead of ‘.cbt’. The 'pairs' file, if requested' is written using the same file name, but with a suffix of '.prs'. |
| *readsFNP*  | List of file names/patterns of reads to be tiled. e.g. MiSeq_Ecoli_DH10B_110721_PF_R?_paired.fastq. There can be any number of these files names for file name patterns. |

The canonical binary kMers are written to a ‘.cbt’ file, and using the ‘-asread’ option causes them to be written to a ‘.abt’ 
file. The format of these files is:
- Int32	‘k’ length used. The top 8 bits of this Int32 are reserved for use as a file format marker (currently 0).  
- Sets of {uint64, int32, int32} triplets, one for each kMer.  
... The uint64 value is a kMer in packed binary format.   
...	Two bits per base (A=00, C=01, G=10, T=11).  
...	kMers are left-adjusted within the 64-bit word.  
...	A kMer can be no bigger than 32 bases long.   
... The first int32 is the number of times the corresponding kMer was seen.  
... The next int32 is the number of times it was seen in its reverse complement form.

The 'pairs' file format is similar:
- Int32	'gap' in bases between the pair of 16-mers making up the 'pair'. The top 8 bits of this Int32 are reserved for use as a file format marker (currently 0).  
- Sets of {uint64, int32} pairs, one for each kMer pair.  
... The uint64 value is a kMer pair in packed binary format. Top 32-bits are the first 16-mer, low 32-bits are the second 16S-mer in the pair.  
... The int32 is the number of times the corresponding kMer pair was seen. 


Tessel is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, and these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. Tessel is ‘installed’ simply by copying its code file to an appropriate directory on your system. 

You can compile Tessel yourself using the `dotnet publish` command. You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk).  
The Tessel code itself is in Program.cs in this directory, and you'll also need to download the files
in WorkingDocsCoreLibrary. Tessel can be built as 'frame-work dependent' code or as 
standalone code (AOT) with necessary run-time code linked into the Tessel executable. The AOT code will run on systems that do not have the 
.NET run-time installed. AOT code generation is only supported from .NET7 onwards.

The type of executable produced by `dotnet publish` is controlled by the `PublishProfile` option. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. Small scripts are provided that will 
build Tessel executables. The AOT builds have to be done on a system that is compatible with the intended execution targets as 
parts of the platform run-time are linked into the executables. Pre-built Tessel code is provided for Windows and Linux, and 
.NET SDKs are available that will allow Tessel to be built for both x64 and ARM macOS systems. The Linux code has been built on 
Ubuntu 18 and tested on Ubuntu 22 and SUSE LES 15. 

The command `dotnet publish ./Tessel.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml` will build a
framework-dependent x64 Linux Tessel executable, and other versions can be built by changing the name of the profile file in the 
publish command.





