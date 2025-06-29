# Kelpie V2.3 

Kelpie V2.3 (2.3.4) is now available. 

* Better use of paired-reads to decide between plausible paths; 
* Multi-threads the final (extension) phase for better performance; 
* Supports multiple forward or reverse primers (e.g.-f AGAGTTTGATCMTGGCTCAG,TTCYGKTTGATCCYGSCRGA); 
* Extended sequences can be written out dereplicated (with size= annotations). 

# Kelpie usage

Kelpie usage is basically unchanged with this release, although there are new options to enable new functionality.
```
usage: Kelpie [-h] [-t #thrds] -f forwardPrimer -r reversePrimer [-filtered|-unfiltered] WGSReadsFNP extendedReadsFN (V2.0.10)
       -h                - Write out this extended help and exit
       -threads nn       - max parallel threads to allow. max will use all available (default)
       -f forwardPrimer e.g. GTGYCAGCMGCCGCGGTAA
       -r reversePrimer e.g. GGACTACNVGGGTWTCTAAT
       -filtered         - (default) WGS reads have been pre-filtered down to the genomic region of interest.
                           This (small) read set will be kept in memory for faster processing during filter construction.
       -unfiltered       - WGS reads were not pre-filtered and will be re-read from files as needed.
       WGSreadsFNP       - List of file names/patterns of reads to be processed. e.g. S1_270_reads_anonymous_R?_16S.fasta
       extendedReadsFN   - File name for extended between-primer reads. e.g. S1_270_reads_anonymous_16S_V4.fasta

Less commonly needed options (most can often be ignored)
       -strict           - kMers used for final reads filtering must be found in both files of a pair. (default)
       -loose            - Opposite of -strict, all kMers are use for filtering.
       -paired|unpaired  - Use -paired if you have paired-end reads, and -unpaired if you want to force pairs of reads files to be processed individually.
                           Default is to assume pairs of reads files are paired. 
       -mismatches m[+n] - Allow up to mm mismatches in the 3' or 5' ends of a primer. Default is 1. Only increase for long imprecise primers. AKA -mm.
       -min nn           - Minimum length accepted after assembly/extension, even if terminal primer not found. Default is that all extended 'amplicons' must finish with a terminating primer.
       -length nn[-mm]   - Expected length of targeted region/amplicon, including primers, e.g. 295 for 16SV4. Either a single length or a range can be specified.
       -mindepth nn      - kMers found fewer than this many times are dropped from the filtering table. Only added to help pull out full length 'amplicons' from a very deep dataset.
       -qualtrim nn      - FASTQ reads are quality trimmed on the RHS using a sliding window until this qual score is reached. Default is 30. AKA -qt & -tq.
       -noLCF            - no low-complexity filter. Low complexity kMers are usually not included in the region filters. This option lets them be included.
       -save primerTag   - Save the filtered/trimmed between-primer reads, and add this tag when building the file names. e.g. v4 --> S1_270_reads_anonymous_R?_16S_v4.fasta
       -primers          - Save the actual primer sequences found in the reads to XXXX_primers.txt.
       -derep            - Dereplicate the extended reads (with size= annotations)
       -tmp/-kept        - Used to improve efficiency of processing unfiltered WGS datasets.
       -log              - Debugging log. Use -trn nnn to generate additional trace on a selected read. 
```
For example: 

`Kelpie_v -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT CAMI_medium\S1-2_270_reads_anonymous_R?.fq CAMI_medium_16S_Kv2UF.fa -unfiltered -loose -mm 2 -tmp c:\SeqTemp\K2` 

`Kelpie_v2 -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT C:\SeqTemp\KelpieTests\Cel119_?.fastq Cel119_16S_v4_v234_US.fa -unfiltered -qt 0 -strict -tmp c:\SeqTemp\K2`

`Kelpie_v2 -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT Cel119_?_16S_20_fz_25_qt.fa Cel119_16S_v4_v2234_FS.fa -filtered -strict` 

`Kelpie_v2 -f TCNACNAAYCAYAARRAYATYGG -r TANACYTCNGGRTGNCCRAARAAYCA PlatesGH_2003\Sample_PRO1747_PlateG_*_R?_val_?_COI.fa PlatesGH_2003_FolmD_FF_strict.fa -filtered -strict`

`Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA PlatesGH_2003\Sample_PRO1747_PlateG_*_R?_val_?.fq PlatesGH_2003_BF3BR2_UF_strict.fa -unfiltered -strict -tmp C:\SeqTemp\K2 -kept I:\K2_kept`

## Installing and building Kelpie

Kelpie is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, and these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. Kelpie is ‘installed’ simply by copying its code file to an appropriate directory on your system. The Linux .NET 8 compilations were done using Ubuntu 22. The .NET 6 and 7 compilations
were done using Ubuntu 18 to target an older glibc version. 

You can compile Kelpie yourself using the `dotnet publish` command. You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk).  
The Kelpie code itself is in Program.cs in this directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs in WorkingDocsCoreLibrary. Kelpie can be built as 'frame-work dependent' code or as 
standalone code (AOT) with necessary run-time code linked into the Kelpie executable. The AOT code will run on systems that do not have the 
.NET run-time installed. AOT code generation is only supported from .NET7 onwards.

The type of executable produced by `dotnet publish` is controlled by the `PublishProfile` option. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. Small scripts are provided that will 
build Kelpie executables. The AOT builds have to be done on a system that is compatible with the intended execution targets as 
parts of the platform run-time are linked into the executables. Pre-built Kelpie code is provided for Windows and Linux, and 
.NET SDKs are available from Microsoft that will allow Kelpie to be built for both x64 and ARM macOS systems. The Linux code has been built on both
Ubuntu 20 (glibc 2.31) and 24 (glibc 2.39) and tested on Ubuntu 24 and SUSE LES 15. 

The command `dotnet publish ./Kelpie_v2.csproj -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml` will build a
framework-dependent x64 Linux Kelpie executable, and other versions can be build by changing the name of the profile file in the 
publish command.

### Advice on primers and filtering
The in-silico primer sequences used by Kelpie work in much the same way as their molecular equivalents. Primers can match into the expected conserved region, and they can 
match elsewhere as well. In a real PCR, the first run with the forward primer can result in off-target matches, but the following reverse primer run will select only those 
sequences that also contain the reverse primer. Kelpie does something similar, it first finds reads containing forward or reverse primers, and then finds reads that overlap these 
primer-containing reads. Because Kelpie can only work with the (short) reads it is given, not full-length amplicons, it can't discard off-target reads as easily, as it doesn't know anything 
about the target region. These off-target reads can result in large numbers of reads being kept during the filtering stage. These off-target reads result in larger kMer 
tables, but are otherwise innocuous. Off-target reads are worth controlling though, as they can degrade performance as they still have to be tiled and de-noised.

Kelpie simulates primer mismatches when finding primer-containing reads, and the number of mismatches allowed is set by the -mismatches (-mm) option. The default is 1 mismatch. 
The right number of mismatches is a balance between off-target matches and primer effectiveness. The -primers option will show what primers sequences were actually found in the WGS reads. 
I often do an initial run with a higher -mm value (2 or 3), and then adjust the primer sequences to better match what's in the data, and then reduce the -mm value to reduce off-target matches. 

The fastest way of running Kelpie is with pre-filtered reads. I use BuildFilter with a curated representative reference set (e.g. RefSeq 16S) and then run FilterReads to extract only 
those WGS reads that appear to be coming from my target region. If you are using FilterReads, you should use the -eitherinpair or -bothinpair options to keep the pairedness of the filtered reads,
now that Kelpie can use paired-reads to resolve alternative paths in the kMer tree. 

### Release notes

### V2.3.4

* Made the read extending/assembling phase multi-threaded and this is now much faster.
* Multiple forward and reverse primers supported. e.g -f AGAGTTTGATCMTGGCTCAG,TTCYGKTTGATCCYGSCRGA -r TACNGNTACCTTGTTACGACTT for extracting both archaeal and bacterial full-length 16S in a single run.
* Improvements to accuracy for low coverage and noisy data. 
* Better use of paired-reads and read coverage to avoid excessive tree exploration and reduce the risk of chimeras with noisy data.
* Paired-read code now ignores trivial pairs, coming from short pairs (poor DNA prep). 

### V2.2.0

* Paired-reads are now used to resolve some potential branches when choosing between multiple 'next' kMers. The list of possible next 
kMers is first trimmed by using long kMers, up to the length of a read, and if there are still multiple possible paths, these
alternative kMers are checked for paired-read support. If there are still multiple choices, these paths are explored recursively as before.
Paired-reads are most useful when extracting amplicons which are much bigger than any single read, such as 27F-1492R 16S. 
The paired-read code will only work with paired reads. If you are using a filter to extract target-region reads from a WGS
dataset prior to running Kelpie, you should make sure the filtered reads have kept their pairings. If you are usig FilterReads, 
you should be using the -bothinpair or -eitherinpair options.
* Adapter sequences are detected and trimmed. Adapter sequences are trimmed from starting reads at the beginning of the filter-building step. 
Subsequent filter-building iterations won't build on these primer+adapter sequences, effectively removiing them from the filter 
and the set of reads used for kMer table construction and extension. 
* The context-checking code has been revised, and it now checks from longest context to shortest. This change results in more branches
being resolved by contexts and slightly improved accuracy.

### V2.1.0

* Kelpie was ported to .NET6/7 with this release and the compilation process has changed as a result. Kelpie can now be built
using `dotnet publish` to produce either a framework-dependent or standalone code files. The standalone file contains all the 
necessary run-time support and no platform installation is needed. The framework-dependent code files are smaller but need to 
have a suitable .NET run-time installed before they can be run.
* .gz sequence data files are now supported directly. No need to unzip first.
* The min/max amplicon length range parameter has been replaced with -length (the expected length of the amplicons). This length is used to set a limit on the number 
of iterations when building the between-primers kMer set. It is really only needed with poorly-behaved primers that match other than the anticipated target region.
* Maximum mismatches can be set separately for the forward and reverse primers (e.g. -mm 1+2 allows 1 mismatch for the forward 
primer and 2 for the reverse)
* A few small bug fixes and tweaks have been included.


### V2.0.11

* Fixes a minor bug that could occur when total number of primer reads was < 100.
* Changes what happens when the -primers option is used. This option used to just write out a file of primers found/used, but 
now also adds the primers associated with each extended read to that read's header.
* Completes the move from kMer 'pairs' to 'contexts'. Pairs and contexts are identical for long kMers <= 64b. Any sequence data between the starting and ending 
kMers in a sequence were ignored with 'pairs', on the assumption that the hashed starting and ending kMers would be suficiently distinct and that any 
collisions would be resolved at the next length. One rare corner case uncovered in testing showed that this assumption could 
cause problems, and fully-inclusive 'contexts' were adopted as a result. 
* Trailing runs of at least 16 Gs are trimmed from reads, regardless of qual scores. Illumina produces runs of Gs as 'default' bases and one dataset had 
such runs with good qual scores. Runs of GGs within reads are untouched, although high-G kMers could still be discarded by the low-complexity filter.
* Trial extensions used to be processed in ACGT order. The highest depth variant is now tried first, and if it reaches a terminal proimer, its length
is used to set a maximum length for the exploration of the other variants. 

### V2.0.10

Kelpie V2.0.10 was the first public release of Kelpie V2. 

The V2 release added one major feature to V1 and substantial parts of the code were multi-threaded to improve performance. 
Kelpie V2 can now work directly with full (unfiltered) WGS datasets, starting with just a pair of primer sequences, rather than 
requiring reads to be pre-filtered. 
Iteratively building the between-primers region filter requires multiple passes over the reads. With unfiltered reads this 
means re-reading the (often large) WGS dataset files, and this process takes some time. 
The 83GB unfiltered MECSM W1 dataset used for testing takes about 10 minutes to process on a system with fast (NVMe SSD) local 
storage and 24 cores, while the same run using pre-filtered reads takes about 6s (processing 73 MB of 16S filtered data).

The '-unfiltered' and '-filtered' options are used to tell Kelpie whether it is processing an unfiltered WGS dataset, or if 
it has been given a much smaller pre-filtered-for-gene dataset that can be kept in memory. '-filtered' is the default, for 
compatibility with Kelpie V1 (and because it's so much faster). The first step in processing an unfiltered WGS dataset is to 
copy/convert it to a set of temporary files. These temporary files are in single-line FASTA format, and each partition (a file)
holds up to 5,000,000 WGS reads. These partitions are processed in parallel in the following processing steps. These temporary 
files should be placed on the fastest available storage (such as a local SSD). The location of the temporary files can be specified
using the '-tmp' option. These files are usually deleted at the end of a Kelpie run, but the '-kept' allows them to be left, 
and used again in subsequent Kelpie runs, saving some time. Setting the -kept option tells Kelpie where to look for the saved 
partitioned fasta files. If these files are found in the specified directory, they are copied over to the temp directory and 
Kelpie will proceed as usual. If the -kept directory is specified but doesn't contain the expected files, Kelpie will 
revert to generating the temp files from the specified WGS data files as usual, but it will not delete these temp files after 
they have been used to generate the set of region-filtered reads. These temp files 
can then be manually moved to the -kept directory in preparation for the next Kelpie run over the same data. 

The default number of mismatches allowed in a primer core (RHS) has been changed from 2 to 1 in Kelpie v2, reducing the number 
of erroneous not-quite-primer matches, with no impact of overall results (apart from greatly reducing the number of starting 
reads that could not be extended fully) and avoiding spurious read matches for unfiltered WGS datasets. The existing '-strict' 
and '-loose' options may also be use of in reducing the number of spurious 'matching' reads. The 
'-strict' option tightens the kMer filter used to produce the final set of 'filtered' reads by removing any kMers that aren't 
found in both the R1 and R2 files in a pair of WGS sequence data files. '-strict' is the default, but '-loose' can pick up some 
additional low-abundance sequences at times. 

Kelpie V2 also offers improved assembly accuracy, returns more fully-extended reads returned and goes further into the 
low-abundance tail of rarer sequences. These improvements come from extensive revisions to kMer table denoising, making use of
kMer contexts and deferred culling to better preserve low-abundance kMers that follow high abundance regions. Kelpie V2 has 
been extensively tested using the dataset mentioned in the Kelpie paper, and with COI data using a number of primers. 
Kelpie V2 has also been used to extract full-length 16S sequences from both cultures and metagenomic samples.

## Kelpie v1 (deprecated)

Kelpie extracts/assembles full-length inter-primer sequences from WGS metagenomic datasets. 
You could think of it as something akin to in-silico PCR, taking a pair of primer sequences 
and returning a set of full-length inter-primer reads. A paper describing Kelpie has been published 
in PeerJ (https://peerj.com/articles/6174/). 

Kelpie is a command-line program and is usually run as follows:
```
Kelpie -f forwardPrimer -r reversePrimer readsToFilterFNP extendedReadsFN  
where   forwardPrimer is the forward primer sequence  
        reversePrimer is the reverse primer sequence  
        readsToFilterFNP is a list of reads file names or a file name pattern  
        extendedReadsFN is the file name for the extended reads  
```
For example: `Kelpie -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT W1_?_16S_20_fz_25.fa W1_16S_v4.fa` 
                     
`Kelpie -help` will give the full set of parameters - but only the basic ones just shown are needed in most cases.

The input to Kelpie can be a set of pre-filtered WGS sequence data files or unfiltered (v2). Filtering here means that the full WGS has been processed 
in some way so as to discard reads from everywhere other than the general region of interest, such as the 16S gene. This filtering greatly 
reduces the size of the files that Kelpie has to process. This type of filtering is done by EBI using Infernal in HMM-only mode
against a library of ribosomal RNA models from Rfam, but can also be done quickly using FilterReads (another member of the 
WorkingDogs pack) that you will find in this GitHub repository.

The Coal Seam Microbiome datasets referenced in the Kelpie paper are available for download from the CSIRO Data Access Portal
at DOI:10.4225/08/5b31ca6373d48. Both the full and 16S-filtered reads are available there for all 3 samples (W1, W2 and W3). 
The command example above can be used to run Kelpie on the W1 dataset, using the EMP 16S V4 primers and producing a set of
extended reads in the file W1_16S_v4.fa.







