# Kelpie V2 (V2.0.11)

Kelpie V2 is now available. This release adds one major feature to V1, and also improves overall efficiency and efficacy. 
Substantial parts of the code are also now multi-threaded to improve performance. 

Kelpie V2 can now work directly with full WGS datasets, starting with just a pair of primer sequences, rather than requiring reads to be pre-filtered. 
Iteratively building the between-primers region filter requires multiple passes over the reads. With unfiltered reads this means re-reading the (often large) WGS dataset files, 
and this process takes some time. 
The 83GB unfiltered MECSM W1 dataset used for testing takes about 10 minutes to process on a system with fast (NVMe SSD) local storage and 24 cores, while the same run using pre-filtered reads
takes about 6s (processing 73 MB of 16S filtered data).

The new '-unfiltered' and '-filtered' options are used to tell Kelpie whether it is processing an unfiltered WGS dataset, or if it has been given a much smaller 
pre-filtered-for-gene dataset that can be kept in memory. '-filtered' is the default, for compatibility with Kelpie V1 (and because it's so much faster). 
The first step in processing an unfiltered WGS dataset is to copy/convert it to a set of temporary files. These temporary files are in single-line FASTA format, and each partition (a file)
holds up to 5,000,000 WGS reads. These partitions are processed in parallel in the following processing steps. These temporary files should be placed on the fastest
available storage (such as a local SSD). The location of the temporary files can be specified using the new '-tmp' option. These files are usually deleted at the end of a Kelpie run, but 
the '-kept' allows them to be left, and used again in subsequent Kelpie runs, saving some time. Setting the -kept option tells Kelpie where to look for the saved partitioned fasta files. If these files 
are found in the specified directory, they are copied over to the temp directory and Kelpie will proceed as usual. If the -kept directory is specified but doesn't contain the expected files, Kelpie will 
revert to generating the temp files from the specified WGS data files as usual, but it will not delete these temp files after they have been used to generate the set of region-filtered reads. These temp files 
can then be manually moved to the -kept directory in preparation for the next Kelpie run over the same data. 

The default number of mismatches allowed in a primer core (RHS) has been changed from 2 to 1 in Kelpie v2, reducing the number of erroneous not-quite-primer matches, with no impact of overall results (apart from
greatly reducing the number of starting reads that could not be extended fully) and avoiding spurious read matches for unfiltered WGS datasets. 
The existing '-strict' and '-loose' options may also be use of in reducing the number of spurious 'matching' reads. The '-strict' option 
condenses the kMer filter used to produce the final set of 'filtered' reads by removing any kMers that aren't found in both the R1 and R2 files in a pair of WGS
sequence data files. '-strict' is the default, but '-loose' can pick up some additional low-abundance sequences at times. 

Kelpie V2 also offers improved assembly accuracy, returns more fully-extended reads returned and goes further into the low-abundance tail of rarer sequences. These improvements come from 
extensive revisions to kMer table denoising, making use of kMer contexts and deferred culling to better preserve low-abundance kMers that follow high abundance regions. Kelpie V2 has been extensively tested
using the dataset mentioned in the Kelpie paper, and with COI data using a number of primers. Kelpie V2 has also been used to extract full-length 16S sequences from both cultures and metagenomic samples. 

Kelpie usage is basically unchanged with this release, although there are now new options to enable new functionality.
```
usage: Kelpie [-h] [-t #thrds] -f forwardPrimer -r reversePrimer WGSReadsFNP extendedReadsFN (V2.0.10)
       -h                - Write out this extended help and exit
       -threads nn       - max parallel threads to allow. max will use all available (default)
       -f forwardPrimer e.g. GTGYCAGCMGCCGCGGTAA
       -r reversePrimer e.g. GGACTACNVGGGTWTCTAAT
       WGSreadsFNP       - List of file names/patterns of reads to be processed. e.g. S1_270_reads_anonymous_R?_16S.fasta
       extendedReadsFN   - File name for extended between-primer reads. e.g. S1_270_reads_anonymous_16S_V4.fasta

Less commonly needed options (most can often be ignored)
       -filtered         - (default) WGS reads have been pre-filtered down to the genomic region of interest.
                           This (small) read set will be kept in memory for faster processing during filter construction.
       -unfiltered       - WGS reads were not pre-filtered and will be re-read from files as needed.
       -strict           - kMers used for final reads filtering must be found in both files of a pair. (default)
       -loose            - Opposite of -strict, all kMers are use for filtering.
       -paired|unpaired  - Use -paired if you have paired-end reads, and -unpaired if you want to force pairs of reads files to be processed individually.
                           Default is to assume pairs of reads files are paired. Only used when cleaning final kMer filter table with -strict.
       -mismatches mm    - Allow up to mm mismatches in the 3' or 5' ends of a primer. Default is 1. Only increase for long imprecise primers. AKA -mm.
       -min nn           - Minimum length required after assembly/extension. Default is that extended 'amplicons' have to finish with a terminating primer.
       -max nn           - Stop extending reads after this length is reached. Default is to continue extending until a terminating primer is reached.
       -mindepth nn      - kMers found fewer than this many times are dropped from the filtering table. Only added to help pull out full length 'amplicons' from a very deep dataset.
       -qualtrim nn      - FASTQ reads are quality trimmed on the RHS using a sliding window until this qual score is reached. Default is 30. AKA -qt & -tq.
       -noLCF            - no low-complexity filter. Low complexity kMers are usually not included in the region filters. This option lets them be included.
       -save primerTag   - Save the filtered/trimmed between-primer reads, and add this tag when building the file names. e.g. v4 --> S1_270_reads_anonymous_R?_16S_v4.fasta
       -primers          - Save the actual primer sequences found in the reads to XXXX_primers.txt.
       -tmp/-kept        - Used to improve efficiency of processing unfiltered WGS datasets.
       -log              - Debugging log from read extension phase.
```
For example: 

`Kelpie_v -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT CAMI_medium\S1-2_270_reads_anonymous_R?.fq CAMI_medium_16S_Kv2UF.fa -unfiltered -loose -mm 2 -tmp c:\SeqTemp\K2` 

`Kelpie_v2 -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT C:\SeqTemp\KelpieTests\Cel119_?.fastq Cel119_16S_v4_v210_US.fa -unfiltered -qt 0 -strict -tmp c:\SeqTemp\K2`

`Kelpie_v2 -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT Cel119_?_16S_20_fz_25_qt.fa Cel119_16S_v4_v210_FS.fa -filtered -strict` 

`Kelpie_v2 -f TCNACNAAYCAYAARRAYATYGG -r TANACYTCNGGRTGNCCRAARAAYCA PlatesGH_2003\Sample_PRO1747_PlateG_*_R?_val_?_COI.fa PlatesGH_2003_FolmD_FF_strict.fa -filtered -strict`

`Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA PlatesGH_2003\Sample_PRO1747_PlateG_*_R?_val_?.fq PlatesGH_2003_BF3BR2_UF_strict.fa -unfiltered -strict -tmp C:\SeqTemp\K2 -kept I:\K2_kept`

Kelpie v2 can be compiled using the steps described below, and binaries created using mono mkbundle are provided in this repository for various Linux and OSX versions. Small changes/fixes
were also made to some of the WorkingDogsCoreLibrary files as part of Kelpie V2 development. (kMers.cs, kMerTables.cs, SeqFiles.cs, kMerCollections.cs). The revised versions of these files need 
to be downloaded if you are going to compile Kelpie_V2 from source. _

### V2.0.10
Kelpie V2.0.10 was the first public release of Kelpie V2. 

### V2.0.11

* Fixes a minor bug that could occur when total number of primer reads was < 100.
* Changes what happens when the -primers option is used. This option used to just write out a file of primers found/used, but now also adds the primers associated 
with each extended read to that read's header.
* Completes the move from kMer 'pairs' to 'contexts'. Pairs and contexts are identical for long kMer <= 64b. Any seqeunce data between the starting and ending 
kMers in a sequence were ignored with 'pairs', on the assumption that the hashed starting and ending kMers would be suficiently distinct and that any 
collisions would be resolved at the next length. One rare corner case
uncovered in testing showed that this assumption could cause problems, and fully-inclusive contexts were adopted as a result. 
* Trailing runs of at leas 16 Gs are trimmed from reads, regardless of qual scores. Illumina produces runs of Gs as 'default' bases and one dataset had 
such runs with good qual scores. Runs of GGs within reads are untouched, although high-G kMers could still be discarded by the low-complexity filter.
* Trial extensions used to be processed in ACGT order. The highest depth variant is now tried first, and if it reaches a terminal proimer, its length
is used to set a maximum length for the exploration of the other variants. 
### Kelpie_v2.1 Preview

The next planned Kelpie release will be v2.1. The only new feature planned for this release is directly supporting gzipped sequence data files. The Kelpie code will be 
ported to .Net 5 for this release, and make use of its cross-platform code capability to generate Linux and OSX code files. This release should be available by the end of 2020. 

## Kelpie v1 Release Notes 
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

The input to Kelpie is currently a set of pre-filtered WGS sequence data files. Filtering here means that the full WGS has been processed 
in some way so as to discard reads from everywhere other than the general region of interest, such as the 16S gene. This filtering greatly 
reduces the size of the files that Kelpie has to process. This type of filtering is done by EBI using Infernal in HMM-only mode
against a library of ribosomal RNA models from Rfam, but can also be done quickly using FilterReads (another member of the 
WorkingDogs pack) that you will find in this GitHub repository.

The Coal Seam Microbiome datasets referenced in the Kelpie paper are available for download from the CSIRO Data Access Portal
at DOI:10.4225/08/5b31ca6373d48. Both the full and 16S-filtered reads are available there for all 3 samples (W1, W2 and W3). 
The command example above can be used to run Kelpie on the W1 dataset, using the EMP 16S V4 primers and producing a set of
extended reads in the file W1_16S_v4.fa.

Kelpie is written in C# and is provided pre-compiled for Windows, OSX and common Linux variants. These pre-compiled code files 
should be stand-alone and should not require the installation of any additional run-time libraries. Kelpie is ‘installed’ simply by copying 
its code file to an appropriate directory on your system. 

You can compile Kelpie yourself using Mono (under OSX and Linux), or under Visual Studio on Windows. You’ll need the mono-devel 
package installed on Linux or OSX. The Kelpie code itself is in Program.cs in this directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs from WorkingDocsCoreLibrary. 
One simple way of compiling Kelpie (if you should need to) is to make a directory containing these 4 .cs files
and then run `csc /out:Kelpie.exe /target:exe *.cs` from this directory. Earlier versions of mono use ‘msc’ rather than ‘csc’ but 
the syntax of the compilation command is otherwise the same. The C# compiler will produce the executable Kelpie.exe. This code file can be run
natively on Windows, and on other platforms by using ‘mono Kelpie’ or by first cross-compiling using mkbundle to generate native code. 

The native Linux and OSX executables provided with the package were produced by cross-compiling Kelpie.exe using mkbundle, 
targeting various Linux and OSX releases (using runtimes downloaded from https://download.mono-project.com/runtimes/raw/). For example:  
	`mkbundle -o ubuntu-18.04\Kelpie Kelpie.exe --simple --cross mono-5.20.1-ubuntu-18.04-x64`   
You can also use mkbundle to generate a native executable for your current system. You’ll need to know where the mono-devel 
installation put the machine.config and the mono libraries. 
	`mkbundle -o Kelpie --simple Kelpie.exe --machine-config /etc/mono/4.5/machine.config -L /apps/mono/5.4.1.7/lib/mono/4.5`  
The mkbundle command has been changing with recent mono releases, and you may need to try other variants, such as:
	`mkbundle -o Kelpie --cross default Kelpie.exe`
	
It is expected that all use of mono will be replaced once .NET Core supports direct multi-platform code generation. 




