# Kelpie

Kelpie extracts/assembles full-length inter-primer sequences from WGS metagenomic datasets. 
You could think of it as something akin to in-silico PCR, taking a pair of primer sequences 
and returning a set of full-length inter-primer reads. A paper describing Kelpie was published in PeerJ in
early 2019 (https://peerj.com/articles/6174/).

**This initial release of Kelpie was superseded by Kelpie_v2 in September 2020. This new release has its own directory in the
WorkingDogs repository.**

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

The input to Kelpie is a set of filtered WGS sequence data files. Filtering here means that the full WGS has been processed 
in some way so as to discard reads from everywhere other than the region of interest, such as 16S. This filtering greatly 
reduces the size of the files that Kelpie has to process. This type of filtering is done by EBI using Infernal in HMM-only mode
against a library of ribosomal RNA models from Rfam, but can also be done quickly using FilterReads (another member of the 
WorkingDogs mob) that you will find in this GitHub repository.

The Coal Seam Microbiome datasets referenced in the Kelpie paper are available for download from the CSIRO Data Access Portal
at DOI:10.4225/08/5b31ca6373d48. Both the full and 16S-filtered reads are available there for all 3 samples (W1, W2 and W3). 
The command example above can be used to run Kelpie on the W1 dataset, using the EMP 16S V4 primers and producing a set of
extended reads in the file W1_16S_v4.fa.

Kelpie is written in C# and is provided pre-compiled for Windows, OSX and common Linux variants. These pre-compiled code files 
should be stand-alone and should not require the installation of any additional run-time libraries. Kelpie is ‘installed’ simply by copying 
its code file to an appropriate directory on your system. 

You can compile Kelpie yourself using Mono (under OSX and Linux), or under Visual Studio on Windows. You’ll need the mono-devel 
package installed on Linux or OSX. The Kelpie code itself is in Program.cs in this directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs 
from WorkingDocsCoreLibrary. One simple way of compiling Kelpie (if you should need to) is to make a directory containing these 4 .cs files
and then run `csc /out:Kelpie.exe /target:exe *.cs` from this directory. Earlier versions of mono use ‘msc’ rather than ‘csc’ but 
the syntax of the compilation command is otherwise the same. The C# compiler will produce the executable Kelpie.exe. This code file can be run
natively on Windows, and on other platforms by using ‘mono Kelpie’ or by first cross-compiling using mkbundle to generate native code. 

The native Linux and OSX executables provided with the package were produced by cross-compiling Kelpie.exe with mkbundle, 
targeting various Linux and OSX releases. For example:  
	`mkbundle -o ubuntu-18.04\Kelpie Kelpie.exe --simple --cross mono-5.12.0-ubuntu-18.04-x64`   
You can also use mkbundle to generate a native executable for your current system. You’ll need to know where the mono-devel 
installation put the machine.config and the mono libraries. 
	`mkbundle -o Kelpie --simple Kelpie.exe --machine-config /etc/mono/4.5/machine.config -L /apps/mono/5.4.1.7/lib/mono/4.5`  
The mkbundle command has been changing with recent mono releases, and you may need to try other variants, such as:
	`mkbundle -o Kelpie --cross default Kelpie.exe`
	
It is expected that all use of mono will be replaced late in 2018 once .NET Core 3 arrives, with direct multi-platform code generation. 




