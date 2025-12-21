# Pup

Pup is a (usually) symmetric scalable gene/genome comparer. It takes sets of gene or genomes and does a full NxN similarity comparison of them by 
counting shared kMers. The primary output is an NxN table of shared kMers counts
but symmetric percentage PIMs and sparse (tabbed pairs) PIM formats are also generated. 
Pup is usually run with 'inexact' matching, allowng for substitutions, insertions and deletions, 
much like the seed+extension matching model used in aligners. It can also do exact kMer matching, or allow substitutions only ('exact' and 'subs'). 

Unlike conventional 'count shared kMers' comparators, Pup effectively corrects each source sequence as it is matched to a target sequence.
As a result, Pup doesn't fall off a cliff as quickly with taxonomic distance. A single base mismatch will result in a single missed match, not 'k' consecutive mismatches. It also is not overly optimisic, in contrast to ANI, 
as it looks at the whole gene/genome, not just at cherry-picked bits already determined to be similar or likely to be similar. 

Pup testing has shown good linearity down to quite low similarity, very comparable to traditional MSA tools used on the 
same (small) datasets, and it can scale from comparing genes to whole genomes. 
```
usage: Pup [-k kMerSize] [-s] [-dirs|-files|-singleFile] [-exact|-inexact|-subs] [-t #threads] [-source sourceFNP] -out countsFN genomeFNs/Dirs

The input (target) sequences can come in 3 confusing forms:
	-dirs	- a set of directories. One directory per genome. Each directory contains a set of sequence data files and each target is
		  all the sequences found in a directory. Useful for comparing draft genomes/metagenomes. 
	-files	- a set of files (default). One file per genome. Each file contains a number of sequences (each file being a genome/draft genome)
	-single	- file(s) containing a number of sequences (used for gene-level comparison). The sequences within the file(s) are cross-compared

The -s option allows for a recursive search (following down subdirectories) when finding matching files. Default is to only look in current directory

The comparisons can be exact or allowing for mismatches (subs or inexact (subs/dels/ins)). The default is '-inexact'

The default behaviour is an NxN comparison from all the (target) genes/genomes onto the same set of target genes/genomes.
The -source option allows for the source and destinations to be different, useful if matching a set of genomes against a (different) set of reference genomes.

The result (-out) is a matrix of counts of shared kMers, between each source organism and each target organism. 
These results come in 3 formats: counts, symmetric similarity PIM (percent identity matrix), sparse similarity PIM. 
The symmetric PIM is comparable to an EBI MSA PIM, and the sparse PIM can be used to generate a tree using 'usearch -cluster_aggd'.

There are many other options, many undocumented, that set parameters for the comparisons. Some of the more useful ones are: ones are
	-k		sets the kMer length (default k=21)
	-sl		the seed length used during scanning for a matching kMer (default 11)
	-mm		the number of mismatches allowed before an extension is abandoned (default 10)
	-t		sets the number of concuurent threads used for comparing sequences (default is half the available processor threads)
	-minlength	Don't compare sequences shorter than this length (default 0). 
```
If you are not using the '-source' option, the command line syntax is quite Unix-friendly and will work as expected even if the shell does file name pattern expansion.
Pup is expecting some number of options (with '-' prefixes) and a list of file names or file name patterns, and is indifferent as to whether the 
shell turns any patterns into lists of file names. 
The '-source' option takes a file name pattern as its parameter and a Unix shell will usually turn this pattern into file names as well. 
The first of these files names will be chosen as the subject and all other files coming from the sourceFNP pattern will be treated as additional target files, 
not source files. The solution is to tell the shell not to expand (glob) the sourceFNP parameter by quoting it ("sourceFNP"). The
file name pattern will then be passed intact (unexpanded) to Pup which will use it the find the names of the source files. 

Some examples:

> Pup -k 21 -t 24 -s -files -out Aspergillus_FO_inexact_GxG.txt Aspergillus_oryzae_genomes/\*.fna Aspergillus_flavus_genomes/\*.fna (takes a day or two to compare the whole genomes)

> Pup -k 21 -t 16 -inexact -singlefile c23oTests.fa -out c23o_k21_inexact_v130_gxG.txt (compares a set of c23o genes, coming from a single file)

> Pup -k 21 -t 16 -inexact -files -out Methanoculleus_k21_inexact_GxG.txt Methano*.fna

## Installing and building Pup

Pup is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, but these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. The Linux code has been built on both
Ubuntu 20 (glibc 2.31) and 24 (glibc 2.39) and tested on Ubuntu 24 and SUSE LES 15.

Pup is ‘installed’ simply by copying the desired code file to an appropriate directory on your system. For example, if you want to 
run the self-contained .Net8 Pup code file on a Linux system compatible with Ubuntu 24 (glibc 2.39), you should copy the 
Pup file from the Linux64DN8FDU24 directory.You may have to set permission bits on the Pup code file before you can execute it on Linix. 

Compilation scripts are provided for both Windows (.ps1) and Linux (.sh). For example, `BuildPup_Linux64DN8AOTU24.sh` 
will build a stand-alone Pup executable targeting .NET8. This script was run on Ubuntu 24 to build the executable, and the resulting Pup code will be 
expecting glibc 2.39 or above to be available when it is run. 

You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk) and these scripts assume the directory structure
found in the GitHub repository, as shown below. 
The Pup code itself is in Pup.cs in the Pup directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs to the WorkingDocsCoreLibrary. The Pup compilation scripts are intended to be run from the 
WorkingDogs/Pup_Clancy/Pup directory.

```
WorkingDogs
	Pup+Clancy
		PupClancyMatching.cs
		Pup.sln
		Pup
			Pup.cs
			Pup.csproj
			Properties
				PublishProfiles
					<XXX>.pubxml
	WorkingDogsCoreLibrary
		kMers.cs
		kMerCollections.cs
		Sequence.cs
		SeqFiles.cs
```

The type of executable produced is controlled by the PublishProfile (.pubxml) parameter supplied to `dotnet publish`. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. 
AOT builds for Linux have to be done on a system that is 
compatible with the intended execution platform as 
the required glibc level is embedded into the executable code.

You can build Pup for other target systems, such as MacOS on Arm, by creating a 
.pubxml file with the appropriate RuntimeIdentifier - see https://learn.microsoft.com/en-us/dotnet/core/rid-catalog. 
The command `dotnet publish ./Pup_v2.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml` 
in the `BuildPup_Linux64DN8AOTU24.sh` script is what builds the
framework-dependent x64 Linux Pup executable, and other versions can be built by changing the name of the profile file in this 
publish command.