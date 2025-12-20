# Clancy

Clancy is an asymmetric gene to genome comparer. It takes a set of genes (e.g. an annotated genome) and compares each gene in
this source set to all of the genomes given as targets. 
The output is a table listing all of genes from the source, and how many kMers each of them share with each of the target genomes. 
Clancy was written to help with finding novel and conserved genes/regions.
For example, comparing an *Aspergillus flavus* reference genome to all available *A. flavus* and *A. oryzae* shows the 
presence/absence of the aflatoxin region across all the flavus/oryzae genomes 

Clancy and Pup share most of their code, so many of the same options are supported. The default is to do seed-based scanning 
and inexact (subs/ins/dels) matching with source 'correction'. 
```
usage: Clancy [-k kmerSize] [-exact|-inexact|-subs] [-t #threads] [-targets genes|genomes] -out countsFN [-s] sourceFN targetFNs

The first file name encountered in the parameter list is used as the source genome (sourceFN) and is expected to be a fasta file of 
annotated gene sequences.. It is often a .ffn file from an annotation pipeline run.  
The remaining file names are the targets. These are typically downloaded full genomes in FASTA format (-targets genomes) but can be files of 
called genes as well (-targets genes).
The -targets option says whether each target file is treated as a single collection of kMers or whether it is split into genes/contigs that are treated as distinct targets.

The -s option allows for a recursive search (following down subdirectories) for matching files.

The comparisons can be exact or allowing for differences (subs or inexact (subs/dels/ins)).

The result (-out) is a table of each source gene compared to each target. 

There are many other options, many undocumented, that set parameters for the comparisons. 
'-t' sets the number of concuurent threads used for comparing seqeunces. The default is half the available processor threads.
'-k' sets the kMer length (default k=21) 
'-sl' the seed length used during scanning for a matching kMer, and 
'-mm' the number of mismatches allowed before an extension is abandoned. The defaults for these parameters generally work well.

```

Some examples:

> clancy NM-KH31_C3.ffn -targets genomes Priestia_aryabhattai_genomes/\*.fna Priestia_megaterium_genomes\*.fna -k 21 -s -inexact -t 32 -out Clancy_NM-KH31_C3_genes_onto_PaPm.txt

> clancy -t 24 Aspergillus_flavus_annotations\ncbi_dataset\data\GCF_009017415.1\cds_from_genomic.ffn -targets genomes -inexact -s Aspergillus_flavus_genomes/\*.fna Aspergillus_oryzae_genomes/\*.fna -out A.flavus_NRRL_3357_to_A.flavus+oryzae_gxG.txt

> clancy -inexact -t 16 Methanoculleus_bourgensis_GCF_000304355.2_Mb_MS2.ffn -targets genomes *.fna -o Mb_MS2_onto_allMc_gxG.txt

## Installing and building Clancy

Clancy is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, but these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. The Linux code has been built on both
Ubuntu 20 (glibc 2.31) and 24 (glibc 2.39) and tested on Ubuntu 24 and SUSE LES 15.

Clancy is ‘installed’ simply by copying the desired code file to an appropriate directory on your system. For example, if you want to 
run the self-contained .Net8 Clancy code file on a Linux system compatible with Ubuntu 24 (glibc 2.39), you should copy the 
Clancy file from the Linux64DN8FDU24 directory.You may have to set permission bits on the Clancy code file before you can execute it on Linix. 

Compilation scripts are provided for both Windows (.ps1) and Linux (.sh). For example, `BuildClancy_Linux64DN8AOTU24.sh` 
will build a stand-alone Clancy executable targeting .NET8. This script was run on Ubuntu 24 to build the executable, and the resulting Clancy code will be 
expecting glibc 2.39 or above to be available when it is run. 

You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk) and these scripts assume the directory structure
found in the GitHub repository, as shown below. 
The Clancy code itself is in Clancy.cs in the Clancy directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs to the WorkingDocsCoreLibrary. The Clancy compilation scripts are intended to be run from the 
WorkingDogs/Clancy_Clancy/Clancy directory.

```
WorkingDogs
	Clancy+Clancy
		ClancyClancyMatching.cs
		Clancy.sln
		Clancy
			Clancy.cs
			Clancy.csproj
			Properties
				PublishProfiles
					<XXX>.pubxml
	WorkingDogsCoreLibrary
		kMers.cs
		Sequence.cs
		SeqFiles.cs
```

The type of executable produced is controlled by the PublishProfile (.pubxml) parameter supplied to `dotnet publish`. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. 
AOT builds for Linux have to be done on a system that is 
compatible with the intended execution platform as 
the required glibc level is embedded into the executable code.

You can build Clancy for other target systems, such as MacOS on Arm, by creating a 
.pubxml file with the appropriate RuntimeIdentifier - see https://learn.microsoft.com/en-us/dotnet/core/rid-catalog. 
The command `dotnet publish ./Clancy_v2.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml` 
in the `BuildClancy_Linux64DN8AOTU24.sh` script is what builds the
framework-dependent x64 Linux Clancy executable, and other versions can be built by changing the name of the profile file in this 
publish command.