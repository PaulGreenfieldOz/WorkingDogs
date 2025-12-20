
# BuildFilter
BuildFilter builds kMer filters (.mer files) for use by FilterReads. These kMer filters are simply sets of (usually) canonical, binary-packed kMers tiled from a set of (reference) sequences.  

BuildFilter is a command-line program with this cryptic usage hint:  
`
BuildFilter [-v1|-v1c|-v2c|-v2a] -k merSize [+/-lcf] [-s] [-mindepth nn] [-minlength nn] kMersFN genesFN 
`

A typical run of BuildFilter is…  
```
BuildFilter -k 20 +lcf RDPv16+RefSeq_5-18_16S_NC_20.mer RDPv16+RefSeq_5-18_16S_NC.fa  
```

…which will tile the 16S reference sequences found in RDPv16+RefSeq_5-18_16S_NC.fa for 20-mers, and turn these into a set of distinct, canonical 20-mers, with low complexity kMers discarded.

Canonical means that the reverse-complement of every kMer is calculated, and only the lexicographically lowest of the kMer and its RC is kept to represent the kMer. For example, 
with a kMer of ACACACACACA and its RC of TGTGTGTGTGT, the ACACACACACA form will be kept and used to represent both variants. Keeping only canonical kMers makes it much easier 
to handle sequence data which could have been derived from either strand of a DNA molecule. 

BuildFilter now can generate either canonical or as-read kMers. This change led to versioning the file format: v1 is the original 
(canonical-only) format; v2 uses parts of the initial header word
in the .mer file to hold a file format type, and a canonical/as-read flag. V1 files will continue to be accepted by FilterReads, and V2 files are correctly parsed and interpreted. The file version 
and canonical/as-read tiling are controlled by the `-v` option. -v1 produces a file in the origainl format; -v1c produces an identical file, but the fact that is in canonical more is more obvious; -v2c 
produces a v2 file holding only canonical kMers (equivalent to a v1c file); and -v2a produces a v2 file holding only as-read kMers.

BuildFilter can also attempt to find and discard low-complexity reads/kMers. These are often not very useful or informative, and can result in overly-accepting filters.

The full set of BuildFilter parameters is:
- -k *kMer_length*  
Length of kMers in the filter. This would normally be at least 20 to take advantage of the distinctiveness of such kMers. kMers cannot be longer than 32 bases.
[-v1|-v1c|-v2c|-v2a]
optional file format and kMer type. Default is -v2c. Only provided for backwards compatibility.
- +lcf   
Turns on low-complexity filtering of the kMers/reads. This is the default.
- -lcf   
Turns off low-complexity filtering of the kMers/reads. Not usually desirable but needed in cases such as iterative gene-filtering.
- -mindepth nn   
kMers found fewer than this many times (across all sequences) will be discarded.
- -minlength nn  
Reads shorter than this length will not be tiled for kMers
- -s   
Search for matching file names in any subdirectories as well as in the current directory. Default is to only look for matching file names in the current directory. 
- *kMersFN*  
(optional) Name for the generated kMers file. If no kMersFn name is given, one is constructed from the first (and only in this case) seqsFN by adding _kMer_length.mer at its end. 
- *seqsFNP*  
List of file names/patterns of sequences to be tiled for kMers. e.g. RDPv16+RefSeq_5-18_16S_NC.fa
There can be any number of file names (or file name patterns), and all of these will be tiled to generate a single kMer file.  

## Installing and building BuildFilter

BuildFilter is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, but these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. The Linux code has been built on both
Ubuntu 20 (glibc 2.31) and 24 (glibc 2.39) and tested on Ubuntu 24 and SUSE LES 15.

BuildFilter is ‘installed’ simply by copying the desired code file to an appropriate directory on your system. For example, if you want to 
run the self-contained .Net8 BuildFilter code file on a Linux system compatible with Ubuntu 24 (glibc 2.39), you should copy the 
BuildFilter file from the Linux64DN8FDU24 directory.You may have to set permission bits on the BuildFilter code file before you can execute it on Linix. 

Compilation scripts are provided for both Windows (.ps1) and Linux (.sh). For example, `BuildBuildFilter_Linux64DN8AOTU24.sh` 
will build a stand-alone BuildFilter executable targeting .NET8. This script was run on Ubuntu 24 to build the executable, and the resulting BuildFilter code will be 
expecting glibc 2.39 or above to be available when it is run. 

You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk) and these scripts assume the directory structure
found in the GitHub repository, as shown below. 
The BuildFilter code itself is in BuildFilter.cs in the BuildFilter directory, and you'll also need to download 
kMers.cs, SeqFiles.cs and Sequence.cs to the WorkingDocsCoreLibrary. The BuildFilter compilation scripts are intended to be run from the 
WorkingDogs/BuildFilter directory.

```
WorkingDogs
	BuildFilter
		BuildFilter.sln
		BuildFilter
			BuildFilter.cs
			BuildFilter.csproj
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

You can build BuildFilter for other target systems, such as MacOS on Arm, by creating a 
.pubxml file with the appropriate RuntimeIdentifier - see https://learn.microsoft.com/en-us/dotnet/core/rid-catalog. 
The command `dotnet publish ./BuildFilter_v2.sln -c release /p:PublishProfile=Linux64DN10AOTFolderProfile.pubxml` 
in the `BuildBuildFilter_Linux64DN8AOTU24.sh` script is what builds the
framework-dependent x64 Linux BuildFilter executable, and other versions can be built by changing the name of the profile file in this 
publish command.

