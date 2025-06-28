
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

BuildFilter is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, and these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. BuildFilter is ‘installed’ simply by copying its code file to an appropriate directory on your system. 

You can compile BuildFilter yourself using the `dotnet publish` command. You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk).  
The BuildFilter code itself is in Program.cs in this directory, and you'll also need to download the files
in WorkingDocsCoreLibrary. BuildFilter can be built as 'frame-work dependent' code or as 
standalone code (AOT) with necessary run-time code linked into the BuildFilter executable. The AOT code will run on systems that do not have the 
.NET run-time installed. AOT code generation is only supported from .NET7 onwards.

The type of executable produced by `dotnet publish` is controlled by the `PublishProfile` option. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. Small scripts are provided that will 
build BuildFilter executables. The AOT builds have to be done on a system that is compatible with the intended execution targets as 
parts of the platform run-time are linked into the executables. Pre-built BuildFilter code is provided for Windows and Linux, and 
.NET SDKs are available that will allow BuildFilter to be built for both x64 and ARM macOS systems. The Linux code has been built on 
Ubuntu 20 (for glibc 2.31) and Ubuntu 24 (for glibc 2.39). and tested on Ubuntu 24 and SUSE LES 15.5. 

The command `dotnet publish ./BuildFilter.csproj -c release /p:PublishProfile=Linux64DN8FDFolderProfile.pubxml` will build a
framework-dependent x64 Linux BuildFilter executable, and other versions can be built by changing the name of the profile file in the 
publish command.

