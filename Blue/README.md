# Blue

Blue is a kMer-based DNA error correction tool. Blue v1 is described in _Blue: correcting sequencing errors using consensus and context_. _Bioinformatics. 2014;30:2723-32_. The significant difference between Blue v1 and Blue v2 is a change in intent, moving from reducing the number of errors in a set of reads to producing only correct reads. Blue v1 would correct whatever errors it could correct, and at times would leave uncorrected errors at the end of reads. The intent was to correct what could be safely corrected, and pass the resulting improved reads on to other tools. The _good_ parameter in Blue v1 mitigated this behaviour somewhat by discarding reads with insufficient numbers of _good_ kMers after correction, ensuring that the most erroneous and uncorrectable reads were dropped rather than being passed on to the next tool in a pipeline.

Blue v2 tries to ensure that all the reads it are completely corrected, and uses tail trimming to remove uncorrectable parts of reads and drops any reads that appear to still contain errors after correction and trimming. Blue is quite averse to rewriting reads. Incremental kMer-based algorithms can often infer which base should follow any given kMer, and such extended reads can quickly diverge from the read actually being corrected. This produces chimeric corrected reads whose start comes from the read being corrected and whose tail comes from elsewhere in the genome. Blue avoids this rewriting by tracking the number of corrections being made to a read, and will abandon a read if too many corrections are being made too close together. Blue will also abandon a read if reaches a point where there are no valid next kMers. Blue v1 would pass on these abandoned but partially corrected reads, after checking that the read had sufficient good kMers. Blue v2 will trim such reads back to the last good kMer, and discard it if it is now too short (as defined by the repurposed _good_ parameter).

The overall result of this change, and numerous others, is much improved accuracy (in relative terms). On the E. coli DH10B dataset referred to in the paper, Blue now has 99.98% of the corrected (and possibly trimmed) reads aligning with zero changes against the reference sequence, up from 99.90% previously (with _–good 80_ used in both cases).

Blue v2.2 is a port onto the .NET platform. This port changes the way Blue is compiled/built, and also allows the direct processing of gzipped sequence data files (no need to unzip first).
v2.2.1 fixed a bug in the writing of .gz files, and imporved support for correcting metagenomes with extremely low coverage reads (use -min 1 in Tessel and Blue). Also added -unzip option as writing .gz files is slow.

The first step in using Blue is to generate a set of kMers (and kMer pairs) using Tessel. Tessel will also generate a histogram of the kMer repetition frequency that you can use to set the _minReps_ parameter for your Blue run. Once you have a set of kMers (a .cbt file), and (preferably) a set of corresponding kMer pairs (a .prs file), you can go ahead and correct your reads. Blue is a command-line program with the following cryptic usage hint:
```
Blue [-help] [-r run] –m minReps [-hp] [-t threads] [-l length] [-fixed]
[-variable] [-good nn%] [-problems] [-extend nn] [-output dir]
[-paired] [-unpaired] cbtFN readsFNs or patterns
```
A typical run of Blue is something like…
```
Blue –r hv80 –m 20 –t 16 –g 80 –mq 25 –v DH10B_25.cbt MiSeq_Ecoli_DH10B_110721_PF_R?_paired.fq
```
… which will correct the pair of reads files that matched MiSeq_Ecoli_DH10B_110721_PF_R?_paired.fq, writing corrected reads to MiSeq_Ecoli_DH10B_110721_PF_R1_paired_hv80.fq and MiSeq_Ecoli_DH10B_110721_PF_R2_paired_hv80.fq, using the kMers found in DH10B_25.cbt and the pairs from DH10B_25.prs. The read pairing in the starting reads will be maintained in the corrected reads. The expected minimum depth for a good read is set to 20. Corrected reads will be trimmed back to the final good kMer (if the read cannot be completely corrected), and the trimmed reads must be at least 80% of their starting length (no more than 20% of the read can be trimmed away). The reads will be left at their trimmed length (if they were trimmed). The error prone tail of the region will be estimated using a qual score of 25.

The full set of Blue parameters is:

| Option | Parameter | Description       |
| --- | --- | --- |
| -r or -run | _runName_ | This is the name of this healing run. It is appended to the name of each of the reads files to produce the name of the corresponding file of corrected reads. For example, with _-r h_ the corrected form of ERR022075_1.fastq will be ERR022075_1_h.fastq. The default value is 'corrected_minReps'. |
| -s or -stats | _statsFN_ | This is the name of the file to be used to save the statistics from this run of Blue. Blue will generate a default statistics file from the first reads file name parameter and the runName, e.g. ERR022075_h_stats.txt. |
| -t or -threads | _noOfThreads_ | This option says is how many parallel threads to use for the tiling. The default is to use just 2 threads. Blue scales well with the number of threads, but will eventually be limited by the time needed to read and write the sequence files. There is also some per-thread memory allocations, so using more threads will also use more memory. |
| -m or -min | _minReps_ | This is the minimum kMer repetition depth expected in good kMers and is used (rarely) to distinguish between _good_ and _bad_ kMers. In general, this should be set to somewhere in the dip between the LHS error spike and the start of the Poisson curve derived from the good reads (more on this later). This value isn't very critical as the good/poor/bad depths are calculated dynamically for all reads if they contain any 'good' kMers at all. |
| -max | _maxReps_ | Reads with an average depth of coverage above this level will not be corrected. This rarely useful option was added to better handle a metagenomic dataset with a very deep-coverage repeat region that was best ignored during correction and assembly. |
| -l or -length | _maxLength_ | Reads are trimmed to this length before correction. This option may be of use when correcting poor quality data with very noisy tails, such as some MiSeq 300-mer reads. The minQual parameter may also prove to be useful with this type of noisy data. The default is to do no length-based trimming. |
| -mq or -minqual | _minQual_ | Sets the minimum qual score to be used when estimating the starting point of the noisy tail in each read. Once Blue reaches this point when correcting a read, it stops trying quite so hard, reducing the recursive depth allowed during tree exploration and only checking for substitution errors for Illumina-style data. The effect of this option is to allow Blue to successfully correct more of the most error-ridden reads. The default value is a qual score of 30. |
| -hp |   | The –hp option tells Blue to try even harder to when finding errors that could be corrected. There are a number of tests done at every base position, all based on depth of coverage. These tests will pick up random indel errors, but indels are so common at the end of homopolymer runs in 454 and IonTorrent data that multiple different hp run lengths all look to be OK as well. For example, if our genome had AAAAAA then with Illumina data this is what we'd see almost all the time. With 454-like data, we_d probably get 5 or 7 As as frequently as 6 As, so depth of coverage test would indicate that none of them were errors. The –hp flag causes Blue to look for the end of hp runs and forces an attempt at correction at that point, regardless of depth of coverage. Blue will always attempt to detect and correct indel errors regardless of the setting of this option, and setting –hp has a performance cost so it should not be set for Illumina data. |
| -g or -good | _%length_ | The meaning of this parameter has changed in Blue v2. It now specifies the minimum read length required after correction and trimming (and possible extension). The default value is 70 which means that corrected reads must be at least 70% of their original length after trimming. Blue maintains the pairedness of its input files, so if one read of a pair fails the _good_ test, both the failing read and its mate will be dropped (but written to a _\_singles_ file in case you want to use them). Those reads that fail this minimum length test will be written to a _problems_ file if you have set the _-problems_ option. |
| -o or -output | _outputDir_ | This specifies the output directory where Blue will write all its files, including the corrected reads. By default, these files are written to the directory where the reads were found. |
| -problems |   | Asks Blue to save any uncorrectable reads in a _\_problems_ file for later examination. |
| -v or -variable |   | With Blue v2, read lengths are allowed to change during correction as a result of insertions, deletions and the trimming of uncorrected read tails. By default, reads are written just as they are at the end of correction/trimming, and so can vary in length. |
| -fixed |   | This option asks Blue to try to extend any trimmed reads to get them closer to their starting lengths. This extension works by starting at the end of a trimmed read and trying to extend it unambiguously, one base at a time. The extension process is stopped if two or more viable next kMers are found at any point. This extended length is the one that is checked against the minimum length set through the _good_ option. |
| -fixPadded |   | The extended reads produced by _-fixed_ can still vary in length as the extension process will stop as soon as a read cannot be unambiguously extended. This option forces all corrected reads to be exactly the same length as they were prior to correction by extending them and then padding them out with Ns to the required length. This option is provided for compatibility with any tools that require all reads to be the same length.   |
| -extend | _No. of bases to extend_ | This option gets Blue to try to unambiguously extend all corrected reads by the given number of bases, regardless of whether they were trimmed or not. This option can result in longer corrected reads, but the extension process is not perfect. |
| -paired |   | By default, if Blue is asked to correct pairs of reads files, it will maintain strict pairing between the reads written to each of the output files. If the _-_good__ option is used, then Blue will discard a complete pair of reads if any either of them fails the goodness test. |
| -unpaired |   | This option tells Blue not to treat pairs of reads files as a pair set. The main use of this option is with the –_good_ option as it allows a set of unpaired files to be corrected in a single run of Blue, without having good reads from one file discarded because the corresponding (but unrelated) read in another file was uncorrected. |
| cbtFN |   | This the name of the kMer file to be used (produced by Tessel). Blue will also look for a .prs file with same name and use it if it finds it. |
| reads patterns or FNs |   | These parameters specify the names of the sequence data files to be corrected. You can either supply a list of space separated files names or a filename pattern. On Windows you would normally use a pattern and let Blue turn it into a set of matching file names. On Linux, the same pattern will normally be turned into a list of file names by the shell, with equivalent results. |


Blue is written in C# and is provided pre-compiled for Windows and Linux (and can be built for macOS). The AOT versions of these code files 
are stand-alone and should not require the installation of any additional run-time libraries. Smaller framework-dependent (FD) code
files are also provided, and these need to have an appropriate .NET run-time installed. See https://learn.microsoft.com/en-gb/dotnet/core/install
for instructions. Blue is ‘installed’ simply by copying its code file to an appropriate directory on your system. The Linux .NET 6 & 7 code files
have been compiled on Ubuntu 18, and the .NET 8 code under Ubuntu 22. If you have glibc version issues, try one of the other versions.

You can compile Blue yourself using the `dotnet publish` command. You’ll need to have installed the appropriate .NET SDK (see https://learn.microsoft.com/en-us/dotnet/core/sdk).  
The Blue code itself is in Program.cs in this directory, and you'll also need to download the files
in WorkingDocsCoreLibrary. Blue can be built as 'frame-work dependent' code or as 
standalone code (AOT) with necessary run-time code linked into the Blue executable. The AOT code will run on systems that do not have the 
.NET run-time installed. AOT code generation is only supported from .NET7 onwards.

The type of executable produced by `dotnet publish` is controlled by the `PublishProfile` option. Profiles are held in the 
Properties/PublishProfiles directory, for both framework-dependent and AOT compilations. Small scripts are provided that will 
build Blue executables. The AOT builds have to be done on a system that is compatible with the intended execution targets as 
parts of the platform run-time are linked into the executables. Pre-built Blue code is provided for Windows and Linux, and 
.NET SDKs are available that will allow Blue to be built for both x64 and ARM macOS systems. The Linux code has been built on 
Ubuntu 18 and tested on Ubuntu 22 and SUSE LES 15. 

The command `dotnet publish ./Blue.csproj -c release /p:PublishProfile=Linux64DN6FDFolderProfile.pubxml` will build a
framework-dependent x64 Linux Blue executable, and other versions can be built by changing the name of the profile file in the 
publish command.



