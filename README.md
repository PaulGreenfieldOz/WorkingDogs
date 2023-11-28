### WorkingDogs

The home repo for the Working Dogs bioinformatics tools: Blue, Pup, Kelpie and other members of the pack. These tools are all kMer-based and make use of the very useful almost-always-unique properties of moderate-length kMers (k=20+).

Blue is a DNA read error corrector based on kMers and their context. It was published as "Blue: correcting sequencing errors using consensus and context" in Bioinformatics. 2014; 30(19), 2723-2732. 
Release 2.0 appeared here around the middle of 2018.

Pup is a kMer-based NxN whole gene/genome comparator uses techniques from error correction to get more-accurate-than-BLAST similarity numbers for sequences ranging from genes to whole genomes.
Clancy is a specialised clone of Pup meant for comparing (1:N) annotated genes to genomes.

Kelpie is a specialised assembler for between-primer regions, taking WGS reads and producing the equivalent of amplicons. It was published
in PeerJ in January 2019. https://peerj.com/articles/6174/. 

The initial release of Kelpie has been replaced by Kelpie_v2, providing both new features and improved performance/efficacy. Unfiltered WGS reads files can now be used directly, although
pre-filtering reads is still recommended where possible for performance reasons. Much of the code is now multi-threaded as well. Accuracy and the ability to extract/assemble
very low abundance sequences have also been improved.

FilterReads is a boringly-named kMer filter than can rapidly find reads of interest in very large datasets. Given a 16S filter, it produces results very similar to HMM-based
tools, but runs at close to IO speed.

Other kMer-based bioinformatics tool will join the pack from time to time. The GHAP amplicon pipeline will appear soon, and Pup will make its appearance soon too.

All code is now being compiled with .NET, producing both framework-dependent and pre-packed (AOT) code files. .NET 6, 7 & 8 are all supported. The framework-dependend
code files need to have the corresponding .NET run-time installed, while the AOT code has any necessary run-time code bundled into
the code file. The usual Linux glibc dependency issues has led to both Ubuntu 18 and Ubuntu 22 being used for compilations. If none of the pre-compiled Linux code files work on your system, you can 
install the latest supported .NET SDK and build your own bespoke code files. Later .NET releases support MacOS, both x64 and ARM) but you will have to compile this yourself.