 # Introduction
Working Dogs is a collection of (mostly) kMer-based tools for doing useful things with DNA sequence data. 

This library contains code common to all the Working Dogs pack of tools, including Blue, Pup, Kelpie and the imaginatively named FilterReads. Much of this code is concerned with kMers (short DNA strings of length 'k'). These kMers are processed in two forms: char strings or packed 2-bits per base into a ulong. The packed forms are much more efficient, but are limited to <= 32-mers. 

The library consists of 9 files:  

*kMers.cs*                
Code for generating and manipulating kMers, in be in either string and packed binary form.  
*kMerTables.cs*  
Loading large kMer tables from files (produced by Tessel).  
*kMerPairs.cs*            
Base code for handling packed pairs of short kMers, derived from the ends of longer DNA sequences.  
*kMerCollections.cs*      
Code for large, fast tables of kMers. Somewhat based on Dictionary and Hashset but lock-free, thread-safe and much more specialised.  
*TrimExtendReads.cs*      
kMer-based code for trimming and extending reads, using kMer depth thresholds.  
*Sequence.cs*             
Sequence is like StringBuilder, but for DNA, and can include qual data as well as bases. These DNA strings are mutable (support for sub, ins, del, truncate etc).  
*SeqFiles.cs*             
Code for fast reading and writing of files of DNA sequence data, in both FASTA and FASTQ format.  
*NCBITaxonomy.cs*         
Support code for loading and following NCBI taxonomy trees.  

# Getting Started
This library is built as a NuGet package and meant to be loaded via a NuGet reference into tools such as Blue and Kelpie. The NuGet package is available at https://www.nuget.org/packages/WorkingDogsCoreLibrary/
Its various files can be added to .Net Framework projects via Add Existing Item (and I add as link rather than copying the code into the project directory).
You can also just download/copy the libraries into a directory, along with the calling program (such as Kelpie) and compile using something like `csc /out:Kelpie.exe /target:exe *.cs`.
There are no prerequisites or dependencies, apart from Mono (currently), Visual Studio or .NET Core (future), that you’ll need in order to get a C# compiler.

# Build and Test
The CoreLibrary code can either be used through a project reference to the NuGet package, or it can be turned into a .Net Core class library that is referred to from with a project, or you can add the source code into existing (or new) DNA tool projects. 
Apart from writing your own tests for any changes, you can test using Blue and Kelpie. Blue uses most of the functionality of this library. Some methods are present but are not used by the currently-released WD tools. Sample code using most of these methods/classes is available on request.

# Contribute
This library is a set of useful functions common to the WD collection of tools. Feel free to make use of it in your own C# bioinformatics tools, or to suggest fixes for any bugs that are still present. 
