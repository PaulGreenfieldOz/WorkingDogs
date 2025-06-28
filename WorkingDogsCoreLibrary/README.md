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
The files in this library are referenced in the .csproj files for each code project. These references are created in Visual Studio via Add Existing Item (and I add as link rather than copying the code into the project directory).
