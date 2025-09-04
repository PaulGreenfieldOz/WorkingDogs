# Pup
Pup is a symmetric gene/genome comparer. It takes sets of gene og genomes and does a full NxN comparison of them. The primary output is shared kMers
but symmetric percentage PIMs and in sparse (tabbed pairs) formats. Pup can do exact kMer matching, or allow for subs or subs/ins/dels. Unlike
a strainght 'count common kMers' comparator, Pup effectively corrects the source seqeunces as it is matched to a target sequence.
As a result, Pup doesn't fall off a cliff as quickly. A single base mismatch will result in a single missed macth, not 'k' mismatches. 

Pup testing has shown good linearity down to quite low similarity, very comparable to traditionl MSA tools used on the same datasets.

usage: Pup -k kMerSize [-s] [-dirs|-files|-singleFile] [-exact|-inexact|-subs] [-t #threads] [-source sourceFNP] -out countsFN genomeFNs/Dirs

	 The input (target) sequences can come in 3 confusing forms:
	 -dirs        a set of directories, each containing a set of files containing one or more sequences
	              a directory for each genome, and a file for each contig/sequence (such as a .fna file)
	 -files       a set of files, each containing a number of sequences (each file being a genome/draft genome)
	              one file per genome, each containing a number of sequences
	 -single      a single file containing a number of sequences (used for gene-level comparison)
	              sequences within the file are compared
	
	 The -s option allws for a recursive search (following down subdirectories) for matching files.
	 
	 The comparisons can be exact or allowing for differences (subs or inexact (subs/dels/ins)).
	
	 The default behaviour is an NxN comparison from all the (target) genomes onto the same set of target genomes.
	 
	 The -source option allows for the source and destinations to be different, useful if matching a set of genomes against a set of reference genomes.
	
	 The result is a matrix of each organism compared to each other organism in 3 formats: counts, symmetric similarity PIM, sparse similarity PIM
	
	