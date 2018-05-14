# Code Task

Given a FASTQ file, produce a list, sorted by frequency, of the "N" most frequent DNA k-mers (a substring of
length k) of length "K" in that file, along with the number of times they appear.

For reference and testing purposes, you can use the following FASTQ files from the 1000 Genomes Project: 
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01595/sequence_read/
```

## Design Decisions & Issues

1.	I used 64bit int variables for keeping char arrays to improve efficieny.
	I did not have time to implement it for 32bit systems.

2.	I did not control exceptions by try and catch method since i had no time
	for those things while i had so many ideas to improve the efficiency.
	Still could not try all ideas.
	
3.	I skipped the lines which have characters other than A,C,G and T.
	Therefore the results will probably be different than non-skipping version
	
4.	I read two articles about K-mer counting:
	First one	:	MSPKmerCounter: A Fast and Memory Efficient Approach for 
					K-mer Counting
	Second one	:	KMC 2:fast and resource-frugal k-mer counting
	
	I used their minimizer, superkmer and partition ideas.
	But i found a method myself which were not presented in articles. 
	Maybe it is a well-known method, did not see anywhere. My idea is about
	keeping three different histograms. Details are in the code. As a result for
	low topcount query, filtering by histograms gives so good performance.
	
5.	My program uses threads as many as available however the main issue is 
	reading the big file which is not a good part to use threads. Threads run
	only in hashtable and sorting part.
	
6.	I used the idea of not selecting minimizers which has prefix AAA,ACA and
	AA except at the beginning, however still not enough to make partitions
	uniform. I did not have time to implement choosing reverse minimizers.
	
### Prerequisites

What things you need to install the software and how to install them

```
Linux 64-bit
Windows 64-bit
```

### Installing

Just run make file to build for Linux 

```
% make
```
Or open project file in Visual Studio for Windows

## How To Run

K: length of substrings
N: most frequent substrings

```
% [executible] [FASTQfile] [K - length of substrings] [N - many most frequent substrings]
```

## Author

* **Mehmet SARIYUCE**



