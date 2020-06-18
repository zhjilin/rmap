# README


## INTRODUCTION

This code enables a relatively fast mapping of structured RNA Binding Protein(RBP) binding motifs to a large genome set.
However, it used a primitive algorithm to examine scores of all possible loci in the given sequence.

Any N occurs in the sub-sequence are not considered, therefore, there are might be some boundary effect.

Two modes can be used for 1) calculation of kmer count-score table 2) search the potential positions in the sequence with high score.

By default the program will perform the task 2), which gives scores to all possible locus. One can switch to the kmer count-score mode by using '-f'.
Keep in mind, it will yield a massive text file if you run it for the whole genome. It might be useful to filter position with low scores.


##DEPENDENCY

The code relies on the following libraries

### 1) BOOST
https://www.boost.org
### 2) kseq.h from Heng Li
https://github.com/attractivechaos/klib/blob/master/kseq.h

## COMPILING

In order to compile this program, one needs:

### for mac:
>=clang-800.0.42.1
BOOST >=106300

### for linux:
>=g++ 4.4.7
BOOST >=105900

```
g++ -std=c++11 main.cpp -o rmap -lboost_program_options -lz
```
This will generate a binary file 'rmap' in the current folder'

—-implemented for macOS and CentOS 6 (other platform were not tested).

INSTALLATION

The binary file rmap can be placed at desired location after compiling successfully.

### ACKNOWLEDGEMENTS

I appreciate the help from Andreas at Björn Högberg lab for discussion, and the access of kseq.h by Heng Li.


## USAGE

```
./rmap --help
```
Usage:rmap [options]

Allowed options:
--help                 produce help message
--fasta arg            input file
--list arg             specify a list of motifs
--motif arg            specify one motif file
-k [ --kmer ] arg (=8) specify the kmer length
-g [ --gap ] arg (=0)  specify the kmer length
-f [ --kmerflag ]      flag to switch to the kmer counting-score procedure
-r [ --revcomp ]       reverse complement sequence
-o [ --output ]        output the kmer count-score
-c [ --corr ]          output the correlation



### An quick example command to run the program:

1) For only one motif

```
./rmap -r -k 8 -g 10 --fasta test.fa --motif rbp.matrix
```

2) For a list of motifs
```
./rmap -r -k 8 -g 10 --fasta test.fa --motif list.rbp.txt
```




## ADDITIONAL SCRIPTS USED FOR ANALYSIS AND VISUALIZATION

Major scripts involved in the analysis visualization were deposited in the directory Scripts-HTR

### 0.FamilyDistribution
The count of RBP in each defined family (based on RBD)
Table downloaded with the following criteria from CISBP-RNA database: Need to check the database again. 

The R script is used to generate the figure presented in Fig.1C.

### 1.MotifPresentation
To visualize the base-pairing of the (predicted) stem-loop motifs.
The perl script used to generate the shade at the base position, which is involved to participate in forming the stem.

### 2.VennDiagram

RNAcompete and RNA bind-n-seq were combined with HTR-SELEX data to indicate the coverage of RBPs with motif characterized in each RBP family.
The R script was used to generate venn diagram presented in Fig.1B.

### 3.LengthDist
The length distribution of stem and loop from structured motifs (STEM and LOOP separately) in Fig.2F .

### 4.MotifComposition 

The dominant base at each position in pwm matrix is used to generate the overall nucleotide composition distribution for each motif (Fig. 4G)

### 5.InformationContent
Information content were calculated and visualized for comparison

### 6.MotifComparison

Scripts to place the motifs from different experiment for the same RBP together. 
For RNAcompete we used the motifs directly. However, due to the limitation of access to pwm matrix for RNBS, we used the available logos (downloaded from ENCODE portal) only.

### 7.RNAfold

The RNA sequence used to predict the secondary structure by RNAfold for the illustration in Fig.2B

### 8.Dimer

The script to visualize the proportion of RBPs with dimeric binding pattern in pie chart.

### 9.MotifMatchesEnrichment
MetaPlot_PerMotif-splicing180419.r. The script to generate metaplot presented in Fig 6 and Supplemental_Data_S1-S2. The required input can be found on Supplemental_Data_S4
Density2GO-tab-sp2.r. The script used to get GO Enrichment (R package clusterProfiler) using Ensembl gene ID queried through biomaRt package.
MetaPlot_PerMotif-distance2020April05.r. The script to generate the distance to feature histrograms Supplemental Data S3.
