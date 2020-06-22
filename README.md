# README


## INTRODUCTION

This code enables a relatively fast mapping of structured RNA Binding Protein(RBP) binding motifs to a large genome set.
However, it uses a relatively primitive algorithm to examine the scores of all possible loci in the given sequence.

Any Ns’ occuring in the sub-sequence are not considered, therefore, there might be some boundary effects.

Two modes can be used for 1) calculation of kmer count-score table 2) search the potential positions in the sequence with high score.

By default, the program will perform the task 2), which gives scores to all possible loci. Option ‘-f’ switches to the kmer count-score mode.
Note: Running for a whole genome will create a massive text file. Only one sequence is allowed for one run. It might be useful to filter positions with low scores.


## DEPENDENCY

The code relies on the following libraries

### 1) BOOST
https://www.boost.org
### 2) kseq.h from Heng Li
https://github.com/attractivechaos/klib/blob/master/kseq.h

## COMPILING

Compiling this program requires:

### for Mac:
>=clang-800.0.42.1
BOOST >=106300

### for Linux:
>=g++ 4.4.7
BOOST >=105900

```
g++ -std=c++11 main.cpp -o rmap -lboost_program_options -lz
```
This will generate a binary file 'rmap' in the current folder

—-implemented for macOS and CentOS 6 (other platform were not tested).

### INSTALLATION

The binary file rmap can be placed at any desired location after compiling successfully.

### ACKNOWLEDGEMENTS

I appreciate the help from Andreas Gådin for discussion, and the access to kseq.h by Heng Li.


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



### An example command to run the program:

1) For only one motif

```
./rmap -r -k 8 -g 10 --fasta test.fa --motif rbp.matrix
```

2) For a list of motifs
```
./rmap -r -k 8 -g 10 --fasta test.fa --motif list.rbp.txt
```




## ADDITIONAL SCRIPTS USED FOR ANALYSIS AND VISUALIZATION

Major scripts used in the analysis and visualization are deposited in the directory Scripts-HTR

### 0.FamilyDistribution
The count of RBP in each defined family (based on RBD, organism: homo sapiens)
Table downloaded from CISBP-RNA database (version 0.6). 

The R script FamDist18Oct.r is used to generate the figure presented in Fig.1C.

### 1.MotifPresentation
To visualize the base-pairing of the (predicted) stem-loop motifs.
The perl script shade_logo2.pl is used to generate the shade at the base position, which is involved to participate in forming the stem.

### 2.VennDiagram

RNAcompete and RNA bind-n-seq were combined with HTR-SELEX data to indicate the coverage of RBPs with the motif characterized in each RBP family.
The R script venn_diagram.R is used to generate Venn diagram presented in Fig.1B.

### 3.LengthDist
The length distribution of stem and loop from structured motifs presented (STEM and LOOP separately) in Fig.2F.

### 4.MotifComposition 

The script DrawTree17May.R is used to present all motifs based on the dendrograms (Fig. 4A-F). MarkTree-stemloop20180508.table uses a binary code to indicate the RBPs that can bind to structured motifs.


The dominant base at each position in pwm matrix is used to generate the overall nucleotide composition distribution for each motif(Fig. 4G). Scripts to generate the above information is count_freq.pl and PieChart-code.r.


### 5.InformationContent
Information content was calculated and visualized for comparison (Supplemental Fig. S8). Example command lines are illustrated in file cmd170517.sh.

### 6.MotifComparison

Scripts to place the motifs from different experiments for the same RBP together. 
For RNAcompete we used the motifs directly. However, due to the limitation of access to the pwm matrix for RNA bind-n-seq (RNBS), we used the available logos (downloaded from ENCODE portal) only.

### 7.RNAfold

The RNA sequence used to predict the secondary structure by RNAfold for the illustration in Fig.2B.

### 8.Dimer

The script drawPie.R is used to to visualize the proportion of RBPs with a dimeric binding pattern in the pie chart (Supplemental Fig. S4) .

### 9.MotifMatchesEnrichment

MetaPlot_PerMotif-splicing180419.r. The script to generates the metaplot presented in Fig 6 and Supplemental_Data_S1-S2. The required input can be found onin Supplemental_Data_S4
Density2GO-tab-sp2.r. The script used to obtain GO Enrichment (R package clusterProfiler) using Ensembl gene ID queried through biomaRt package.
MetaPlot_PerMotif-distance2020April05.r. The script to generate the distance to feature histrograms Supplemental Data S3.
