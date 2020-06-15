
#folder prime_motifs contains all motifs in PFMs format.
perl count_freq.pl ../dominating17Apr28/prime_motifs/|awk '$3 ~/mono/' > mono-20180425-linear.table
perl count_freq.pl ../dominating17Apr28/prime_motifs/|awk '$3 ~/di/' > di-20180425-linear.table
Rscript --vanilla code.r di-20180425-linear.table
Rscript --vanilla code.r mono-20180425-linear.table
