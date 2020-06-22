#Information content must be calculated from tree parts. 
#Note: The SLM model generates frequency table, while the normal pwm is the count table. 

perl IC_calculator17Apr29.pl --m f ../dominating17Apr28/slm_replacement 
perl IC_calculator17Apr29.pl --m d ../dominating17Apr28/dimatrix
perl IC_calculator17Apr29.pl --m n ../dominating17Apr28/prime_motifs


#the following part is for the stem IC correlation  (Supplemetal Figs)
perl  stem_loop_IC-parser.pl --f 0 Bits_d_2017May24.table fordi_slm.length.table > Bits_d_2017May24.table.labeled
perl stem_loop_IC-parser.pl --f 1 Bits_f_2017May24.table forpwm_slm.length.table > Bits_f_2017May24.table.labeled
Rscript --vanilla StemIC_corrPlot.r Bits_f_2017May24.table.labeled Bits_d_2017May24.table.labeled

#the following section is for the STEM and LOOP comparison, also overall landscape

perl stem_loop_IC-parser.pl --f 0 Bits_f_2017May24.table forpwm_slm.length.table  > Bits_f_2017May24.table.SLM_stemloop


cat Bits_n_2017May24.table Bits_f_2017May24.table > All_Bits2017May24.table
Rscript --vanilla draw_overall_perbase.R All_Bits2017May24.table pubavail.bits.table 0.05