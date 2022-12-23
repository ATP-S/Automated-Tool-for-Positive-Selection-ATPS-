We were able to use our python-based automation tool for positive selection pipeline (ATPS) to perform detailed comparative evolutionary genomic analyses and then estimate the nonsynonymous to synonymous nucleotide substitution ratios.Our tool ATPS can fetch sequences of species from NCBI. In addition, the users can input their data using arguments to scan for positively selected genes based on sequence comparison and a codon-based method. The program is capable of dealing with different software packages that operate together via Python coding. Additionally, the program can handle the input and output files in each step and provides  the user with the results for each gene, which will then be saved in a directory. The ATPS was successfully developed to detect positively selected genes (PSGs), allow the analysis of specific evolutionary branches, and print out the Bayes Empirical Bayes (BEB) analysis, which figures out the positive amino acid sites.




Test N_1
-----------
we have tested ATPS on nine genes in 19 species using this command line 
python3 ATPS.py -G TP53,ATM,CDX2,FOXA2,NF1,NKX2–1,RB1,STK11,APC -S Homo_sapiens,Felis_catus,Pan_troglodytes,Equus_caballus,Canis_lupus,Mesocricetus_auratus,Rattus_norvegicus,Gorilla_gorilla,Sus_scrofa,Tupaia_chinensis,Cavia_porcellus,Heterocephalus_glaber,Bubalus_bubalis,Ovis_aries,Macaca_mulatta,Oryctolagus_cuniculus,Bos_taurus,Cricetulus_griseus,Macaca_fascicularis,Loxodonta_africana -I Homo_sapiens -A mu
