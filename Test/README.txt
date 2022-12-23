#sammary
We were able to use our python-based automation tool for positive selection pipeline (ATPS) to perform detailed comparative evolutionary genomic analyses and then estimate the nonsynonymous to synonymous nucleotide substitution ratios.Our tool ATPS can fetch sequences of species from NCBI. In addition, the users can input their data using arguments to scan for positively selected genes based on sequence comparison and a codon-based method. The program is capable of dealing with different software packages that operate together via Python coding. Additionally, the program can handle the input and output files in each step and provides  the user with the results for each gene, which will then be saved in a directory. The ATPS was successfully developed to detect positively selected genes (PSGs), allow the analysis of specific evolutionary branches, and print out the Bayes Empirical Bayes (BEB) analysis, which figures out the positive amino acid sites.




Test N_1
-----------
we have tested ATPS on nine genes in 19 species using this command line 
python3 ATPS.py -G TP53,ATM,CDX2,FOXA2,NF1,NKX2–1,RB1,STK11,APC -S Homo_sapiens,Felis_catus,Pan_troglodytes,Equus_caballus,Canis_lupus,Mesocricetus_auratus,Rattus_norvegicus,Gorilla_gorilla,Sus_scrofa,Tupaia_chinensis,Cavia_porcellus,Heterocephalus_glaber,Bubalus_bubalis,Ovis_aries,Macaca_mulatta,Oryctolagus_cuniculus,Bos_taurus,Cricetulus_griseus,Macaca_fascicularis,Loxodonta_africana -I Homo_sapiens -A mu


# Summary

PosiGene is a tool that (i) detects positively selected genes on genome-scale, 
(ii) allows analysis of specific evolutionary branches, (iii) can be used in 
arbitrary species contexts and (iv) offers visualization of the candidates. As 
data input the program requires only the coding sequences of your chosen species
set in fasta or genbank format. From them, orthologs, alignments and a 
phylogenetic tree are reconstructed to finally apply the branch-site test of 
positive selection. Filtering mechanisms are implemented to minimize the 
occurrence of false positives. PosiGene was tested on simulated as well as real
data to ensure the reliability of the predicted positively selected genes.

# Installation

After unpacking, no further installation steps are needed.

# Documentation

To learn how to use PosiGene please read the user guide that can be found under 
doc/user_guide.pdf.

# Test run

To test whether the package works please execute:
perl PosiGene.pl -o=test  -as=Harpegnathos_saltator  -tn=10  -rs=Acromyrmex_echinatior:test_data/Acromyrmex_echinatior_sample.fasta  -nhsbr=Acromyrmex_echinatior:test_data/Acromyrmex_echinatior_sample.fasta,Atta_cephalotes:test_data/Atta_cephalotes_sample.fasta,Camponotus_floridanus:test_data/Camponotus_floridanus_sample.fasta,Harpegnathos_saltator:test_data/Harpegnathos_saltator_sample.fasta,Linepithema_humile:test_data/Linepithema_humile_sample.fasta,Pogonomyrmex_barbatus:test_data/Pogonomyrmex_barbatus_sample.fasta,Solenopsis_invicta:test_data/Solenopsis_invicta_sample.fasta

It should be finished after some minutes, telling you where you can find a 
result table. If the program runs through and the produced result table equals 
that at test_data/Harpegnathos_saltator_results_short.tsv everything is fine.
