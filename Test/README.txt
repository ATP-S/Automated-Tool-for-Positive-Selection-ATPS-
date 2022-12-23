#sammary
We were able to use our python-based automation tool for positive selection pipeline (ATPS) to perform detailed comparative evolutionary genomic analyses and then estimate the nonsynonymous to synonymous nucleotide substitution ratios.Our tool ATPS can fetch sequences of species from NCBI. In addition, the users can input their data using arguments to scan for positively selected genes based on sequence comparison and a codon-based method. The program is capable of dealing with different software packages that operate together via Python coding. Additionally, the program can handle the input and output files in each step and provides  the user with the results for each gene, which will then be saved in a directory. The ATPS was successfully developed to detect positively selected genes (PSGs), allow the analysis of specific evolutionary branches, and print out the Bayes Empirical Bayes (BEB) analysis, which figures out the positive amino acid sites.




Test N_1
-----------
we have tested ATPS on nine genes in 19 species using this command line 
python3 ATPS.py -G TP53,ATM,CDX2,FOXA2,NF1,NKX2–1,RB1,STK11,APC -S Homo_sapiens,Felis_catus,Pan_troglodytes,Equus_caballus,Canis_lupus,Mesocricetus_auratus,Rattus_norvegicus,Gorilla_gorilla,Sus_scrofa,Tupaia_chinensis,Cavia_porcellus,Heterocephalus_glaber,Bubalus_bubalis,Ovis_aries,Macaca_mulatta,Oryctolagus_cuniculus,Bos_taurus,Cricetulus_griseus,Macaca_fascicularis,Loxodonta_africana -I Homo_sapiens -A mu

# qiime2 (the QIIME 2 framework)

![](https://github.com/qiime2/qiime2/workflows/ci/badge.svg)

Source code repository for the QIIME 2 framework.

QIIME 2™ is a powerful, extensible, and decentralized microbiome bioinformatics
platform that is free, open source, and community developed. With a focus on
data and analysis transparency, QIIME 2 enables researchers to start an
analysis with raw DNA sequence data and finish with publication-quality figures
and statistical results.

Visit [https://qiime2.org](https://qiime2.org) to learn more about the QIIME 2
project.

## Installation

Detailed instructions are available in the
[documentation](https://docs.qiime2.org/).

## Users

Head to the [user docs](https://docs.qiime2.org/) for help getting started,
core concepts, tutorials, and other resources.

Just have a question? Please ask it in our
[forum](https://forum.qiime2.org/c/user-support).

## Developers

Please visit the [contributing page](https://dev.qiime2.org) for more
information on contributions, documentation links, and more.

## Citing QIIME 2
