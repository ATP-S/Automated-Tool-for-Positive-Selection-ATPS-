# Introduction

We were able to use our python-based automation tool for positive selection pipeline (ATPS) to perform detailed comparative evolutionary genomic analyses and then estimate the nonsynonymous to synonymous nucleotide substitution ratios.Our tool ATPS can fetch sequences of species from NCBI. In addition, the users can input their data using arguments to scan for positively selected genes based on sequence comparison and a codon-based method. The program is capable of dealing with different software packages that operate together via Python coding. Additionally, the program can handle the input and output files in each step and provides  the user with the results for each gene, which will then be saved in a directory. The ATPS was successfully developed to detect positively selected genes (PSGs), allow the analysis of specific evolutionary branches, and print out the Bayes Empirical Bayes (BEB) analysis, which figures out the positive amino acid sites.



# Wiki (How to use the program)

A detailed documentation including: How to use the tool, its arguments, and a detailed description of what's included can be found at the [GitHub Wiki](https://github.com/APS-P/APSP/wiki). In short, it contains a detailed guide on how to use ATPS

   # Who is it for?

The package was designed with the sole aim of facilitating positive selection analysis on geneticists, whether the user has specific genes in mind that they would like to see if it's under positive selection, or they have a generated file that they want to use.

Another goal that came along the way, was to combine all popular tools in a positive selection pipeline, efficiently under one functional code, saving time for users and organizing all the expected outputs in the evolutionary analysis pipeline under one folder. Users no longer need to manually insert output files from each independent tool, now they either let it automatically fetch and work, or give it the generated file of their choice and let the automation do the job.

   # How it works?

The flowchart below explains briefly how the program works:
![Image for flowchart](https://github.com/APS-P/APSP/blob/main/images/work_flow.jpeg)

- In the figure, we start with filling arguments from the table below depending on the preferred pipeline, the user chooses specie/s, gene/s, target specie/s, and the preferred MSA tool. The user could also state that they have a ready input file, in such case the fetch step will be skipped and it will start directly at the MSA step. 

- In the figure, you will also notice the word "Parallel" these are two Bio-python packages responsible for the calculation of MLD file and calculation of PhloP score consecutively and in parallel with the other steps in the diagram.

- The program does all of this on every single gene.

- Note: Multiple sequence alignment (MSA) uses one of three most used software aligners (MAFFT, MUSCLE, ClustalO)

   # Input

* To input your own **.FASTA** file manually, you need it in a specific format to work, we will provide you with it in this >link<.
* If you prefer to conduct your research using a fetched input, then the program will fetch the most accurate, highest-scoring, and longest(Not spliced) isoform of your gene.

   # Output (Per gene)

- Each package and tool within our program results in its own output.
- CodeML spreadsheet output for each model.
- wigscore spreadsheet, and a picture (.png) of it.
- BEB spreadsheet

   # Arguments

Check the table below:
- Note: Arguments are case sensitive, only write in CAPSLOCK mode.

| Argument  |  Function |
|---|---|
| -G | Enter Gene or list of genes separated by a comma (,) |
|  -S |  Pass specie name or a list of species |
| -I  |  Target species for selection |
| -T  |  Thread count (Default is full utilization)|
|  -A |  Type of alignment tool (MAFFT = mf, MUSCLE = mu, ClustalO = cl)  |
|  -IF |  The manual input if the user chooses to, must be in the format we provided above |
| - O | output directory 


# Installation instructions

These instructions will guide you through the process of installing and running the [marwanjs/atps].

   # Prerequisites

Before you begin, you will need to have [Docker](https://www.docker.com/get-started/) installed on your system.
   
   # Installation

1. Open a terminal window.
2. Run the following command to pull the Docker image from Docker Hub:

```console
$ docker pull marwanjs/atps:latest
```

3. Once the image has been downloaded, run the following command to start the container:


```console
$ docker run -it -v $HOME:$HOME c4dae898af64
```

Optional- for ubuntu users. to custumize the CPUs and RAM inside the docker image, run the following command:

```console
$ docker run --memory=4g --cpus=8 -it -v $HOME:$HOME c4dae898af64
```

   # Example for a test run
   
```console
$ python3 ATPS.py -G TP53,ATM,CDX2,FOXA2,NF1 -S Homo_sapiens,Felis_catus,Rattus_norvegicus -I Homo_sapiens -A mu -R 100 -O home/Users/Desktop
```

# Dependencies

   # softwares 
   
1- MUSCLE

2- MAFFT

3- clustalo

4- Gblocks

8- jmodeltest

9- phyml

8- phylofit

9- codeml

   # python packages 
   
1- biopython

2- pandas

3- glob

4- matplotlib 

5- os 

6- bokeh

7- mne

8- numpy

9- gc 

10- tkinter 


# Technical support

ATPS is actively supported by our team. Please use the [Issues section](https://github.com/APS-P/APSP/issues) if you happen to find any errors.

## Citation
Will be added after publication and official release.
