Analysis software for multi-target CRISPR
====
Companion code for:

*Zou, R.S., Marin-Gonzalez, A., Liu, Y., Liu, H.B., Shen, L., Dveirin, R., Luo, J.X., Kalhor, R. and Ha, T. Massively parallel genomic perturbations with multi-target CRISPR reveal new insights on Cas9 activity and DNA damage responses at endogenous sites. bioRxiv (2022). https://www.biorxiv.org/content/10.1101/2022.01.18.476836v1*

## Software requirements
- [Anaconda Python 3.7](https://www.anaconda.com/distribution/) (Anaconda's python distribution comes with the required numpy and scipy libraries)
- [pysam](https://pysam.readthedocs.io/en/latest/installation.html)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://www.htslib.org/download/)
- Ensure that both `samtools` and `bowtie2` are added to path and can be called directly from bash

## Data requirements
- https://www.ncbi.nlm.nih.gov/bioproject/PRJNA733683

## Installation
1. Download sequencing reads in FASTQ format from SRA
2. Download the prebuilt bowtie2 indices for human hg19 and hg38 genome assemblies
    - [Human hg38](https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip)
    - [Human hg19](https://genome-idx.s3.amazonaws.com/bt/hg19.zip)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/` and `hg19_bowtie2/`
3. Download two human hg19 and hg38 genome assemblies in FASTA format
    - [hg38.fa](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
    - [hg19.fa](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
    - Extract from archive, move to the corresponding folders named `hg38_bowtie2/` and `hg19_bowtie2/`
4. Generate FASTA file indices
    - `samtools faidx hg38_bowtie2/hg38.fa`
    - `samtools faidx hg19_bowtie2/hg19.fa`

## Usage
Bash scripts are used to automate the processing of sequencing data.

Python scripts are used to perform analysis of various data featured in the manuscript.
They are labeled `script_*_*.py`, such as `script_1_putative.py`


## List of mgRNA sequences
List of target sites for the 'GG', 'CT', and 'TA' mgRNA sequences, along with its 
hg38 genomic context, are included in  `mgRNA_target_sites.xlsx`.