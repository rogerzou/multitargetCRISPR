"""
Script for analysis of Cas9 and MRE11 ChIP-seq in iPSCs
"""

import src.mtss as m
import src.msa as msa
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """

iPSCcasGGr1 = datadir + "210301_chipseq/A01_hg38_final.bam"
iPSCcasGGr2 = datadir + "210301_chipseq/A02_hg38_final.bam"
iPSCmreGGr1 = datadir + "210301_chipseq/A03_hg38_final.bam"
iPSCmreGGr2 = datadir + "210301_chipseq/A04_hg38_final.bam"
HEKmreGGr1 = datadir + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
HEKcasGGr1 = datadir + "200206_chipseq/AluGG-Cas9_hg38_final.bam"
alnpath_hg38 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"

""" macs2 peak detection """
casGGnpk = datadir + "200206_chipseq/macs/AluGG-Cas9_hg38_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" Set analysis path """
ana = datadir + "Alu_ana_8_ipsc/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets_counts/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None


""" ############################################################################################ """
""" For all putative on-target sites, determine paired-end read subsets for Cas9 and MRE11,
    for both HEK and iPSC datasets. """
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCmreGGr1, ana_1 + "iPSCmreGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCmreGGr2, ana_1 + "iPSCmreGGr2_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCcasGGr1, ana_1 + "iPSCcasGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCcasGGr2, ana_1 + "iPSCcasGGr2_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               HEKmreGGr1, ana_1 + "HEKmreGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               HEKcasGGr1, ana_1 + "HEKcasGGr1_rs")
