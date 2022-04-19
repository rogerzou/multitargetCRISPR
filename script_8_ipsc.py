"""
Script for analysis of Cas9 and MRE11 ChIP-seq in iPSCs
"""

import src.mtss as m
import src.msa as msa
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
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

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
casGG3h_npk2 = datadir + "200804_chipseq/macs/A09_hg38_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" Set analysis path """
ana = datadir + "Alu_ana_8_ipsc/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets_counts/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_peak_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" For all putative on-target sites, determine paired-end read subsets for Cas9 and MRE11,
    for both HEK and iPSC datasets. """
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCmreGGr1, ana_1 + "target_iPSCmreGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCmreGGr2, ana_1 + "target_iPSCmreGGr2_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCcasGGr1, ana_1 + "target_iPSCcasGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               iPSCcasGGr2, ana_1 + "target_iPSCcasGGr2_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               HEKmreGGr1, ana_1 + "target_HEKmreGGr1_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               HEKcasGGr1, ana_1 + "target_HEKcasGGr1_rs")


""" ############################################################################################ """
""" For all putative target sites, determine paired-end read subsets for Cas9 and MRE11,
    for both HEK and iPSC datasets. """
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               iPSCmreGGr1, ana_1 + "macs_iPSCmreGGr1_rs")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               iPSCmreGGr2, ana_1 + "macs_iPSCmreGGr2_rs")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               iPSCcasGGr1, ana_1 + "macs_iPSCcasGGr1_rs")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               iPSCcasGGr2, ana_1 + "macs_iPSCcasGGr2_rs")


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq. """
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             iPSCcasGGr1, ana_2 + "iPSCcasGGr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             iPSCcasGGr2, ana_2 + "iPSCcasGGr2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             iPSCmreGGr1, ana_2 + "iPSCmreGGr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             iPSCmreGGr2, ana_2 + "iPSCmreGGr2")
