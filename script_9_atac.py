"""
Script for analysis of ATAC-seq results
"""

import src.mtss as m
import src.msa as msa
import src.hic as hic
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu - Lab computer)
    desktop = "/mnt/c/Users/Roger/Desktop/"
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
atacWT = labhome + "201110_atac/N03_sorted.bam"
atacGG3h = labhome + "201110_atac/N02_sorted.bam"
alnpath = labhome + "Alu_ana_1_putative/1_protosearch/psearch_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
GG_mre11 = labhome + "200206_chipseq/macs/AluGG-MRE11_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = labhome + "Alu_ana_9_atac/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Obtain read counts in kb scale to mb scale for ATAC-seq """
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacWT, ana_1 + "WT-ON_atac_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacGG3h, ana_1 + "GG-ON_atac_3h_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              atacWT, ana_1 + "WT-ON_atac_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              atacGG3h, ana_1 + "GG-ON_atac_3h_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5E4, AluGG),
              atacWT, ana_1 + "WT-ON_atac_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5E4, AluGG),
              atacGG3h, ana_1 + "GG-ON_atac_3h_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5E5, AluGG),
              atacWT, ana_1 + "WT-ON_atac_500000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5E5, AluGG),
              atacGG3h, ana_1 + "GG-ON_atac_3h_500000_rc.csv")


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites """
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG),
                    hg38, atacWT, ana_2 + "WT_atac", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG),
                    hg38, atacGG3h, ana_2 + "GG-3h_atac", span_rad=1500, res=1, wind_rad=2)
hic.get_span_width(msa.target_gen(alnpath, hg38, 1500, AluGG),
                   hg38, atacGG3h, atacWT, ana_2 + "GG-ON-atac_3h", w_rad=20, skip=1, false_ct=10)
