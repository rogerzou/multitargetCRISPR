"""
Script for analysis of BLISS results
"""

import src.mtss as m
import src.msa as msa
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
blissWT = datadir + "201022_bliss/A36_mod_sorted.bam"
blissGG3h = datadir + "201022_bliss/A37_mod_sorted.bam"
blissGG00m_1 = datadir + "201022_bliss/A38_mod_sorted.bam"
blissGG10m_1 = datadir + "201022_bliss/A39_mod_sorted.bam"
blissGG30m_1 = datadir + "201022_bliss/A40_mod_sorted.bam"
blissGG00m_2 = datadir + "201022_bliss/A41_mod_sorted.bam"
blissGG10m_2 = datadir + "201022_bliss/A42_mod_sorted.bam"
blissGG30m_2 = datadir + "201022_bliss/A43_mod_sorted.bam"
alnpath = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
casGG3h_npk2 = datadir + "200804_chipseq/macs/A09_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = datadir + "Alu_ana_6_bliss/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Get read subsets for vfCRISPR time-resolved BLISS for AluGG """
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, blissWT, ana_1 + "WT_bliss_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG3h, ana_1 + "GG_bliss_3h_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG00m_1, ana_1 + "GG_bliss_00m_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG10m_1, ana_1 + "GG_bliss_10m_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG30m_1, ana_1 + "GG_bliss_30m_rs_1")
subset_bliss_1 = [ana_1 + "WT_bliss_rs_1.csv",
                  ana_1 + "GG_bliss_00m_rs_1.csv",
                  ana_1 + "GG_bliss_10m_rs_1.csv",
                  ana_1 + "GG_bliss_30m_rs_1.csv",
                  ana_1 + "GG_bliss_3h_rs_1.csv"]
m.read_kinetics(subset_bliss_1, ana_1 + "GG_bliss_kin_1", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_bliss_1, ana_1 + "GG_bliss_kin_1", endname='RefSeq sense', hname='Cend')


m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8),
               hg38, blissWT, ana_1 + "WT_bliss_rs_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG3h, ana_1 + "GG_bliss_3h_rs_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG00m_2, ana_1 + "GG_bliss_00m_rs_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG10m_2, ana_1 + "GG_bliss_10m_rs_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8),
               hg38, blissGG30m_2, ana_1 + "GG_bliss_30m_rs_2")
subset_bliss_2 = [ana_1 + "WT_bliss_rs_2.csv",
                  ana_1 + "GG_bliss_00m_rs_2.csv",
                  ana_1 + "GG_bliss_10m_rs_2.csv",
                  ana_1 + "GG_bliss_30m_rs_2.csv",
                  ana_1 + "GG_bliss_3h_rs_2.csv"]
m.read_kinetics(subset_bliss_2, ana_1 + "GG_bliss_kin_2", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_bliss_2, ana_1 + "GG_bliss_kin_2", endname='RefSeq sense', hname='Cend')

""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from BLISS """
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissWT, ana_2 + "WT-ON_bliss")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG3h, ana_2 + "GG-ON_bliss_3h")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG00m_1, ana_2 + "GG-ON_bliss_00m_1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG10m_1, ana_2 + "GG-ON_bliss_10m_1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG30m_1, ana_2 + "GG-ON_bliss_30m_1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG00m_2, ana_2 + "GG-ON_bliss_00m_2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG10m_2, ana_2 + "GG-ON_bliss_10m_2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38[0], 1500, AluGG),
                             blissGG30m_2, ana_2 + "GG-ON_bliss_30m_2")
