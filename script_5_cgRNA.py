"""
Script for:
(1) Determination of time-resolved MRE11 enrichment after activation of cgRNA.
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
h2WTin = datadir + "200212_chipseq_WT1/A18_gh2ax_hg38_final.bam"
bpWTin = datadir + "200212_chipseq_WT1/A15_53bp1_hg38_final.bam"
mreGG00m_cg = datadir + "201012_chipseq/A01_hg38_final.bam"
mreGG10m_cg = datadir + "201012_chipseq/A02_hg38_final.bam"
mreGG30m_cg = datadir + "201012_chipseq/A03_hg38_final.bam"
bpGG00m_cg = datadir + "201012_chipseq/A04_hg38_final.bam"
bpGG10m_cg = datadir + "201012_chipseq/A14_hg38_final.bam"
bpGG30m_cg = datadir + "201012_chipseq/A15_hg38_final.bam"
h2GG00m_cg = datadir + "201012_chipseq/A16_hg38_final.bam"
h2GG10m_cg = datadir + "201012_chipseq/A17_hg38_final.bam"
h2GG30m_cg = datadir + "201012_chipseq/A18_hg38_final.bam"
alnpath_hg38 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_Alu_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = datadir + "Alu_ana_5_cgRNA/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_peak_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Get read subsets for vfCRISPR time-resolved MRE11 ChIP-seq for AluGG """
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, mreGG00m_cg, ana_1 + "GG_cgMRE_00m_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, mreGG10m_cg, ana_1 + "GG_cgMRE_10m_rs_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8),
               hg38, mreGG30m_cg, ana_1 + "GG_cgMRE_30m_rs_1")
subset_mre11_1 = [ana_1 + "GG_cgMRE_00m_rs_1.csv",
                  ana_1 + "GG_cgMRE_10m_rs_1.csv",
                  ana_1 + "GG_cgMRE_30m_rs_1.csv"]
m.read_kinetics(subset_mre11_1, ana_1 + "GG_cgMRE_kin_1", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_mre11_1, ana_1 + "GG_cgMRE_kin_1", endname='RefSeq sense', hname='Cspan')
m.read_kinetics(subset_mre11_1, ana_1 + "GG_cgMRE_kin_1", endname='RefSeq sense', hname='Cend')


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq. """
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38[0], 1500, AluGG),
                             mreGG00m_cg, ana_2 + "GG_cgON_mre11_00m")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38[0], 1500, AluGG),
                             mreGG10m_cg, ana_2 + "GG_cgON_mre11_10m")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38[0], 1500, AluGG),
                             mreGG30m_cg, ana_2 + "GG_cgON_mre11_30m")
