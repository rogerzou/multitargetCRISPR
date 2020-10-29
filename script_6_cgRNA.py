"""
Script for:
(1) Determination of time-resolved MRE11, 53BP1, and gH2AX enrichment after activation of cgRNA.
"""

import src.mtss as m
import src.chipseq as c
import src.hic as hic
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu - Lab computer)
    desktop = "/mnt/c/Users/Roger/Desktop/"
    hg38 = "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
h2WTin = labhome + "200212_chipseq_WT1/A18_WT1H_final.bam"
bpWTin = labhome + "200212_chipseq_WT1/A15_WT1B_final.bam"
mreGG00m_cg = labhome + "201012_chipseq/A01_final.bam"
mreGG10m_cg = labhome + "201012_chipseq/A02_final.bam"
mreGG30m_cg = labhome + "201012_chipseq/A03_final.bam"
bpGG00m_cg = labhome + "201012_chipseq/A04_final.bam"
bpGG10m_cg = labhome + "201012_chipseq/A14_final.bam"
bpGG30m_cg = labhome + "201012_chipseq/A15_final.bam"
h2GG00m_cg = labhome + "201012_chipseq/A16_final.bam"
h2GG10m_cg = labhome + "201012_chipseq/A17_final.bam"
h2GG30m_cg = labhome + "201012_chipseq/A18_final.bam"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_1_new_peaks.narrowPeak"
GG_mre11 = labhome + "200206_chipseq/macs/16_AluGG-MRE11_peaks.narrowPeak"

""" Set analysis path """
ana = labhome + "Alu_ana_6_cgRNA/"
os.makedirs(ana) if not os.path.exists(ana) else None


""" ############################################################################################ """
""" Get read subsets for vfCRISPR time-resolved MRE11 ChIP-seq for AluGG """
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG00m_cg,
               ana_1 + "GG-C9_cgMRE_00m_rs_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG10m_cg,
               ana_1 + "GG-C9_cgMRE_10m_rs_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG30m_cg,
               ana_1 + "GG-C9_cgMRE_30m_rs_1")


""" ############################################################################################ """
""" Get vfCRISPR AluGG 53BP1 and gH2AX ChIP-seq averaged enrichment as wiggle files. """
ana_2 = ana + "2_wiggle_windows/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
c.to_wiggle_windows(h2GG00m_cg, ana_2 + "GG_cgH2_00m_1", 500)
c.to_wiggle_windows(h2GG10m_cg, ana_2 + "GG_cgH2_10m_1", 500)
c.to_wiggle_windows(h2GG30m_cg, ana_2 + "GG_cgH2_30m_1", 500)
c.to_wiggle_windows(bpGG00m_cg, ana_2 + "GG_cgBP_00m_1", 500)
c.to_wiggle_windows(bpGG10m_cg, ana_2 + "GG_cgBP_10m_1", 500)
c.to_wiggle_windows(bpGG30m_cg, ana_2 + "GG_cgBP_30m_1", 500)


""" ############################################################################################ """
""" Get 53BP1 and gH2AX span and profiles """
ana_3 = ana + "3_span_profiles/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, h2GG10m_cg, h2GG00m_cg, ana_3 + "GG_cgH2_10m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, h2GG30m_cg, h2GG00m_cg, ana_3 + "GG_cgH2_30m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, bpGG10m_cg, bpGG00m_cg, ana_3 + "GG_cgBP_10m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, bpGG30m_cg, bpGG00m_cg, ana_3 + "GG_cgBP_30m_1_span")
