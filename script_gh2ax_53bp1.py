"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.mtss as m
import src.hic as hic
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
enc, enc_a = datadir + "public/", datadir + "public/analysis/"
h2WTin = datadir + "200212_chipseq_WT1/A18_gh2ax_hg38_final.bam"
bpWTin = datadir + "200212_chipseq_WT1/A15_53bp1_hg38_final.bam"
h2GGin = datadir + "200206_chipseq/AluGG-gH2AX_hg38_final.bam"
bpGGin = datadir + "200206_chipseq/AluGG-53BP1_hg38_final.bam"
h2TAin = datadir + "200316_chipseq/AluTA-gh2ax-rep1_hg38_final.bam"
bpTAin = datadir + "200316_chipseq/AluTA-53bp1-rep1_hg38_final.bam"
h2CTin = datadir + "200316_chipseq/AluCT-gh2ax-rep1_hg38_final.bam"
bpCTin = datadir + "200316_chipseq/AluCT-53bp1-rep1_hg38_final.bam"
bpGGin_noD = datadir + "200316_chipseq/AluGG-53bp1-noD-rep1_hg38_final.bam"
bpGGin_PKi = datadir + "200316_chipseq/AluGG-53bp1-PKi-rep1_hg38_final.bam"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

""" macs2 output """
# for on-target analysis at GG, CT, TA
GG_cas9 = datadir + "200206_chipseq/macs/AluGG-Cas9_hg38_final_peaks.narrowPeak"
CT_cas9 = datadir + "200316_chipseq/macs/AluCT-cas9-rep1_hg38_final_peaks.narrowPeak"
TA_cas9 = datadir + "200316_chipseq/macs/AluTA-cas9-rep1_hg38_final_peaks.narrowPeak"
GG_mre11 = datadir + "200206_chipseq/macs/AluGG-MRE11_hg38_final_peaks.narrowPeak"
CT_mre11 = datadir + "200316_chipseq/macs/AluCT-mre11-rep1_hg38_final_peaks.narrowPeak"
TA_mre11 = datadir + "200316_chipseq/macs/AluTA-mre11-rep1_hg38_final_peaks.narrowPeak"
# for on- and off-target analysis at GG
GG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
GG3h_npk2 = datadir + "200804_chipseq/macs/A09_hg38_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = datadir + "Alu_ana_4_gh2ax_53bp1/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_counts/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_wiggle_windows/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_filt_count_span/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_filt_epigenetics/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
ana_5 = ana + "5_merged/"
os.makedirs(ana_5) if not os.path.exists(ana_5) else None


# """ ############################################################################################ """
# """ Get all enrichment counts in 10kb window for AluGG, AluTA, and AluCT """
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), h2GGin, ana_1 + "GG-C9_gh2ax_10000_rc.csv")
# m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), h2CTin, ana_1 + "CT-C9_gh2ax_10000_rc.csv")
# m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), h2TAin, ana_1 + "TA-C9_gh2ax_10000_rc.csv")
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin, ana_1 + "GG-C9_53bp1_10000_rc.csv")
# m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), bpCTin, ana_1 + "CT-C9_53bp1_10000_rc.csv")
# m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), bpTAin, ana_1 + "TA-C9_53bp1_10000_rc.csv")
#
#
# """ ############################################################################################ """
# """ Obtain peak profiles of gH2AX and 53BP1 data """
# c.to_wiggle_windows(hg38, h2GGin, ana_2 + "GG_gh2ax", 500)
# c.to_wiggle_windows(hg38, bpGGin, ana_2 + "GG_53bp1", 500)
# c.to_wiggle_windows(hg38, h2CTin, ana_2 + "CT_gh2ax", 500)
# c.to_wiggle_windows(hg38, bpCTin, ana_2 + "CT_53bp1", 500)
# c.to_wiggle_windows(hg38, h2TAin, ana_2 + "TA_gh2ax", 500)
# c.to_wiggle_windows(hg38, bpTAin, ana_2 + "TA_53bp1", 500)
# c.to_wiggle_windows(hg38, h2WTin, ana_2 + "WT_gh2ax", 500)
# c.to_wiggle_windows(hg38, bpWTin, ana_2 + "WT_53bp1", 500)
#
#
# """ ############################################################################################ """
# """ For target sites separated by more than 2MB, determine 53BP1 and gH2AX span width and counts """
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
# hic.get_span_width(gen, hg38, h2GGin, h2WTin, ana_3 + "GG-Mfilt_gh2ax_width")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
# hic.get_span_width(gen, hg38, h2CTin, h2WTin, ana_3 + "CT-Mfilt_gh2ax_width")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
# hic.get_span_width(gen, hg38, h2TAin, h2WTin, ana_3 + "TA-Mfilt_gh2ax_width")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
# hic.get_span_width(gen, hg38, bpGGin, bpWTin, ana_3 + "GG-Mfilt_53bp1_width")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
# hic.get_span_width(gen, hg38, bpCTin, bpWTin, ana_3 + "CT-Mfilt_53bp1_width")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
# hic.get_span_width(gen, hg38, bpTAin, bpWTin, ana_3 + "TA-Mfilt_53bp1_width")
#
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 10000, hg38, AluGG), 2000000)
# m.read_counts(gen, h2GGin, ana_3 + "GG-Mfilt_gh2ax_10000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 10000, hg38, AluCT), 2000000)
# m.read_counts(gen, h2CTin, ana_3 + "CT-Mfilt_gh2ax_10000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 10000, hg38, AluTA), 2000000)
# m.read_counts(gen, h2TAin, ana_3 + "TA-Mfilt_gh2ax_10000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 10000, hg38, AluGG), 2000000)
# m.read_counts(gen, bpGGin, ana_3 + "GG-Mfilt_53bp1_10000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 10000, hg38, AluCT), 2000000)
# m.read_counts(gen, bpCTin, ana_3 + "CT-Mfilt_53bp1_10000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 10000, hg38, AluTA), 2000000)
# m.read_counts(gen, bpTAin, ana_3 + "TA-Mfilt_53bp1_10000_rc.csv")
#
#
# """ ############################################################################################ """
# """ Get enrichment counts and epigenetic, mismatch, and chromatin information, FILTERED for
#     MRE11 enrichment sites separated by 2mb for AluGG, AluTA, and AluCT """
# # AluGG datasets
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, h3k4me1_1, ana_4 + "GG-Mfilt_h3k4me1_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, h3k4me3_1, ana_4 + "GG-Mfilt_h3k4me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, h3k9me3_1, ana_4 + "GG-Mfilt_h3k9me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, h3k27ac_1, ana_4 + "GG-Mfilt_h3k27ac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, h3k36me3_1, ana_4 + "GG-Mfilt_h3k36me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50, hg38, AluGG), 2000000)
# m.read_counts(gen, dnasei_1, ana_4 + "GG-Mfilt_dnasei_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50, hg38, AluGG), 2000000)
# m.read_counts(gen, mnase_1, ana_4 + "GG-Mfilt_mnase_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, atac_1, ana_4 + "GG-Mfilt_atac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 50000, hg38, AluGG), 2000000)
# m.read_counts(gen, rna_3, ana_4 + "GG-Mfilt_rna_50000_rc.csv")
# # AluCT datasets
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, h3k4me1_1, ana_4 + "CT-Mfilt_h3k4me1_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, h3k4me3_1, ana_4 + "CT-Mfilt_h3k4me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, h3k9me3_1, ana_4 + "CT-Mfilt_h3k9me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, h3k27ac_1, ana_4 + "CT-Mfilt_h3k27ac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, h3k36me3_1, ana_4 + "CT-Mfilt_h3k36me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50, hg38, AluCT), 2000000)
# m.read_counts(gen, dnasei_1, ana_4 + "CT-Mfilt_dnasei_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50, hg38, AluCT), 2000000)
# m.read_counts(gen, mnase_1, ana_4 + "CT-Mfilt_mnase_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, atac_1, ana_4 + "CT-Mfilt_atac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 50000, hg38, AluCT), 2000000)
# m.read_counts(gen, rna_3, ana_4 + "CT-Mfilt_rna_50000_rc.csv")
# # AluTA datasets
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, h3k4me1_1, ana_4 + "TA-Mfilt_h3k4me1_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, h3k4me3_1, ana_4 + "TA-Mfilt_h3k4me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, h3k9me3_1, ana_4 + "TA-Mfilt_h3k9me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, h3k27ac_1, ana_4 + "TA-Mfilt_h3k27ac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, h3k36me3_1, ana_4 + "TA-Mfilt_h3k36me3_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50, hg38, AluTA), 2000000)
# m.read_counts(gen, dnasei_1, ana_4 + "TA-Mfilt_dnasei_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50, hg38, AluTA), 2000000)
# m.read_counts(gen, mnase_1, ana_4 + "TA-Mfilt_mnase_50_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, atac_1, ana_4 + "TA-Mfilt_atac_50000_rc.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 50000, hg38, AluTA), 2000000)
# m.read_counts(gen, rna_3, ana_4 + "TA-Mfilt_rna_50000_rc.csv")
# # One-mismatch annotation
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
# m.read_mismatch(gen, ana_4 + "GG-Mfilt_mismatch.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
# m.read_mismatch(gen, ana_4 + "CT-Mfilt_mismatch.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
# m.read_mismatch(gen, ana_4 + "TA-Mfilt_mismatch.csv")
# # ChromHMM epigenetic chromatin state annotation
# gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
# m.read_chromhmm(gen, hg38, ana_4 + "GG-Mfilt_chromhmm.csv")
# gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
# m.read_chromhmm(gen, hg38, ana_4 + "CT-Mfilt_chromhmm.csv")
# gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
# m.read_chromhmm(gen, hg38, ana_4 + "TA-Mfilt_chromhmm.csv")
#
# B_GG = [m.load_nparray(ana_3 + "GG-Mfilt_gh2ax_10000_rc.csv"),
#         m.load_nparray(ana_3 + "GG-Mfilt_53bp1_width.csv"),
#         m.load_nparray(ana_3 + "GG-Mfilt_gh2ax_width.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_h3k4me1_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_h3k4me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_h3k9me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_h3k27ac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_h3k36me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_dnasei_50_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_mnase_50_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_atac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_rna_50000_rc.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_chromhmm.csv"),
#         m.load_nparray(ana_4 + "GG-Mfilt_mismatch.csv")]
# B_CT = [m.load_nparray(ana_3 + "CT-Mfilt_gh2ax_10000_rc.csv"),
#         m.load_nparray(ana_3 + "CT-Mfilt_53bp1_width.csv"),
#         m.load_nparray(ana_3 + "CT-Mfilt_gh2ax_width.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_h3k4me1_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_h3k4me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_h3k9me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_h3k27ac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_h3k36me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_dnasei_50_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_mnase_50_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_atac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_rna_50000_rc.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_chromhmm.csv"),
#         m.load_nparray(ana_4 + "CT-Mfilt_mismatch.csv")]
# B_TA = [m.load_nparray(ana_3 + "TA-Mfilt_gh2ax_10000_rc.csv"),
#         m.load_nparray(ana_3 + "TA-Mfilt_53bp1_width.csv"),
#         m.load_nparray(ana_3 + "TA-Mfilt_gh2ax_width.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_h3k4me1_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_h3k4me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_h3k9me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_h3k27ac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_h3k36me3_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_dnasei_50_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_mnase_50_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_atac_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_rna_50000_rc.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_chromhmm.csv"),
#         m.load_nparray(ana_4 + "TA-Mfilt_mismatch.csv")]
# num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
# head = "region string, coordinate, sense, sequence, mismatches, " \
#           "53bp1_rc, gh2ax_rc, 53bp1_w, gh2ax_w, " \
#           "h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, " \
#           "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"
# """ generate merged datasets """
# GG = m.load_nparray(ana_3 + "GG-Mfilt_53bp1_10000_rc.csv")
# GG_m = m.mergesubsetcounts(GG, B_GG, num_cols, ana_5 + "GG-Mfilt_bp-h2_merged.csv", head)
# CT = m.load_nparray(ana_3 + "CT-Mfilt_53bp1_10000_rc.csv")
# CT_m = m.mergesubsetcounts(CT, B_CT, num_cols, ana_5 + "CT-Mfilt_bp-h2_merged.csv", head)
# TA = m.load_nparray(ana_3 + "TA-Mfilt_53bp1_10000_rc.csv")
# TA_m = m.mergesubsetcounts(TA, B_TA, num_cols, ana_5 + "TA-Mfilt_bp-h2_merged.csv", head)
# m.mergerows([GG_m, CT_m, TA_m], ana_5 + "ALL-Mfilt_bp-h2_merged.csv", head)
