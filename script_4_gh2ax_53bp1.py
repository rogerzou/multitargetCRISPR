"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.mtss as m
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
enc, enc_a = labhome + "public/", labhome + "public/analysis/"
h2WTin = labhome + "200212_chipseq_WT1/A18_WT1H_final.bam"
bpWTin = labhome + "200212_chipseq_WT1/A15_WT1B_final.bam"
h2GGin = labhome + "200206_chipseq/13_AluGG-gH2AX_final.bam"
bpGGin = labhome + "200206_chipseq/14_AluGG-53BP1_final.bam"
h2TAin = labhome + "200316_chipseq/AluTA-gh2ax-rep1_rmdup.bam"
bpTAin = labhome + "200316_chipseq/AluTA-53bp1-rep1_rmdup.bam"
h2CTin = labhome + "200316_chipseq/AluCT-gh2ax-rep1_rmdup.bam"
bpCTin = labhome + "200316_chipseq/AluCT-53bp1-rep1_rmdup.bam"
bpGGin_noD = labhome + "200316_chipseq/AluGG-53bp1-noD-rep1_rmdup.bam"
bpGGin_PKi = labhome + "200316_chipseq/AluGG-53bp1-PKi-rep1_rmdup.bam"
GG_cas9 = labhome + "200206_chipseq/macs/15_AluGG-Cas9_peaks.narrowPeak"
CT_cas9 = labhome + "200316_chipseq/macs/AluCT-Cas9-rep1_peaks.narrowPeak"
TA_cas9 = labhome + "200316_chipseq/macs/AluTA-Cas9-rep1_peaks.narrowPeak"
GG_mre11 = labhome + "200206_chipseq/macs/16_AluGG-MRE11_peaks.narrowPeak"
CT_mre11 = labhome + "200316_chipseq/macs/AluCT-MRE11-rep1_peaks.narrowPeak"
TA_mre11 = labhome + "200316_chipseq/macs/AluTA-MRE11-rep1_peaks.narrowPeak"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
# rna_1 = enc + "RNAseq_HEK293_SRR1264355.bam"              # RNA-seq #1 (omitted)
# rna_2 = enc + "RNAseq_HEK293_SRR1630838.bam"              # RNA-seq #2 (omitted)
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_1_new_peaks.narrowPeak"
GG3h_npk2 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_2_new_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = labhome + "Alu_ana_4_gh2ax_53bp1/"
os.makedirs(ana) if not os.path.exists(ana) else None


""" ############################################################################################ """
""" For target sites separated by more than 2MB, determine 53BP1 and gH2AX span width """
ana_1 = ana + "1_span_width/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
hic.get_span_width(gen, h2GGin, h2WTin, ana_1 + "GG_gh2ax_width")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG), 2000000)
hic.get_span_width(gen, bpGGin, bpWTin, ana_1 + "GG_53bp1_width")

gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
hic.get_span_width(gen, h2CTin, h2WTin, ana_1 + "CT_gh2ax_width")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg38, AluCT), 2000000)
hic.get_span_width(gen, bpCTin, bpWTin, ana_1 + "CT_53bp1_width")

gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
hic.get_span_width(gen, h2TAin, h2WTin, ana_1 + "TA_gh2ax_width")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg38, AluTA), 2000000)
hic.get_span_width(gen, bpTAin, bpWTin, ana_1 + "TA_53bp1_width")


""" ############################################################################################ """
""" Obtain peak profiles of gH2AX and 53BP1 data """
ana_2 = ana + "2_wiggle_windows/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
c.to_wiggle_windows(h2GGin, ana_2 + "GG_gh2ax", 500)
c.to_wiggle_windows(bpGGin, ana_2 + "GG_53bp1", 500)
c.to_wiggle_windows(h2CTin, ana_2 + "CT_gh2ax", 500)
c.to_wiggle_windows(bpCTin, ana_2 + "CT_53bp1", 500)
c.to_wiggle_windows(h2TAin, ana_2 + "TA_gh2ax", 500)
c.to_wiggle_windows(bpTAin, ana_2 + "TA_53bp1", 500)
c.to_wiggle_windows(h2WTin, ana_2 + "WT_gh2ax", 500)
c.to_wiggle_windows(bpWTin, ana_2 + "WT_53bp1", 500)


""" ############################################################################################ """
""" Get all training and testing data from AluGG, AluTA, and AluCT """
ana_3 = ana + "3_counts/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), h2GGin, ana_3 + "GG-C9_gh2ax_rc.csv")
m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin, ana_3 + "GG-C9_53bp1_rc.csv")
m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), h2TAin, ana_3 + "TA-C9_gh2ax_rc.csv")
m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), bpTAin, ana_3 + "TA-C9_53bp1_rc.csv")
m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), h2CTin, ana_3 + "CT-C9_gh2ax_rc.csv")
m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), bpCTin, ana_3 + "CT-C9_53bp1_rc.csv")















""" """
# path_h5 = desktop + "SRRrun123_merged_corrected_ICE_-2_5.h5"
# m.h5_fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_h5, "chr7", 5529660, 5000000)
# m.h5_fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_h5, "chr8", 127736258, 5000000)

""" Obtain 4C-seq profiles from Hi-C data """
# path_hic = "/Users/rogerzou/Downloads/K562"
# path_out = "/Users/rogerzou/Downloads/4Cseq_GG"
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr7", 5529660, 5000000)     # ACTB cleavage site
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr8", 127736258, 5000000)   # MYC cleavage site
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr1", 89600000, 5000000)
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr1", 90400000, 5000000)
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr17", 57500000, 5000000)
# m.fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_hic, 5, "chr21", 25500000, 5000000)

""" """
# dist = 3000000
# kb_res = 10
# gg_mre11_filt = m.gen_filter_dist(m.macs_gen(GG_mre11, 10000, hg38, AluGG), dist)
# m.fourCseq_gen(gg_mre11_filt, path_out, path_hic, kb_res, dist)
