"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.mtss as m


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana, enc, enc_a = labhome + "Alu_analysis/", labhome + "public/", labhome + "public/analysis/"
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

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

path_h5 = desktop + "SRRrun123_merged_corrected_ICE_-2_5.h5"
m.h5_fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_h5, "chr7", 5529660, 5000000)
m.h5_fourCseq_single("/Users/rogerzou/Downloads/4Cseq", path_h5, "chr8", 127736258, 5000000)

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


""" Obtain peak profiles of gH2AX and 53BP1 data """
# c.to_wiggle_windows(h2GGin, desktop + "GG_gh2ax", 5000)
# c.to_wiggle_windows(bpGGin, desktop + "GG_53bp1", 5000)



""" get all training and testing data from AluGG, AluTA, and AluCT """
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), h2GGin, ana + "GG-C9_gh2ax_rc.csv")
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin, ana + "GG-C9_53bp1_rc.csv")
# m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), h2TAin, ana + "TA-C9_gh2ax_rc.csv")
# m.read_counts(m.macs_gen(TA_cas9, 10000, hg38, AluTA), bpTAin, ana + "TA-C9_53bp1_rc.csv")
# m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), h2CTin, ana + "CT-C9_gh2ax_rc.csv")
# m.read_counts(m.macs_gen(CT_cas9, 10000, hg38, AluCT), bpCTin, ana + "CT-C9_53bp1_rc.csv")
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin_noD, ana + "GG-C9-noD_53bp1_rc.csv")
# m.read_counts(m.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin_PKi, ana + "GG-C9-PKi_53bp1_rc.csv")
