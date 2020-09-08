"""
Script for:
(1) Determination of genome-wide target sites (both on- and off-target) for each gRNA.
(2) Determine enrichment of MRE11, Cas9, and ENCODE epigenetic data at each target site.
(3) Develop machine learning models to predict MRE11/Cas9 enrichment solely using ENCODE.

"""

import src.mtss as m


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana, enc, enc_a = labhome + "Alu_analysis2/", labhome + "public/", labhome + "public/analysis/"
mreGGbam = labhome + "200206_chipseq/16_AluGG-MRE11_final.bam"
casGGbam = labhome + "200206_chipseq/15_AluGG-Cas9_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_rmdup.bam"
casTAbam = labhome + "200316_chipseq/AluTA-cas9-rep1_rmdup.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_rmdup.bam"
casCTbam = labhome + "200316_chipseq/AluCT-cas9-rep1_rmdup.bam"
mreGGbam_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_rmdup.bam"
mreGGbam_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_rmdup.bam"
casGGnpk = labhome + "200206_chipseq/macs/15_AluGG-Cas9_peaks.narrowPeak"
casCTnpk = labhome + "200316_chipseq/macs/AluCT-Cas9-rep1_peaks.narrowPeak"
casTAnpk = labhome + "200316_chipseq/macs/AluTA-Cas9-rep1_peaks.narrowPeak"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"


""" get all training and testing data from AluGG, AluTA, and AluCT """
m.read_subsets(m.macs_gen(casGGnpk, 750, hg38, AluGG), casGGbam, ana + "GG-C9_cas9_rs")
m.read_subsets(m.macs_gen(casTAnpk, 750, hg38, AluTA), casTAbam, ana + "TA-C9_cas9_rs")
m.read_subsets(m.macs_gen(casCTnpk, 750, hg38, AluCT), casCTbam, ana + "CT-C9_cas9_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam, ana + "GG-C9_mre11_rs")
m.read_subsets(m.macs_gen(casTAnpk, 1250, hg38, AluTA), mreTAbam, ana + "TA-C9_mre11_rs")
m.read_subsets(m.macs_gen(casCTnpk, 1250, hg38, AluCT), mreCTbam, ana + "CT-C9_mre11_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam_nD, ana + "GG-C9-noD_mre11_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam_PK, ana + "GG-C9-PKi_mre11_rs")

m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), h3k4me1_1, ana + "GG-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), h3k4me3_1, ana + "GG-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), h3k9me3_1, ana + "GG-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), h3k27ac_1, ana + "GG-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), h3k36me3_1, ana + "GG-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 200, hg38, AluGG), dnasei_1, ana + "GG-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 200, hg38, AluGG), mnase_1, ana + "GG-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 5000, hg38, AluGG), atac_1, ana + "GG-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 2500, hg38, AluGG), rna_3, ana + "GG-C9_rna_rc.csv")

m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), h3k4me1_1, ana + "TA-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), h3k4me3_1, ana + "TA-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), h3k9me3_1, ana + "TA-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), h3k27ac_1, ana + "TA-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), h3k36me3_1, ana + "TA-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 200, hg38, AluTA), dnasei_1, ana + "TA-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 200, hg38, AluTA), mnase_1, ana + "TA-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 5000, hg38, AluTA), atac_1, ana + "TA-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 2500, hg38, AluTA), rna_3, ana + "TA-C9_rna_rc.csv")

m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), h3k4me1_1, ana + "CT-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), h3k4me3_1, ana + "CT-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), h3k9me3_1, ana + "CT-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), h3k27ac_1, ana + "CT-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), h3k36me3_1, ana + "CT-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 200, hg38, AluCT), dnasei_1, ana + "CT-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 200, hg38, AluCT), mnase_1, ana + "CT-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 5000, hg38, AluCT), atac_1, ana + "CT-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 2500, hg38, AluCT), rna_3, ana + "CT-C9_rna_rc.csv")

""" One-mismatch annotation """
m.read_mismatch(m.macs_gen(casGGnpk, 750, hg38, AluGG), ana + "GG-C9_mismatch.csv")
m.read_mismatch(m.macs_gen(casTAnpk, 750, hg38, AluTA), ana + "TA-C9_mismatch.csv")
m.read_mismatch(m.macs_gen(casCTnpk, 750, hg38, AluCT), ana + "CT-C9_mismatch.csv")

""" ChromHMM epigenetic chromatin state annotation """
m.read_chromhmm(m.macs_gen(casGGnpk, 750, hg38, AluGG), ana + "GG-C9_chromhmm.csv")
m.read_chromhmm(m.macs_gen(casTAnpk, 750, hg38, AluTA), ana + "TA-C9_chromhmm.csv")
m.read_chromhmm(m.macs_gen(casCTnpk, 750, hg38, AluCT), ana + "CT-C9_chromhmm.csv")

""" Set arrays for epigenetic data """
B1 = m.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
B4 = m.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
B5 = m.load_nparray(ana + "GG-C9_h3k9me3_rc.csv")
B6 = m.load_nparray(ana + "GG-C9_h3k27ac_rc.csv")
B7 = m.load_nparray(ana + "GG-C9_h3k36me3_rc.csv")
B8 = m.load_nparray(ana + "GG-C9_dnasei_rc.csv")
B9 = m.load_nparray(ana + "GG-C9_mnase_rc.csv")
B10 = m.load_nparray(ana + "GG-C9_atac_rc.csv")
B11 = m.load_nparray(ana + "GG-C9_rna_rc.csv")
B12 = m.load_nparray(ana + "GG-C9_chromhmm.csv")
B13 = m.load_nparray(ana + "GG-C9_mismatch.csv")
B_GG = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13]

B1 = m.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
B4 = m.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
B5 = m.load_nparray(ana + "GG-C9_h3k9me3_rc.csv")
B6 = m.load_nparray(ana + "GG-C9_h3k27ac_rc.csv")
B7 = m.load_nparray(ana + "GG-C9_h3k36me3_rc.csv")
B8 = m.load_nparray(ana + "GG-C9_dnasei_rc.csv")
B9 = m.load_nparray(ana + "GG-C9_mnase_rc.csv")
B10 = m.load_nparray(ana + "GG-C9_atac_rc.csv")
B11 = m.load_nparray(ana + "GG-C9_rna_rc.csv")
B12 = m.load_nparray(ana + "GG-C9_chromhmm.csv")
B13 = m.load_nparray(ana + "GG-C9_mismatch.csv")
B_TA = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13]

B1 = m.load_nparray(ana + "CT-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana + "CT-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana + "CT-C9_h3k9me3_rc.csv")
B4 = m.load_nparray(ana + "CT-C9_h3k27ac_rc.csv")
B5 = m.load_nparray(ana + "CT-C9_h3k36me3_rc.csv")
B6 = m.load_nparray(ana + "CT-C9_dnasei_rc.csv")
B7 = m.load_nparray(ana + "CT-C9_mnase_rc.csv")
B8 = m.load_nparray(ana + "CT-C9_atac_rc.csv")
B9 = m.load_nparray(ana + "CT-C9_rna_rc.csv")
B10 = m.load_nparray(ana + "CT-C9_chromhmm.csv")
B11 = m.load_nparray(ana + "CT-C9_mismatch.csv")
B_CT = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13]

num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"

""" generate merged datasets for MRE11 """
head = m.load_npheader(ana + "GG-C9_mre11_rs.csv") + rc_head
GG_mre = m.load_nparray(ana + "GG-C9_mre11_rs.csv")
GG_mre_m = m.mergesubsetcounts(GG_mre, B_GG, num_cols, ana + "GG-mre-merged.csv", head)
TA_mre = m.load_nparray(ana + "TA-C9_mre11_rs.csv")
TA_mre_m = m.mergesubsetcounts(TA_mre, B_TA, num_cols, ana + "TA-mre-merged.csv", head)
CT_mre = m.load_nparray(ana + "CT-C9_mre11_rs.csv")
CT_mre_m = m.mergesubsetcounts(CT_mre, B_CT, num_cols, ana + "CT-mre-merged.csv", head)
m.mergerows([GG_mre_m, TA_mre_m, CT_mre_m], ana + "ALL-mre-merged.csv", head)

""" generate merged datasets for Cas9 """
head = m.load_npheader(ana + "GG-C9_cas9_rs.csv") + rc_head
GG_cas = m.load_nparray(ana + "GG-C9_cas9_rs.csv")
GG_cas_m = m.mergesubsetcounts(GG_cas, B_GG, num_cols, ana + "GG-cas-merged.csv", head)
TA_cas = m.load_nparray(ana + "TA-C9_cas9_rs.csv")
TA_cas_m = m.mergesubsetcounts(TA_cas, B_TA, num_cols, ana + "TA-cas-merged.csv", head)
CT_cas = m.load_nparray(ana + "CT-C9_cas9_rs.csv")
CT_cas_m = m.mergesubsetcounts(CT_cas, B_CT, num_cols, ana + "CT-cas-merged.csv", head)
m.mergerows([GG_cas_m, TA_cas_m, CT_cas_m], ana + "ALL-cas-merged.csv", head)

