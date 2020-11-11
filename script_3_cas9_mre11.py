"""
Script for:
(1) Determination of genome-wide target sites (both on- and off-target) for each gRNA.
(2) Determination of MRE11 and Cas9 enrichment, and ENCODE epigenetic data at each target site.
"""

import src.mtss as m
import src.msa as msa
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
enc, enc_a = labhome + "public/", labhome + "public/analysis/"
mreGGbam = labhome + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
casGGbam = labhome + "200206_chipseq/AluGG-Cas9_hg38_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_hg38_final.bam"
casTAbam = labhome + "200316_chipseq/AluTA-cas9-rep1_hg38_final.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_hg38_final.bam"
casCTbam = labhome + "200316_chipseq/AluCT-cas9-rep1_hg38_final.bam"
mreGGbam_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_hg38_final.bam"
mreGGbam_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_hg38_final.bam"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR6418075.bam"              # ATAC-seq (medium deep)
atac_2 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq (bad)
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

""" macs2 peak detection """
# for on-target analysis at GG, CT, TA
casGGnpk = labhome + "200206_chipseq/macs/AluGG-Cas9_hg38_final_peaks.narrowPeak"
casCTnpk = labhome + "200316_chipseq/macs/AluCT-cas9-rep1_hg38_final_peaks.narrowPeak"
casTAnpk = labhome + "200316_chipseq/macs/AluTA-cas9-rep1_hg38_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = labhome + "Alu_ana_3_cas9_mre11/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_msa/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_subsets/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_epigenetics/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_merged/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


""" ############################################################################################ """
""" Determine multiple sequence alignments (up to 300) for each ChIP-seq paired-end reads around
    both on- and off-target sites (Figure 1X) """
msa.get_bamfile_pe_reads(m.macs_gen(casGGnpk, 750, hg38, AluGG), casGGbam, ana_1 + "GG-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casCTnpk, 750, hg38, AluCT), casCTbam, ana_1 + "CT-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casTAnpk, 750, hg38, AluTA), casTAbam, ana_1 + "TA-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam, ana_1 + "GG-C9_mre11")
msa.get_bamfile_pe_reads(m.macs_gen(casCTnpk, 1250, hg38, AluCT), mreCTbam, ana_1 + "CT-C9_mre11")
msa.get_bamfile_pe_reads(m.macs_gen(casTAnpk, 1250, hg38, AluTA), mreTAbam, ana_1 + "TA-C9_mre11")
msa.bowtie2_msa_paired(ana_1 + "GG-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "GG-C9_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-C9_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-C9_mre11", hg38[1])
msa.parse_msa_sam_paired(ana_1 + "GG-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "GG-C9_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-C9_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "GG-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "CT-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "TA-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "GG-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "CT-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "TA-C9_mre11_msa")


""" ############################################################################################ """
""" get all training and testing data from AluGG, AluTA, and AluCT """
m.read_subsets(m.macs_gen(casGGnpk, 750, hg38, AluGG), casGGbam, ana_2 + "GG-C9_cas9_rs")
m.read_subsets(m.macs_gen(casCTnpk, 750, hg38, AluCT), casCTbam, ana_2 + "CT-C9_cas9_rs")
m.read_subsets(m.macs_gen(casTAnpk, 750, hg38, AluTA), casTAbam, ana_2 + "TA-C9_cas9_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam, ana_2 + "GG-C9_mre11_rs")
m.read_subsets(m.macs_gen(casCTnpk, 1250, hg38, AluCT), mreCTbam, ana_2 + "CT-C9_mre11_rs")
m.read_subsets(m.macs_gen(casTAnpk, 1250, hg38, AluTA), mreTAbam, ana_2 + "TA-C9_mre11_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam_nD, ana_2 + "GG-C9-noD_mre11_rs")
m.read_subsets(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam_PK, ana_2 + "GG-C9-PKi_mre11_rs")


""" ############################################################################################ """
""" Determine epigenetic, mismatch, and chromatin information """
# AluGG datasets
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), h3k4me1_1, ana_3 + "GG-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), h3k4me3_1, ana_3 + "GG-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), h3k9me3_1, ana_3 + "GG-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), h3k27ac_1, ana_3 + "GG-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), h3k36me3_1, ana_3 + "GG-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50, hg38, AluGG), dnasei_1, ana_3 + "GG-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50, hg38, AluGG), mnase_1, ana_3 + "GG-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), atac_1, ana_3 + "GG-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casGGnpk, 50000, hg38, AluGG), rna_3, ana_3 + "GG-C9_rna_rc.csv")
# AluCT datasets
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), h3k4me1_1, ana_3 + "CT-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), h3k4me3_1, ana_3 + "CT-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), h3k9me3_1, ana_3 + "CT-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), h3k27ac_1, ana_3 + "CT-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), h3k36me3_1, ana_3 + "CT-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50, hg38, AluCT), dnasei_1, ana_3 + "CT-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50, hg38, AluCT), mnase_1, ana_3 + "CT-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), atac_1, ana_3 + "CT-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casCTnpk, 50000, hg38, AluCT), rna_3, ana_3 + "CT-C9_rna_rc.csv")
# AluTA datasets
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), h3k4me1_1, ana_3 + "TA-C9_h3k4me1_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), h3k4me3_1, ana_3 + "TA-C9_h3k4me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), h3k9me3_1, ana_3 + "TA-C9_h3k9me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), h3k27ac_1, ana_3 + "TA-C9_h3k27ac_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), h3k36me3_1, ana_3 + "TA-C9_h3k36me3_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50, hg38, AluTA), dnasei_1, ana_3 + "TA-C9_dnasei_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50, hg38, AluTA), mnase_1, ana_3 + "TA-C9_mnase_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), atac_1, ana_3 + "TA-C9_atac_rc.csv")
m.read_counts(m.macs_gen(casTAnpk, 50000, hg38, AluTA), rna_3, ana_3 + "TA-C9_rna_rc.csv")
# One-mismatch annotation
m.read_mismatch(m.macs_gen(casGGnpk, 750, hg38, AluGG), ana_3 + "GG-C9_mismatch.csv")
m.read_mismatch(m.macs_gen(casCTnpk, 750, hg38, AluCT), ana_3 + "CT-C9_mismatch.csv")
m.read_mismatch(m.macs_gen(casTAnpk, 750, hg38, AluTA), ana_3 + "TA-C9_mismatch.csv")
# ChromHMM epigenetic chromatin state annotation
m.read_chromhmm(m.macs_gen(casGGnpk, 750, hg38, AluGG), ana_3 + "GG-C9_chromhmm.csv")
m.read_chromhmm(m.macs_gen(casCTnpk, 750, hg38, AluCT), ana_3 + "CT-C9_chromhmm.csv")
m.read_chromhmm(m.macs_gen(casTAnpk, 750, hg38, AluTA), ana_3 + "TA-C9_chromhmm.csv")

""" ############################################################################################ """
""" Generate merged datasets """
# AluGG datasets
B1 = m.load_nparray(ana_3 + "GG-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana_3 + "GG-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana_3 + "GG-C9_h3k9me3_rc.csv")
B4 = m.load_nparray(ana_3 + "GG-C9_h3k27ac_rc.csv")
B5 = m.load_nparray(ana_3 + "GG-C9_h3k36me3_rc.csv")
B6 = m.load_nparray(ana_3 + "GG-C9_dnasei_rc.csv")
B7 = m.load_nparray(ana_3 + "GG-C9_mnase_rc.csv")
B8 = m.load_nparray(ana_3 + "GG-C9_atac_rc.csv")
B9 = m.load_nparray(ana_3 + "GG-C9_rna_rc.csv")
B10 = m.load_nparray(ana_3 + "GG-C9_chromhmm.csv")
B11 = m.load_nparray(ana_3 + "GG-C9_mismatch.csv")
B_GG = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
# AluCT datasets
B1 = m.load_nparray(ana_3 + "CT-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana_3 + "CT-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana_3 + "CT-C9_h3k9me3_rc.csv")
B4 = m.load_nparray(ana_3 + "CT-C9_h3k27ac_rc.csv")
B5 = m.load_nparray(ana_3 + "CT-C9_h3k36me3_rc.csv")
B6 = m.load_nparray(ana_3 + "CT-C9_dnasei_rc.csv")
B7 = m.load_nparray(ana_3 + "CT-C9_mnase_rc.csv")
B8 = m.load_nparray(ana_3 + "CT-C9_atac_rc.csv")
B9 = m.load_nparray(ana_3 + "CT-C9_rna_rc.csv")
B10 = m.load_nparray(ana_3 + "CT-C9_chromhmm.csv")
B11 = m.load_nparray(ana_3 + "CT-C9_mismatch.csv")
B_CT = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
# AluTA datasets
B1 = m.load_nparray(ana_3 + "TA-C9_h3k4me1_rc.csv")
B2 = m.load_nparray(ana_3 + "TA-C9_h3k4me3_rc.csv")
B3 = m.load_nparray(ana_3 + "TA-C9_h3k9me3_rc.csv")
B4 = m.load_nparray(ana_3 + "TA-C9_h3k27ac_rc.csv")
B5 = m.load_nparray(ana_3 + "TA-C9_h3k36me3_rc.csv")
B6 = m.load_nparray(ana_3 + "TA-C9_dnasei_rc.csv")
B7 = m.load_nparray(ana_3 + "TA-C9_mnase_rc.csv")
B8 = m.load_nparray(ana_3 + "TA-C9_atac_rc.csv")
B9 = m.load_nparray(ana_3 + "TA-C9_rna_rc.csv")
B10 = m.load_nparray(ana_3 + "TA-C9_chromhmm.csv")
B11 = m.load_nparray(ana_3 + "TA-C9_mismatch.csv")
B_TA = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]

num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"

""" generate merged datasets for MRE11 """
head = m.load_npheader(ana_2 + "GG-C9_mre11_rs.csv") + rc_head
GG_mre = m.load_nparray(ana_2 + "GG-C9_mre11_rs.csv")
GG_mre_m = m.mergesubsetcounts(GG_mre, B_GG, num_cols, ana_4 + "GG-mre-merged.csv", head)
CT_mre = m.load_nparray(ana_2 + "CT-C9_mre11_rs.csv")
CT_mre_m = m.mergesubsetcounts(CT_mre, B_CT, num_cols, ana_4 + "CT-mre-merged.csv", head)
TA_mre = m.load_nparray(ana_2 + "TA-C9_mre11_rs.csv")
TA_mre_m = m.mergesubsetcounts(TA_mre, B_TA, num_cols, ana_4 + "TA-mre-merged.csv", head)
m.mergerows([GG_mre_m, CT_mre_m, TA_mre_m], ana_4 + "ALL-mre-merged.csv", head)

""" generate merged datasets for Cas9 """
head = m.load_npheader(ana_2 + "GG-C9_cas9_rs.csv") + rc_head
GG_cas = m.load_nparray(ana_2 + "GG-C9_cas9_rs.csv")
GG_cas_m = m.mergesubsetcounts(GG_cas, B_GG, num_cols, ana_4 + "GG-cas-merged.csv", head)
CT_cas = m.load_nparray(ana_2 + "CT-C9_cas9_rs.csv")
CT_cas_m = m.mergesubsetcounts(CT_cas, B_CT, num_cols, ana_4 + "CT-cas-merged.csv", head)
TA_cas = m.load_nparray(ana_2 + "TA-C9_cas9_rs.csv")
TA_cas_m = m.mergesubsetcounts(TA_cas, B_TA, num_cols, ana_4 + "TA-cas-merged.csv", head)
m.mergerows([GG_cas_m, CT_cas_m, TA_cas_m], ana_4 + "ALL-cas-merged.csv", head)
