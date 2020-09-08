"""
Script for:
(1) Determine Cas9 binding and genome-editing kinetics genome-wide

"""

import src.mtss as m
import numpy as np

""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana, enc, enc_a = labhome + "Alu_kinetics/", labhome + "public/", labhome + "public/analysis/"
casneg = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/200212_chipseq_WT1/A16_WT1C_final.bam"
mreneg = "/Volumes/Lab-Home/rzou4/NGS_data/2_pcl/200212_chipseq_WT1/A17_WT1M_final.bam"
casGG30m_1 = labhome + "200804_chipseq/A01_rmdup.bam"
casGG1h_1 = labhome + "200804_chipseq/A02_rmdup.bam"
casGG3h_1 = labhome + "200804_chipseq/A03_rmdup.bam"
mreGG30m_1 = labhome + "200804_chipseq/A04_rmdup.bam"
mreGG1h_1 = labhome + "200804_chipseq/A05_rmdup.bam"
mreGG3h_1 = labhome + "200804_chipseq/A06_rmdup.bam"
casGG30m_2 = labhome + "200804_chipseq/A07_rmdup.bam"
casGG1h_2 = labhome + "200804_chipseq/A08_rmdup.bam"
casGG3h_2 = labhome + "200804_chipseq/A09_rmdup.bam"
mreGG30m_2 = labhome + "200804_chipseq/A10_rmdup.bam"
mreGG1h_2 = labhome + "200804_chipseq/A11_rmdup.bam"
mreGG3h_2 = labhome + "200804_chipseq/A12_rmdup.bam"
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

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_1_new_peaks.narrowPeak"
GG3h_npk2 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_2_new_peaks.narrowPeak"

""" Get read subsets for time-resolved dCas9/MRE11 ChIP-seq for AluGG """
m.read_subsets(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), casneg, ana + "GG-C9_cas9_rs_00m_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), casGG30m_1, ana + "GG-C9_cas9_rs_30m_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), casGG1h_1, ana + "GG-C9_cas9_rs_1h_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), casGG3h_1, ana + "GG-C9_cas9_rs_3h_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreneg, ana + "GG-C9_mre11_rs_00m_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG30m_1, ana + "GG-C9_mre11_rs_30m_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG1h_1, ana + "GG-C9_mre11_rs_1h_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG3h_1, ana + "GG-C9_mre11_rs_3h_1")

m.read_subsets(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), casneg, ana + "GG-C9_cas9_rs_00m_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), casGG30m_2, ana + "GG-C9_cas9_rs_30m_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), casGG1h_2, ana + "GG-C9_cas9_rs_1h_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), casGG3h_2, ana + "GG-C9_cas9_rs_3h_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 1250, hg38, AluGG, fenr=8), mreneg, ana + "GG-C9_mre11_rs_00m_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 1250, hg38, AluGG, fenr=8), mreGG30m_2, ana + "GG-C9_mre11_rs_30m_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 1250, hg38, AluGG, fenr=8), mreGG1h_2, ana + "GG-C9_mre11_rs_1h_2")
m.read_subsets(m.macs_gen(GG3h_npk2, 1250, hg38, AluGG, fenr=8), mreGG3h_2, ana + "GG-C9_mre11_rs_3h_2")

""" Quantify kinetics for each set """
subset_cas9_1 = [ana + "GG-C9_cas9_rs_00m_1.csv",
                 ana + "GG-C9_cas9_rs_30m_1.csv",
                 ana + "GG-C9_cas9_rs_1h_1.csv",
                 ana + "GG-C9_cas9_rs_3h_1.csv"]
m.read_kinetics(subset_cas9_1, ana + "GG-C9_cas9_rc_kin_1")
subset_mre11_1 = [ana + "GG-C9_mre11_rs_00m_1.csv",
                  ana + "GG-C9_mre11_rs_30m_1.csv",
                  ana + "GG-C9_mre11_rs_1h_1.csv",
                  ana + "GG-C9_mre11_rs_3h_1.csv"]
m.read_kinetics(subset_mre11_1, ana + "GG-C9_mre11_rc_kin_1")

subset_cas9_2 = [ana + "GG-C9_cas9_rs_00m_2.csv",
                 ana + "GG-C9_cas9_rs_30m_2.csv",
                 ana + "GG-C9_cas9_rs_1h_2.csv",
                 ana + "GG-C9_cas9_rs_3h_2.csv"]
m.read_kinetics(subset_cas9_2, ana + "GG-C9_cas9_rc_kin_2")
subset_mre11_2 = [ana + "GG-C9_mre11_rs_00m_2.csv",
                  ana + "GG-C9_mre11_rs_30m_2.csv",
                  ana + "GG-C9_mre11_rs_1h_2.csv",
                  ana + "GG-C9_mre11_rs_3h_2.csv"]
m.read_kinetics(subset_mre11_2, ana + "GG-C9_mre11_rc_kin_2")


""" One to Two mismatches annotations """
m.read_mismatch(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), ana + "GG-C9_mismatch_1.csv")
m.read_mismatch(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), ana + "GG-C9_mismatch_2.csv")

""" ChromHMM epigenetic chromatin state annotation """
m.read_chromhmm(m.macs_gen(GG3h_npk1, 750, hg38, AluGG, fenr=8), ana + "GG-C9_chromhmm_1.csv")
m.read_chromhmm(m.macs_gen(GG3h_npk2, 750, hg38, AluGG, fenr=8), ana + "GG-C9_chromhmm_2.csv")

""" Get epigenetic info for replicate 1 - AluGG """
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), h3k4me1_1, ana + "GG-C9_h3k4me1_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), h3k4me3_1, ana + "GG-C9_h3k4me3_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), h3k9me3_1, ana + "GG-C9_h3k9me3_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), h3k27ac_1, ana + "GG-C9_h3k27ac_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), h3k36me3_1, ana + "GG-C9_h3k36me3_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50, hg38, AluGG, fenr=8), dnasei_1, ana + "GG-C9_dnasei_50_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50, hg38, AluGG, fenr=8), mnase_1, ana + "GG-C9_mnase_50_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), atac_1, ana + "GG-C9_atac_50000_1.csv")
m.read_counts(m.macs_gen(GG3h_npk1, 50000, hg38, AluGG, fenr=8), rna_3, ana + "GG-C9_rna_50000_1.csv")

""" Get epigenetic info for replicate 2 - AluGG """
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), h3k4me1_1, ana + "GG-C9_h3k4me1_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), h3k4me3_1, ana + "GG-C9_h3k4me3_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), h3k9me3_1, ana + "GG-C9_h3k9me3_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), h3k27ac_1, ana + "GG-C9_h3k27ac_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), h3k36me3_1, ana + "GG-C9_h3k36me3_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50, hg38, AluGG, fenr=8), dnasei_1, ana + "GG-C9_dnasei_50_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50, hg38, AluGG, fenr=8), mnase_1, ana + "GG-C9_mnase_50_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), atac_1, ana + "GG-C9_atac_50000_2.csv")
m.read_counts(m.macs_gen(GG3h_npk2, 50000, hg38, AluGG, fenr=8), rna_3, ana + "GG-C9_rna_50000_2.csv")

""" Set arrays for final epigenetic data """
B1 = m.load_nparray(ana + "GG-C9_h3k4me1_50000_1.csv")
B2 = m.load_nparray(ana + "GG-C9_h3k4me3_50000_1.csv")
B3 = m.load_nparray(ana + "GG-C9_h3k9me3_50000_1.csv")
B4 = m.load_nparray(ana + "GG-C9_h3k27ac_50000_1.csv")
B5 = m.load_nparray(ana + "GG-C9_h3k36me3_50000_1.csv")
B6 = m.load_nparray(ana + "GG-C9_dnasei_50_1.csv")
B7 = m.load_nparray(ana + "GG-C9_mnase_50_1.csv")
B8 = m.load_nparray(ana + "GG-C9_atac_50000_1.csv")
B9 = m.load_nparray(ana + "GG-C9_rna_50000_1.csv")
B10 = m.load_nparray(ana + "GG-C9_chromhmm_1.csv")
B11 = m.load_nparray(ana + "GG-C9_mismatch_1.csv")
B_GG_1 = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]

B1 = m.load_nparray(ana + "GG-C9_h3k4me1_50000_2.csv")
B2 = m.load_nparray(ana + "GG-C9_h3k4me3_50000_2.csv")
B3 = m.load_nparray(ana + "GG-C9_h3k9me3_50000_2.csv")
B4 = m.load_nparray(ana + "GG-C9_h3k27ac_50000_2.csv")
B5 = m.load_nparray(ana + "GG-C9_h3k36me3_50000_2.csv")
B6 = m.load_nparray(ana + "GG-C9_dnasei_50_2.csv")
B7 = m.load_nparray(ana + "GG-C9_mnase_50_2.csv")
B8 = m.load_nparray(ana + "GG-C9_atac_50000_2.csv")
B9 = m.load_nparray(ana + "GG-C9_rna_50000_2.csv")
B10 = m.load_nparray(ana + "GG-C9_chromhmm_2.csv")
B11 = m.load_nparray(ana + "GG-C9_mismatch_2.csv")
B_GG_2 = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]

num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"

""" generate merged datasets for MRE11 and Cas9 """
head = m.load_npheader(ana + "GG-C9_mre11_rc_kin_1.csv") + rc_head
GG_mre = m.load_nparray(ana + "GG-C9_mre11_rc_kin_1.csv")
m.mergesubsetcounts(GG_mre, B_GG_1, num_cols, ana + "GG-mre-merged_kinetics_1.csv", head)
head = m.load_npheader(ana + "GG-C9_cas9_rc_kin_1.csv") + rc_head
GG_cas = m.load_nparray(ana + "GG-C9_cas9_rc_kin_1.csv")
m.mergesubsetcounts(GG_cas, B_GG_1, num_cols, ana + "GG-cas-merged_kinetics_1.csv", head)

head = m.load_npheader(ana + "GG-C9_mre11_rc_kin_2.csv") + rc_head
GG_mre = m.load_nparray(ana + "GG-C9_mre11_rc_kin_2.csv")
m.mergesubsetcounts(GG_mre, B_GG_2, num_cols, ana + "GG-mre-merged_kinetics_2.csv", head)
head = m.load_npheader(ana + "GG-C9_cas9_rc_kin_2.csv") + rc_head
GG_cas = m.load_nparray(ana + "GG-C9_cas9_rc_kin_2.csv")
m.mergesubsetcounts(GG_cas, B_GG_2, num_cols, ana + "GG-cas-merged_kinetics_2.csv", head)


""" ################# RANDOM FOREST MRE11 2 MISMATCHES OR FEWER ################# """
""" Random Forest - MRE11 - 2 mismatches or fewer - YES epigenetics - YES exact mismatches"""
X1, y1, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_1.csv", epi=True, mm=2)
X2, y2, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_2.csv", epi=True, mm=2)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_mre_2orless_TT.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_mre_2orless_TT.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_mre_2orless_TT.sav", alllabels, count=20)

""" Random Forest - MRE11 - 2 mismatches or fewer - YES epigenetics - NO exact mismatches"""
X1, y1, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_1.csv", epi=True, mm=1)
X2, y2, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_2.csv", epi=True, mm=1)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_mre_2orless_TF.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_mre_2orless_TF.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_mre_2orless_TF.sav", alllabels, count=20)

""" Random Forest - MRE11 - 2 mismatches or fewer - YES epigenetics - NO mismatches"""
X1, y1, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_1.csv", epi=True, mm=0)
X2, y2, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_2.csv", epi=True, mm=0)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_mre_2orless_T-.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_mre_2orless_T-.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_mre_2orless_T-.sav", alllabels, count=20)

""" Random Forest - MRE11 - 2 mismatches or fewer - NO epigenetics - YES exact mismatches"""
X1, y1, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_1.csv", epi=False, mm=2)
X2, y2, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_2.csv", epi=False, mm=2)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_mre_2orless_FT.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_mre_2orless_FT.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_mre_2orless_FT.sav", alllabels, count=20)

""" Random Forest - MRE11 - 2 mismatches or fewer - NO epigenetics - NO exact mismatches"""
X1, y1, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_1.csv", epi=False, mm=1)
X2, y2, alllabels = m.getXy_2orLess(ana + "GG-mre-merged_kinetics_2.csv", epi=False, mm=1)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_mre_2orless_FF.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_mre_2orless_FF.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_mre_2orless_FF.sav", alllabels, count=20)

""" ################# RANDOM FOREST CAS9 ALL DATA ################# """
""" Random Forest - Cas9 - all data - YES epigenetics - YES exact mismatches"""
X1, y1, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_1.csv", epi=True, mm=2)
X2, y2, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_2.csv", epi=True, mm=2)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_cas_all_TT.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_cas_all_TT.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_cas_all_TT.sav", alllabels, count=20)

""" Random Forest - Cas9 - all data - YES epigenetics - NO exact mismatches"""
X1, y1, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_1.csv", epi=True, mm=1)
X2, y2, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_2.csv", epi=True, mm=1)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_cas_all_TF.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_cas_all_TF.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_cas_all_TF.sav", alllabels, count=20)

""" Random Forest - Cas9 - all data - YES epigenetics - NO mismatches"""
X1, y1, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_1.csv", epi=True, mm=0)
X2, y2, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_2.csv", epi=True, mm=0)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_cas_all_T-.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_cas_all_T-.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_cas_all_T-.sav", alllabels, count=20)

""" Random Forest - Cas9 - all data - NO epigenetics - YES exact mismatches"""
X1, y1, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_1.csv", epi=False, mm=2)
X2, y2, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_2.csv", epi=False, mm=2)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_cas_all_FT.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_cas_all_FT.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_cas_all_FT.sav", alllabels, count=20)

""" Random Forest - Cas9 - all data - NO epigenetics - NO exact mismatches"""
X1, y1, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_1.csv", epi=False, mm=1)
X2, y2, alllabels = m.getXy_all(ana + "GG-cas-merged_kinetics_2.csv", epi=False, mm=1)
X = np.concatenate((X1, X2))
y = np.concatenate((y1, y2))
X_train, X_test, y_train, y_test = m.data_split(X, y)
# m.RandomForestTrainGridCV(X_train, y_train, ana + "gridRF_GG_cas_all_FF.sav")
m.ModelTest(X_test, y_test, ana + "gridRF_GG_cas_all_FF.sav")
m.FeatureImportance(X_test, y_test, ana + "gridRF_GG_cas_all_FF.sav", alllabels, count=20)


""" Split each time point into separate csv files, with each column corresponding to a different 
    value of under column 'mmtype' """
collist = ['timepoint_1', 'timepoint_2', 'timepoint_3', 'timepoint_4']
m.sort_mmtype_column(ana + "GG-cas-merged_kinetics_1.csv", collist)
m.sort_mmtype_column(ana + "GG-cas-merged_kinetics_2.csv", collist)
m.sort_mmtype_column(ana + "GG-mre-merged_kinetics_1.csv", collist)
m.sort_mmtype_column(ana + "GG-mre-merged_kinetics_2.csv", collist)
