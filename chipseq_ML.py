# MH2AX ChIP-seq


import argparse
import src.alu as alu
import src.chipseq as c


""" File paths """
gg = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/200206_chipseq/"
ct = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/200316_chipseq/"
ta = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/200316_chipseq/"
ana = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/Alu_analysis/"

enc = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/public/"
enc_a = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/public/analysis/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

GG_cas9 = gg + "macs/15_AluGG-Cas9_peaks.narrowPeak"
# gg_mre11 = gg + "macs/16_AluGG-MRE11_peaks.narrowPeak"
CT_cas9 = ct + "macs/AluCT-Cas9-rep1_peaks.narrowPeak"
# ct_mre11 = ct + "macs/AluCT-MRE11_peaks.narrowPeak"
TA_cas9 = ta + "macs/AluTA-Cas9-rep1_peaks.narrowPeak"
# ta_mre11 = ta + "macs/AluTA-MRE11_peaks.narrowPeak"


h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"  # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"  # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"  # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"  # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"  # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"
atac_1 = enc + "ATACseq_HEK293_SRR5627157.bam"
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"
# rna_1 = enc + "RNAseq_HEK293_SRR1264355.bam"
# rna_2 = enc + "RNAseq_HEK293_SRR1630838.bam"
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"

mreGGin = gg + "16_AluGG-MRE11_final.bam"
casGGin = gg + "15_AluGG-Cas9_final.bam"
mreTAin = ta + "AluTA-mre11-rep1_rmdup.bam"
casTAin = ta + "AluTA-cas9-rep1_rmdup.bam"
mreCTin = ct + "AluCT-mre11-rep1_rmdup.bam"
casCTin = ct + "AluCT-cas9-rep1_rmdup.bam"
mreGGin_noD = ta + "AluTA-mre11-rep1_rmdup.bam"

""" get all training and testing data from AluGG, AluTA, and AluCT """
# c.refseq_initialize()
# alu.alu_read_subsets(alu.macs_generator(GG_cas9, 2500, hg38, AluGG), casGGin, ana + "GG-C9_cas9_rs.csv")
# alu.alu_read_subsets(alu.macs_generator(TA_cas9, 2500, hg38, AluTA), casTAin, ana + "TA-C9_cas9_rs.csv")
# alu.alu_read_subsets(alu.macs_generator(CT_cas9, 2500, hg38, AluCT), casCTin, ana + "CT-C9_cas9_rs.csv")
# alu.alu_read_subsets(alu.macs_generator(GG_cas9, 2500, hg38, AluGG), mreGGin, ana + "GG-C9_mre11_rs.csv")
# alu.alu_read_subsets(alu.macs_generator(TA_cas9, 2500, hg38, AluTA), mreTAin, ana + "TA-C9_mre11_rs.csv")
# alu.alu_read_subsets(alu.macs_generator(CT_cas9, 2500, hg38, AluCT), mreCTin, ana + "CT-C9_mre11_rs.csv")

# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), h3k4me1_1, ana + "GG-C9_h3k4me1_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), h3k4me3_1, ana + "GG-C9_h3k4me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), h3k9me3_1, ana + "GG-C9_h3k9me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), h3k27ac_1, ana + "GG-C9_h3k27ac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), h3k36me3_1, ana + "GG-C9_h3k36me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 200, hg38, AluGG), dnasei_1, ana + "GG-C9_dnasei_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 200, hg38, AluGG), mnase_1, ana + "GG-C9_mnase_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 5000, hg38, AluGG), atac_1, ana + "GG-C9_atac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(GG_cas9, 2500, hg38, AluGG), rna_3, ana + "GG-C9_rna_rc.csv")

# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), h3k4me1_1, ana + "TA-C9_h3k4me1_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), h3k4me3_1, ana + "TA-C9_h3k4me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), h3k9me3_1, ana + "TA-C9_h3k9me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), h3k27ac_1, ana + "TA-C9_h3k27ac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), h3k36me3_1, ana + "TA-C9_h3k36me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 200, hg38, AluTA), dnasei_1, ana + "TA-C9_dnasei_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 200, hg38, AluTA), mnase_1, ana + "TA-C9_mnase_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 5000, hg38, AluTA), atac_1, ana + "TA-C9_atac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(TA_cas9, 2500, hg38, AluTA), rna_3, ana + "TA-C9_rna_rc.csv")

# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), h3k4me1_1, ana + "CT-C9_h3k4me1_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), h3k4me3_1, ana + "CT-C9_h3k4me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), h3k9me3_1, ana + "CT-C9_h3k9me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), h3k27ac_1, ana + "CT-C9_h3k27ac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), h3k36me3_1, ana + "CT-C9_h3k36me3_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 200, hg38, AluCT), dnasei_1, ana + "CT-C9_dnasei_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 200, hg38, AluCT), mnase_1, ana + "CT-C9_mnase_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 5000, hg38, AluCT), atac_1, ana + "CT-C9_atac_rc.csv")
# alu.alu_read_counts(alu.macs_generator(CT_cas9, 2500, hg38, AluCT), rna_3, ana + "CT-C9_rna_rc.csv")

alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), h3k4me1_1, ana + "chak-C9_h3k4me1_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), h3k4me3_1, ana + "chak-C9_h3k4me3_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), h3k9me3_1, ana + "chak-C9_h3k9me3_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), h3k27ac_1, ana + "chak-C9_h3k27ac_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), h3k36me3_1, ana + "chak-C9_h3k36me3_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(200, hg38), dnasei_1, ana + "chak-C9_dnasei_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(200, hg38), mnase_1, ana + "chak-C9_mnase_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(5000, hg38), atac_1, ana + "chak-C9_atac_rc.csv")
alu.alu_read_counts(alu.chakrabarti_generator(2500, hg38), rna_3, ana + "chak-C9_rna_rc.csv")


""" Set arrays for epigenetic data """
# B1 = alu.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
# B2 = alu.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
# B3 = alu.load_nparray(ana + "GG-C9_h3k9me3_rc.csv")
# B4 = alu.load_nparray(ana + "GG-C9_h3k27ac_rc.csv")
# B5 = alu.load_nparray(ana + "GG-C9_h3k36me3_rc.csv")
# B6 = alu.load_nparray(ana + "GG-C9_dnasei_rc.csv")
# B7 = alu.load_nparray(ana + "GG-C9_mnase_rc.csv")
# B8 = alu.load_nparray(ana + "GG-C9_atac_rc.csv")
# B9 = alu.load_nparray(ana + "GG-C9_rna_rc.csv")
# B_GG = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
#
# B1 = alu.load_nparray(ana + "TA-C9_h3k4me1_rc.csv")
# B2 = alu.load_nparray(ana + "TA-C9_h3k4me3_rc.csv")
# B3 = alu.load_nparray(ana + "TA-C9_h3k9me3_rc.csv")
# B4 = alu.load_nparray(ana + "TA-C9_h3k27ac_rc.csv")
# B5 = alu.load_nparray(ana + "TA-C9_h3k36me3_rc.csv")
# B6 = alu.load_nparray(ana + "TA-C9_dnasei_rc.csv")
# B7 = alu.load_nparray(ana + "TA-C9_mnase_rc.csv")
# B8 = alu.load_nparray(ana + "TA-C9_atac_rc.csv")
# B9 = alu.load_nparray(ana + "TA-C9_rna_rc.csv")
# B_TA = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
#
# B1 = alu.load_nparray(ana + "CT-C9_h3k4me1_rc.csv")
# B2 = alu.load_nparray(ana + "CT-C9_h3k4me3_rc.csv")
# B3 = alu.load_nparray(ana + "CT-C9_h3k9me3_rc.csv")
# B4 = alu.load_nparray(ana + "CT-C9_h3k27ac_rc.csv")
# B5 = alu.load_nparray(ana + "CT-C9_h3k36me3_rc.csv")
# B6 = alu.load_nparray(ana + "CT-C9_dnasei_rc.csv")
# B7 = alu.load_nparray(ana + "CT-C9_mnase_rc.csv")
# B8 = alu.load_nparray(ana + "CT-C9_atac_rc.csv")
# B9 = alu.load_nparray(ana + "CT-C9_rna_rc.csv")
# B_CT = [B1, B2, B3, B4, B5, B6, B7, B8, B9]


""" generate merged datasets for MRE11 """
# GG_mre = alu.load_nparray(ana + "GG-C9_mre11_rs.csv")
# GG_mre_m = alu.mergesubsetcounts(GG_mre, B_GG, ana + "GG-mre-merged.csv")
# TA_mre = alu.load_nparray(ana + "TA-C9_mre11_rs.csv")
# TA_mre_m = alu.mergesubsetcounts(TA_mre, B_TA, ana + "TA-mre-merged.csv")
# CT_mre = alu.load_nparray(ana + "CT-C9_mre11_rs.csv")
# CT_mre_m = alu.mergesubsetcounts(CT_mre, B_CT, ana + "CT-mre-merged.csv")
# alu.mergerows([GG_mre_m, TA_mre_m, CT_mre_m], ana + "ALL-mre-merged.csv")

""" MRE11 output for ML algorithms """
# X, y, alllabels = alu.getXy_all(ana + "ALL-mre-merged.csv")
# X_train, X_test, y_train, y_test = alu.data_split(X, y)
# alu.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_mre11-all.sav",
#                               solver='lbfgs', alpha=0.1, hidden_layer_sizes=(4,))
# alu.ModelTest(X_test, y_test, "baseNN_mre11-all.sav")
# alu.FeatureImportance(X_test, y_test, "baseNN_mre11-all.sav", alllabels)
# alu.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_mre11-all.sav")
# alu.ModelTest(X_test, y_test, "bestNN_mre11-all.sav")

# X, y, nomlabels = alu.getXy_nomismatch(ana + "ALL-mre-merged.csv")
# X_train, X_test, y_train, y_test = alu.data_split(X, y)
# alu.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_mre11-nom.sav",
#                               solver='lbfgs', alpha=2, hidden_layer_sizes=(5,))
# alu.ModelTest(X_test, y_test, "baseNN_mre11-nom.sav")
# alu.FeatureImportance(X_test, y_test, "baseNN_mre11-nom.sav", nomlabels)
# alu.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_mre11-nom.sav")
# alu.ModelTest(X_test, y_test, "bestNN_mre11-nom.sav")

# alu.RandomForestTrainDefault(X_train, y_train, "baseRF_mre11-nom.sav")
# alu.ModelTest(X_test, y_test, "baseRF_mre11-nom.sav")
# alu.RandomForestTrainGridCV(X_train, y_train, "bestRF_mre11-nom.sav")
# alu.ModelTest(X_test, y_test, "bestRF_mre11-nom.sav")


""" generate merged datasets for Cas9 """
# GG_cas = alu.load_nparray(ana + "GG-C9_cas9_rs.csv")
# GG_cas_m = alu.mergesubsetcounts(GG_cas, B_GG, ana + "GG-cas-merged.csv")
# TA_cas = alu.load_nparray(ana + "TA-C9_cas9_rs.csv")
# TA_cas_m = alu.mergesubsetcounts(TA_cas, B_TA, ana + "TA-cas-merged.csv")
# CT_cas = alu.load_nparray(ana + "CT-C9_cas9_rs.csv")
# CT_cas_m = alu.mergesubsetcounts(CT_cas, B_CT, ana + "CT-cas-merged.csv")
# alu.mergerows([GG_cas_m, TA_cas_m, CT_cas_m], ana + "ALL-cas-merged.csv")

""" Cas9 output for  """
# X, y, alllabels = alu.getXy_all(ana + "ALL-cas-merged.csv")
# X_train, X_test, y_train, y_test = alu.data_split(X, y)
# alu.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_cas9-all.sav",
#                               solver='lbfgs', alpha=2, hidden_layer_sizes=(6,))
# alu.ModelTest(X_test, y_test, "baseNN_cas9-all.sav")
# alu.FeatureImportance(X_test, y_test, "baseNN_cas9-all.sav", alllabels)
# alu.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_cas9-all.sav")
# alu.ModelTest(X_test, y_test, "bestNN_cas9-all.sav")

# alu.RandomForestTrainDefault(X_train, y_train, "baseRF_cas9-all.sav")
# alu.ModelTest(X_test, y_test, "baseRF_cas9-all.sav")
# alu.RandomForestTrainGridCV(X_train, y_train, "bestRF_cas9-all.sav")
# alu.ModelTest(X_test, y_test, "bestRF_cas9-all.sav")
