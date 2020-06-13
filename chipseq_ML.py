# MH2AX ChIP-seq


import argparse
import src.mtss as mtss
import src.chipseq as c


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana = labhome + "Alu_analysis/"
enc = labhome + "public/"
enc_a = labhome + "public/analysis/"
mreGGin = labhome + "200206_chipseq/16_AluGG-MRE11_final.bam"
casGGin = labhome + "200206_chipseq/15_AluGG-Cas9_final.bam"
mreTAin = labhome + "200316_chipseq/AluTA-mre11-rep1_rmdup.bam"
casTAin = labhome + "200316_chipseq/AluTA-cas9-rep1_rmdup.bam"
mreCTin = labhome + "200316_chipseq/AluCT-mre11-rep1_rmdup.bam"
casCTin = labhome + "200316_chipseq/AluCT-cas9-rep1_rmdup.bam"
mreGGin_noD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_rmdup.bam"
mreGGin_PKi = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_rmdup.bam"
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
# GG_mre11 = labhome + "200206_chipseq/macs/16_AluGG-MRE11_peaks.narrowPeak"
# CT_mre11 = labhome + "200316_chipseq/macs/AluCT-MRE11_peaks.narrowPeak"
# TA_mre11 = labhome + "200316_chipseq/macs/AluTA-MRE11_peaks.narrowPeak"
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


""" Determine multiple sequence alignments from paired-end reads around all sites of enrichment """
# mtss.find_msa(mtss.macs_gen(GG_cas9, 750, hg38, AluGG), casGGin, ana + "GG-C9_cas9_msa", hg38)


""" get all training and testing data from AluGG, AluTA, and AluCT """
# mtss.read_subsets(mtss.macs_gen(GG_cas9, 750, hg38, AluGG), casGGin, ana + "GG-C9_cas9_rs.csv")
# mtss.read_subsets(mtss.macs_gen(TA_cas9, 750, hg38, AluTA), casTAin, ana + "TA-C9_cas9_rs.csv")
# mtss.read_subsets(mtss.macs_gen(CT_cas9, 750, hg38, AluCT), casCTin, ana + "CT-C9_cas9_rs.csv")
# mtss.read_subsets(mtss.macs_gen(GG_cas9, 1250, hg38, AluGG), mreGGin, ana + "GG-C9_mre11_rs.csv")
# mtss.read_subsets(mtss.macs_gen(TA_cas9, 1250, hg38, AluTA), mreTAin, ana + "TA-C9_mre11_rs.csv")
# mtss.read_subsets(mtss.macs_gen(CT_cas9, 1250, hg38, AluCT), mreCTin, ana + "CT-C9_mre11_rs.csv")
# mtss.read_subsets(mtss.macs_gen(GG_cas9, 1250, hg38, AluGG), mreGGin_nD, ana + "GG-C9-noD_mre11_rs.csv")
# mtss.read_subsets(mtss.macs_gen(GG_cas9, 1250, hg38, AluGG), mreGGin_PK, ana + "GG-C9-PKi_mre11_rs.csv")

# mtss.read_counts(mtss.macs_gen(GG_cas9, 10000, hg38, AluGG), h2GGin, ana + "GG-C9_gh2ax_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin, ana + "GG-C9_53bp1_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 10000, hg38, AluTA), h2TAin, ana + "TA-C9_gh2ax_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 10000, hg38, AluTA), bpTAin, ana + "TA-C9_53bp1_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 10000, hg38, AluCT), h2CTin, ana + "CT-C9_gh2ax_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 10000, hg38, AluCT), bpCTin, ana + "CT-C9_53bp1_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin_noD, ana + "GG-C9-noD_53bp1_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 10000, hg38, AluGG), bpGGin_PKi, ana + "GG-C9-PKi_53bp1_rc.csv")

# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), h3k4me1_1, ana + "GG-C9_h3k4me1_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), h3k4me3_1, ana + "GG-C9_h3k4me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), h3k9me3_1, ana + "GG-C9_h3k9me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), h3k27ac_1, ana + "GG-C9_h3k27ac_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), h3k36me3_1, ana + "GG-C9_h3k36me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 200, hg38, AluGG), dnasei_1, ana + "GG-C9_dnasei_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 200, hg38, AluGG), mnase_1, ana + "GG-C9_mnase_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 5000, hg38, AluGG), atac_1, ana + "GG-C9_atac_rc.csv")
# mtss.read_counts(mtss.macs_gen(GG_cas9, 2500, hg38, AluGG), rna_3, ana + "GG-C9_rna_rc.csv")
#
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), h3k4me1_1, ana + "TA-C9_h3k4me1_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), h3k4me3_1, ana + "TA-C9_h3k4me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), h3k9me3_1, ana + "TA-C9_h3k9me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), h3k27ac_1, ana + "TA-C9_h3k27ac_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), h3k36me3_1, ana + "TA-C9_h3k36me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 200, hg38, AluTA), dnasei_1, ana + "TA-C9_dnasei_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 200, hg38, AluTA), mnase_1, ana + "TA-C9_mnase_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 5000, hg38, AluTA), atac_1, ana + "TA-C9_atac_rc.csv")
# mtss.read_counts(mtss.macs_gen(TA_cas9, 2500, hg38, AluTA), rna_3, ana + "TA-C9_rna_rc.csv")

# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), h3k4me1_1, ana + "CT-C9_h3k4me1_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), h3k4me3_1, ana + "CT-C9_h3k4me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), h3k9me3_1, ana + "CT-C9_h3k9me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), h3k27ac_1, ana + "CT-C9_h3k27ac_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), h3k36me3_1, ana + "CT-C9_h3k36me3_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 200, hg38, AluCT), dnasei_1, ana + "CT-C9_dnasei_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 200, hg38, AluCT), mnase_1, ana + "CT-C9_mnase_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 5000, hg38, AluCT), atac_1, ana + "CT-C9_atac_rc.csv")
# mtss.read_counts(mtss.macs_gen(CT_cas9, 2500, hg38, AluCT), rna_3, ana + "CT-C9_rna_rc.csv")
#
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k4me1_1, ana + "chak-C9_h3k4me1_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k4me3_1, ana + "chak-C9_h3k4me3_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k9me3_1, ana + "chak-C9_h3k9me3_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k27ac_1, ana + "chak-C9_h3k27ac_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k36me3_1, ana + "chak-C9_h3k36me3_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(200, hg38), dnasei_1, ana + "chak-C9_dnasei_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(200, hg38), mnase_1, ana + "chak-C9_mnase_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), atac_1, ana + "chak-C9_atac_rc.csv")
# mtss.read_counts(mtss.chakrabarti_generator(2500, hg38), rna_3, ana + "chak-C9_rna_rc.csv")

""" Set arrays for epigenetic data """
B1 = mtss.load_nparray(ana + "GG-C9_h3k4me1_rc.csv")
B2 = mtss.load_nparray(ana + "GG-C9_h3k4me3_rc.csv")
B3 = mtss.load_nparray(ana + "GG-C9_h3k9me3_rc.csv")
B4 = mtss.load_nparray(ana + "GG-C9_h3k27ac_rc.csv")
B5 = mtss.load_nparray(ana + "GG-C9_h3k36me3_rc.csv")
B6 = mtss.load_nparray(ana + "GG-C9_dnasei_rc.csv")
B7 = mtss.load_nparray(ana + "GG-C9_mnase_rc.csv")
B8 = mtss.load_nparray(ana + "GG-C9_atac_rc.csv")
B9 = mtss.load_nparray(ana + "GG-C9_rna_rc.csv")
B_GG = [B1, B2, B3, B4, B5, B6, B7, B8, B9]

B1 = mtss.load_nparray(ana + "TA-C9_h3k4me1_rc.csv")
B2 = mtss.load_nparray(ana + "TA-C9_h3k4me3_rc.csv")
B3 = mtss.load_nparray(ana + "TA-C9_h3k9me3_rc.csv")
B4 = mtss.load_nparray(ana + "TA-C9_h3k27ac_rc.csv")
B5 = mtss.load_nparray(ana + "TA-C9_h3k36me3_rc.csv")
B6 = mtss.load_nparray(ana + "TA-C9_dnasei_rc.csv")
B7 = mtss.load_nparray(ana + "TA-C9_mnase_rc.csv")
B8 = mtss.load_nparray(ana + "TA-C9_atac_rc.csv")
B9 = mtss.load_nparray(ana + "TA-C9_rna_rc.csv")
B_TA = [B1, B2, B3, B4, B5, B6, B7, B8, B9]

B1 = mtss.load_nparray(ana + "CT-C9_h3k4me1_rc.csv")
B2 = mtss.load_nparray(ana + "CT-C9_h3k4me3_rc.csv")
B3 = mtss.load_nparray(ana + "CT-C9_h3k9me3_rc.csv")
B4 = mtss.load_nparray(ana + "CT-C9_h3k27ac_rc.csv")
B5 = mtss.load_nparray(ana + "CT-C9_h3k36me3_rc.csv")
B6 = mtss.load_nparray(ana + "CT-C9_dnasei_rc.csv")
B7 = mtss.load_nparray(ana + "CT-C9_mnase_rc.csv")
B8 = mtss.load_nparray(ana + "CT-C9_atac_rc.csv")
B9 = mtss.load_nparray(ana + "CT-C9_rna_rc.csv")
B_CT = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
#
# B1 = mtss.load_nparray(ana + "chak-C9_h3k4me1_rc.csv")
# B2 = mtss.load_nparray(ana + "chak-C9_h3k4me3_rc.csv")
# B3 = mtss.load_nparray(ana + "chak-C9_h3k9me3_rc.csv")
# B4 = mtss.load_nparray(ana + "chak-C9_h3k27ac_rc.csv")
# B5 = mtss.load_nparray(ana + "chak-C9_h3k36me3_rc.csv")
# B6 = mtss.load_nparray(ana + "chak-C9_dnasei_rc.csv")
# B7 = mtss.load_nparray(ana + "chak-C9_mnase_rc.csv")
# B8 = mtss.load_nparray(ana + "chak-C9_atac_rc.csv")
# B9 = mtss.load_nparray(ana + "chak-C9_rna_rc.csv")
# B_chak = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
# mtss.mergesubsetcounts(mtss.load_chakrabarti(), B_chak, ana + "chak-merged.csv")

""" generate merged datasets for MRE11 """
GG_mre = mtss.load_nparray(ana + "GG-C9_mre11_rs.csv")
GG_mre_m = mtss.mergesubsetcounts(GG_mre, B_GG, ana + "GG-mre-merged.csv")
TA_mre = mtss.load_nparray(ana + "TA-C9_mre11_rs.csv")
TA_mre_m = mtss.mergesubsetcounts(TA_mre, B_TA, ana + "TA-mre-merged.csv")
CT_mre = mtss.load_nparray(ana + "CT-C9_mre11_rs.csv")
CT_mre_m = mtss.mergesubsetcounts(CT_mre, B_CT, ana + "CT-mre-merged.csv")
mtss.mergerows([GG_mre_m, TA_mre_m, CT_mre_m], ana + "ALL-mre-merged.csv")

""" generate merged datasets for Cas9 """
GG_cas = mtss.load_nparray(ana + "GG-C9_cas9_rs.csv")
GG_cas_m = mtss.mergesubsetcounts(GG_cas, B_GG, ana + "GG-cas-merged.csv")
TA_cas = mtss.load_nparray(ana + "TA-C9_cas9_rs.csv")
TA_cas_m = mtss.mergesubsetcounts(TA_cas, B_TA, ana + "TA-cas-merged.csv")
CT_cas = mtss.load_nparray(ana + "CT-C9_cas9_rs.csv")
CT_cas_m = mtss.mergesubsetcounts(CT_cas, B_CT, ana + "CT-cas-merged.csv")
mtss.mergerows([GG_cas_m, TA_cas_m, CT_cas_m], ana + "ALL-cas-merged.csv")



""" MRE11 output for ML algorithms """
# X, y, alllabels = mtss.getXy_all(ana + "ALL-mre-merged.csv")
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_mre11-all.sav",
#                               solver='lbfgs', alpha=0.1, hidden_layer_sizes=(4,))
# mtss.ModelTest(X_test, y_test, "baseNN_mre11-all.sav")
# mtss.FeatureImportance(X_test, y_test, "baseNN_mre11-all.sav", alllabels)
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_mre11-all.sav")
# mtss.ModelTest(X_test, y_test, "bestNN_mre11-all.sav")

# X, y, nomlabels = mtss.getXy_nomismatch(ana + "ALL-mre-merged.csv")
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_mre11-nom.sav",
#                               solver='lbfgs', alpha=2, hidden_layer_sizes=(5,))
# mtss.ModelTest(X_test, y_test, "baseNN_mre11-nom.sav")
# mtss.FeatureImportance(X_test, y_test, "baseNN_mre11-nom.sav", nomlabels)
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_mre11-nom.sav")
# mtss.ModelTest(X_test, y_test, "bestNN_mre11-nom.sav")

# mtss.RandomForestTrainDefault(X_train, y_train, "baseRF_mre11-nom.sav")
# mtss.ModelTest(X_test, y_test, "baseRF_mre11-nom.sav")
# mtss.RandomForestTrainGridCV(X_train, y_train, "bestRF_mre11-nom.sav")
# mtss.ModelTest(X_test, y_test, "bestRF_mre11-nom.sav")


""" Cas9 output for  """
# X, y, alllabels = mtss.getXy_all(ana + "ALL-cas-merged.csv")
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_cas9-all.sav",
#                               solver='lbfgs', alpha=2, hidden_layer_sizes=(6,))
# mtss.ModelTest(X_test, y_test, "baseNN_cas9-all.sav")
# mtss.FeatureImportance(X_test, y_test, "baseNN_cas9-all.sav", alllabels)
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_cas9-all.sav")
# mtss.ModelTest(X_test, y_test, "bestNN_cas9-all.sav")

# mtss.RandomForestTrainDefault(X_train, y_train, "baseRF_cas9-all.sav")
# mtss.ModelTest(X_test, y_test, "baseRF_cas9-all.sav")
# mtss.RandomForestTrainGridCV(X_train, y_train, "bestRF_cas9-all.sav")
# mtss.ModelTest(X_test, y_test, "bestRF_cas9-all.sav")


""" ML models for Chak """
params = {
    'solver': ['lbfgs'],
    'alpha': [0.001, 0.01, 0.1, 1, 2, 4, 5, 6, 8, 10],
    'hidden_layer_sizes': [(4, 4), (8, 4), (16, 8), (10,), (20,), (25,), (30,), (40,), (100,)]
}
# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=0)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y, test_size=0.2)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-0.sav", alpha=4, hidden_layer_sizes=(40,))
# mtss.ModelTest(X_test, y_test, "baseNN_chak-0-0.sav")
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-0.sav", params=params)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-0-0.sav")
# mtss.FeatureImportance(X_test, y_test, "baseNN_chak-0-0.sav", alllabels, count=20)

# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=0)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y, test_size=0.2)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-0.sav", alpha=3, hidden_layer_sizes=(4,4))
# mtss.ModelTest(X_test, y_test, "baseNN_chak-2-0.sav")
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-0.sav", params=params)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-2-0.sav")
# mtss.FeatureImportance(X_test, y_test, "baseNN_chak-2-0.sav", alllabels, count=20)

# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=1)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-1.sav", alpha=5, hidden_layer_sizes=(40,))
# mtss.ModelTest(X_test, y_test, "baseNN_chak-0-1.sav")
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-1.sav", params=params)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-0-1.sav")
# mtss.FeatureImportance(X_test, y_test, "bestNN_chak-0-1.sav", alllabels, count=20)

# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=1)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-1.sav", alpha=3, hidden_layer_sizes=(20,))
# mtss.ModelTest(X_test, y_test, "baseNN_chak-2-1.sav")
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-1.sav", params=params)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-2-1.sav")
# mtss.FeatureImportance(X_test, y_test, "bestNN_chak-2-1.sav", alllabels, count=20)

# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=2)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-0-2.sav", alpha=10, hidden_layer_sizes=(4, 4), classifier=True)
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-2.sav", params=params, classifier=True)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-0-2.sav", classifier=True)

# X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=2)
# X_train, X_test, y_train, y_test = mtss.data_split(X, y)
# mtss.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-2-2.sav", alpha=10, hidden_layer_sizes=(4, 4), classifier=True)
# mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-2.sav", params=params, classifier=True)
# mtss.ModelTest(X_test, y_test, "bestNN_chak-2-2.sav", classifier=True)
