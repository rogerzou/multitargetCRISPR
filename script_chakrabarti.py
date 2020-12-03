"""
Script for:
(1) Analysis of data from Chakrabarti et al

"""

import src.mtss as mtss
import src.chak as chak
import src.ml as ml
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
h3k4me1_1 = enc + "hg38/H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "hg38/H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "hg38/H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "hg38/H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "hg38/H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "hg38/DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "hg38/ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq
mnase_1 = enc + "hg38/MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_1 = enc + "hg38/RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

ana = labhome + "Alu_ana_10_chak/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_counts/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None

""" get all training and testing data from AluGG, AluTA, and AluCT """
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), h3k4me1_1, ana_1 + "chak-C9_h3k4me1_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), h3k4me3_1, ana_1 + "chak-C9_h3k4me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), h3k9me3_1, ana_1 + "chak-C9_h3k9me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), h3k27ac_1, ana_1 + "chak-C9_h3k27ac_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), h3k36me3_1, ana_1 + "chak-C9_h3k36me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50, hg38[1]), dnasei_1, ana_1 + "chak-C9_dnasei_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50, hg38[1]), mnase_1, ana_1 + "chak-C9_mnase_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), atac_1, ana_1 + "chak-C9_atac_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38[1]), rna_1, ana_1 + "chak-C9_rna_rc.csv")

""" Set arrays for epigenetic data """
B1 = mtss.load_nparray(ana_1 + "chak-C9_h3k4me1_rc.csv")
B2 = mtss.load_nparray(ana_1 + "chak-C9_h3k4me3_rc.csv")
B3 = mtss.load_nparray(ana_1 + "chak-C9_h3k9me3_rc.csv")
B4 = mtss.load_nparray(ana_1 + "chak-C9_h3k27ac_rc.csv")
B5 = mtss.load_nparray(ana_1 + "chak-C9_h3k36me3_rc.csv")
B6 = mtss.load_nparray(ana_1 + "chak-C9_dnasei_rc.csv")
B7 = mtss.load_nparray(ana_1 + "chak-C9_mnase_rc.csv")
B8 = mtss.load_nparray(ana_1 + "chak-C9_atac_rc.csv")
B9 = mtss.load_nparray(ana_1 + "chak-C9_rna_rc.csv")
B_chak = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna"
head = mtss.load_npheader(ana_1 + "GGk_mre_rc_kin_1_Ctotal.csv") + rc_head
mtss.mergesubsetcounts(chak.load_chakrabarti(), B_chak, num_cols, ana_1 + "chak-merged.csv", head)


""" ML models for Chak """
params = {
    'solver': ['lbfgs'],
    'alpha': [0.001, 0.01, 0.1, 1, 2, 4, 5, 6, 8, 10],
    'hidden_layer_sizes': [(4, 4), (8, 4), (16, 8), (10,), (20,), (25,), (30,), (40,), (100,)]
}
X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=0, index=0)
X_train, X_test, y_train, y_test = ml.data_split(X, y, test_size=0.2)
ml.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-0.sav")
ml.ModelTest(X_test, y_test, "baseNN_chak-0-0.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-0.sav")
ml.ModelTest(X_test, y_test, "bestNN_chak-0-0.sav")
ml.FeatureImportance(X_test, y_test, "baseNN_chak-0-0.sav", alllabels, count=20)

X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=2, index=0)
X_train, X_test, y_train, y_test = ml.data_split(X, y, test_size=0.2)
ml.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-0.sav")
ml.ModelTest(X_test, y_test, "baseNN_chak-2-0.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-0.sav")
ml.ModelTest(X_test, y_test, "bestNN_chak-2-0.sav")
ml.FeatureImportance(X_test, y_test, "baseNN_chak-2-0.sav", alllabels, count=20)

X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=0, index=1)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-1.sav")
ml.ModelTest(X_test, y_test, "baseNN_chak-0-1.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-1.sav")
ml.ModelTest(X_test, y_test, "bestNN_chak-0-1.sav")
ml.FeatureImportance(X_test, y_test, "bestNN_chak-0-1.sav", alllabels, count=20)

X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=2, index=1)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-1.sav")
ml.ModelTest(X_test, y_test, "baseNN_chak-2-1.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-1.sav")
ml.ModelTest(X_test, y_test, "bestNN_chak-2-1.sav")
ml.FeatureImportance(X_test, y_test, "bestNN_chak-2-1.sav", alllabels, count=20)

X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=0, index=2)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-0-2.sav", classifier=True)
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-2.sav", classifier=True)
ml.ModelTest(X_test, y_test, "bestNN_chak-0-2.sav", classifier=True)

X, y, alllabels = chak.getXy_chak(ana_1 + "chak-merged.csv", epi=2, index=2)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-2-2.sav", classifier=True)
ml.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-2.sav", classifier=True)
ml.ModelTest(X_test, y_test, "bestNN_chak-2-2.sav", classifier=True)
