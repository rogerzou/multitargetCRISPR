"""
Script for:
(1) Analysis of data from Chakrabarti et al

"""

import src.mtss as mtss


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana, enc, enc_a = labhome + "Alu_analysis/", labhome + "public/", labhome + "public/analysis/"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3


""" get all training and testing data from AluGG, AluTA, and AluCT """
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k4me1_1, ana + "chak-C9_h3k4me1_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k4me3_1, ana + "chak-C9_h3k4me3_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k9me3_1, ana + "chak-C9_h3k9me3_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k27ac_1, ana + "chak-C9_h3k27ac_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), h3k36me3_1, ana + "chak-C9_h3k36me3_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(200, hg38), dnasei_1, ana + "chak-C9_dnasei_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(200, hg38), mnase_1, ana + "chak-C9_mnase_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(5000, hg38), atac_1, ana + "chak-C9_atac_rc.csv")
mtss.read_counts(mtss.chakrabarti_generator(2500, hg38), rna_3, ana + "chak-C9_rna_rc.csv")

""" Set arrays for epigenetic data """
B1 = mtss.load_nparray(ana + "chak-C9_h3k4me1_rc.csv")
B2 = mtss.load_nparray(ana + "chak-C9_h3k4me3_rc.csv")
B3 = mtss.load_nparray(ana + "chak-C9_h3k9me3_rc.csv")
B4 = mtss.load_nparray(ana + "chak-C9_h3k27ac_rc.csv")
B5 = mtss.load_nparray(ana + "chak-C9_h3k36me3_rc.csv")
B6 = mtss.load_nparray(ana + "chak-C9_dnasei_rc.csv")
B7 = mtss.load_nparray(ana + "chak-C9_mnase_rc.csv")
B8 = mtss.load_nparray(ana + "chak-C9_atac_rc.csv")
B9 = mtss.load_nparray(ana + "chak-C9_rna_rc.csv")
B_chak = [B1, B2, B3, B4, B5, B6, B7, B8, B9]
mtss.mergesubsetcounts(mtss.load_chakrabarti(), B_chak, ana + "chak-merged.csv")


""" ML models for Chak """
params = {
    'solver': ['lbfgs'],
    'alpha': [0.001, 0.01, 0.1, 1, 2, 4, 5, 6, 8, 10],
    'hidden_layer_sizes': [(4, 4), (8, 4), (16, 8), (10,), (20,), (25,), (30,), (40,), (100,)]
}
X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=0)
X_train, X_test, y_train, y_test = mtss.data_split(X, y, test_size=0.2)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-0.sav", alpha=4, hidden_layer_sizes=(40,))
mtss.ModelTest(X_test, y_test, "baseNN_chak-0-0.sav")
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-0.sav", params=params)
mtss.ModelTest(X_test, y_test, "bestNN_chak-0-0.sav")
mtss.FeatureImportance(X_test, y_test, "baseNN_chak-0-0.sav", alllabels, count=20)

X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=0)
X_train, X_test, y_train, y_test = mtss.data_split(X, y, test_size=0.2)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-0.sav", alpha=3, hidden_layer_sizes=(4,4))
mtss.ModelTest(X_test, y_test, "baseNN_chak-2-0.sav")
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-0.sav", params=params)
mtss.ModelTest(X_test, y_test, "bestNN_chak-2-0.sav")
mtss.FeatureImportance(X_test, y_test, "baseNN_chak-2-0.sav", alllabels, count=20)

X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=1)
X_train, X_test, y_train, y_test = mtss.data_split(X, y)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-0-1.sav", alpha=5, hidden_layer_sizes=(40,))
mtss.ModelTest(X_test, y_test, "baseNN_chak-0-1.sav")
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-1.sav", params=params)
mtss.ModelTest(X_test, y_test, "bestNN_chak-0-1.sav")
mtss.FeatureImportance(X_test, y_test, "bestNN_chak-0-1.sav", alllabels, count=20)

X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=1)
X_train, X_test, y_train, y_test = mtss.data_split(X, y)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "baseNN_chak-2-1.sav", alpha=3, hidden_layer_sizes=(20,))
mtss.ModelTest(X_test, y_test, "baseNN_chak-2-1.sav")
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-1.sav", params=params)
mtss.ModelTest(X_test, y_test, "bestNN_chak-2-1.sav")
mtss.FeatureImportance(X_test, y_test, "bestNN_chak-2-1.sav", alllabels, count=20)

X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=0, index=2)
X_train, X_test, y_train, y_test = mtss.data_split(X, y)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-0-2.sav", alpha=10, hidden_layer_sizes=(4, 4), classifier=True)
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-0-2.sav", params=params, classifier=True)
mtss.ModelTest(X_test, y_test, "bestNN_chak-0-2.sav", classifier=True)

X, y, alllabels = mtss.getXy_chak(ana + "chak-merged.csv", epi=2, index=2)
X_train, X_test, y_train, y_test = mtss.data_split(X, y)
mtss.NeuralNetworkTrainDefault(X_train, y_train, "bestNN_chak-2-2.sav", alpha=10, hidden_layer_sizes=(4, 4), classifier=True)
mtss.NeuralNetworkTrainGridCV(X_train, y_train, "bestNN_chak-2-2.sav", params=params, classifier=True)
mtss.ModelTest(X_test, y_test, "bestNN_chak-2-2.sav", classifier=True)
