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
    desktop = "/mnt/c/Users/rzou4/Desktop/"
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
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
atac_1 = enc + "hg38/ATACseq_HEK293_SRR6418075.bam"              # ATAC-seq (medium deep)
mnase_1 = enc + "hg38/MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_1 = enc + "hg38/RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3
pol2_1 = enc + "hg38/POLR2A_HEK293_SRR442119.bam"                # POLR2A ChIP-seq

""" Set analysis path """
ana = labhome + "Alu_ana_11_chak/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_counts/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None


""" ############################################################################################ """
""" get all  """
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), h3k4me1_1, ana_1 + "chak-C9_h3k4me1_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), h3k4me3_1, ana_1 + "chak-C9_h3k4me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), h3k9me3_1, ana_1 + "chak-C9_h3k9me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), h3k27ac_1, ana_1 + "chak-C9_h3k27ac_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), h3k36me3_1, ana_1 + "chak-C9_h3k36me3_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50, hg38), dnasei_1, ana_1 + "chak-C9_dnasei_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50, hg38), mnase_1, ana_1 + "chak-C9_mnase_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), atac_1, ana_1 + "chak-C9_atac_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), rna_1, ana_1 + "chak-C9_rna_rc.csv")
mtss.read_counts(chak.chakrabarti_generator(50000, hg38), pol2_1, ana_1 + "chak-C9_pol2_rc.csv")
