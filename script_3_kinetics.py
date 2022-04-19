"""
Script for:
(1) Determining kinetics of dCas9 binding and MRE11 damage responses.
(2) Determining epigenetic enrichment at all target sites.
(3) Developing machine learning models to predict dCas9 binding and MRE11 damage responses.
"""

import src.mtss as m
import src.ml as ml
import numpy as np
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
enc, enc_a = datadir + "public/", datadir + "public/analysis/"
casneg = datadir + "200212_chipseq_WT1/A16_cas9_hg38_final.bam"
mreneg = datadir + "200212_chipseq_WT1/A17_mre11_hg38_final.bam"
casGG30m_1 = datadir + "200804_chipseq/A01_hg38_final.bam"
casGG1h_1 = datadir + "200804_chipseq/A02_hg38_final.bam"
casGG3h_1 = datadir + "200804_chipseq/A03_hg38_final.bam"
mreGG30m_1 = datadir + "200804_chipseq/A04_hg38_final.bam"
mreGG1h_1 = datadir + "200804_chipseq/A05_hg38_final.bam"
mreGG3h_1 = datadir + "200804_chipseq/A06_hg38_final.bam"
casGG30m_2 = datadir + "200804_chipseq/A07_hg38_final.bam"
casGG1h_2 = datadir + "200804_chipseq/A08_hg38_final.bam"
casGG3h_2 = datadir + "200804_chipseq/A09_hg38_final.bam"
mreGG30m_2 = datadir + "200804_chipseq/A10_hg38_final.bam"
mreGG1h_2 = datadir + "200804_chipseq/A11_hg38_final.bam"
mreGG3h_2 = datadir + "200804_chipseq/A12_hg38_final.bam"
h3k4me1_1 = enc + "hg38/H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "hg38/H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "hg38/H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "hg38/H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "hg38/H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "hg38/DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "hg38/ATACseq_HEK293_SRR6418075.bam"              # ATAC-seq (medium deep)
mnase_1 = enc + "hg38/MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "hg38/RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3
pol2_1 = enc + "hg38/POLR2A_HEK293_SRR442119.bam"                # POLR2A ChIP-seq

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
casGG3h_npk2 = datadir + "200804_chipseq/macs/A09_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = datadir + "Alu_ana_3_kinetics/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_epigenetics/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_merged/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_ml/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


""" ############################################################################################ """
""" Get read subsets for time-resolved dCas9 and MRE11 ChIP-seq for AluGG """
m.read_subsets(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8), hg38,
               casneg, ana_1 + "GGk_cas_rs_00m_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8), hg38,
               casGG30m_1, ana_1 + "GGk_cas_rs_30m_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8), hg38,
               casGG1h_1, ana_1 + "GGk_cas_rs_1h_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8), hg38,
               casGG3h_1, ana_1 + "GGk_cas_rs_3h_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               mreneg, ana_1 + "GGk_mre_rs_00m_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG30m_1, ana_1 + "GGk_mre_rs_30m_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG1h_1, ana_1 + "GGk_mre_rs_1h_1")
m.read_subsets(m.macs_gen(casGG3h_npk1, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG3h_1, ana_1 + "GGk_mre_rs_3h_1")

m.read_subsets(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8), hg38,
               casneg, ana_1 + "GGk_cas_rs_00m_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8), hg38,
               casGG30m_2, ana_1 + "GGk_cas_rs_30m_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8), hg38,
               casGG1h_2, ana_1 + "GGk_cas_rs_1h_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8), hg38,
               casGG3h_2, ana_1 + "GGk_cas_rs_3h_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               mreneg, ana_1 + "GGk_mre_rs_00m_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG30m_2, ana_1 + "GGk_mre_rs_30m_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG1h_2, ana_1 + "GGk_mre_rs_1h_2")
m.read_subsets(m.macs_gen(casGG3h_npk2, 1250, hg38, AluGG, fenr=8), hg38,
               mreGG3h_2, ana_1 + "GGk_mre_rs_3h_2")

""" Quantify kinetics for each set """
subset_cas9_1 = [ana_1 + "GGk_cas_rs_00m_1.csv",
                 ana_1 + "GGk_cas_rs_30m_1.csv",
                 ana_1 + "GGk_cas_rs_1h_1.csv",
                 ana_1 + "GGk_cas_rs_3h_1.csv"]
m.read_kinetics(subset_cas9_1, ana_1 + "GGk_cas_rc_kin_1", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_cas9_1, ana_1 + "GGk_cas_rc_kin_1", endname='RefSeq sense', hname='Cend')
subset_mre11_1 = [ana_1 + "GGk_mre_rs_00m_1.csv",
                  ana_1 + "GGk_mre_rs_30m_1.csv",
                  ana_1 + "GGk_mre_rs_1h_1.csv",
                  ana_1 + "GGk_mre_rs_3h_1.csv"]
m.read_kinetics(subset_mre11_1, ana_1 + "GGk_mre_rc_kin_1", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_mre11_1, ana_1 + "GGk_mre_rc_kin_1", endname='RefSeq sense', hname='Cend')
subset_cas9_2 = [ana_1 + "GGk_cas_rs_00m_2.csv",
                 ana_1 + "GGk_cas_rs_30m_2.csv",
                 ana_1 + "GGk_cas_rs_1h_2.csv",
                 ana_1 + "GGk_cas_rs_3h_2.csv"]
m.read_kinetics(subset_cas9_2, ana_1 + "GGk_cas_rc_kin_2", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_cas9_2, ana_1 + "GGk_cas_rc_kin_2", endname='RefSeq sense', hname='Cend')
subset_mre11_2 = [ana_1 + "GGk_mre_rs_00m_2.csv",
                  ana_1 + "GGk_mre_rs_30m_2.csv",
                  ana_1 + "GGk_mre_rs_1h_2.csv",
                  ana_1 + "GGk_mre_rs_3h_2.csv"]
m.read_kinetics(subset_mre11_2, ana_1 + "GGk_mre_rc_kin_2", endname='RefSeq sense', hname='Ctotal')
m.read_kinetics(subset_mre11_2, ana_1 + "GGk_mre_rc_kin_2", endname='RefSeq sense', hname='Cend')


""" ############################################################################################ """
""" Get mismatch and epigenetic information at each cut site """

""" One to Two mismatches annotations """
m.read_mismatch(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8), ana_2 + "GGk_mismatch_1.csv")
m.read_mismatch(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8), ana_2 + "GGk_mismatch_2.csv")

""" ChromHMM epigenetic chromatin state annotation """
m.read_chromhmm(m.macs_gen(casGG3h_npk1, 750, hg38, AluGG, fenr=8),
                hg38, ana_2 + "GGk_chromhmm_1.csv")
m.read_chromhmm(m.macs_gen(casGG3h_npk2, 750, hg38, AluGG, fenr=8),
                hg38, ana_2 + "GGk_chromhmm_2.csv")

""" Get epigenetic info for replicate 1 - AluGG """
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              h3k4me1_1, ana_2 + "GGk_h3k4me1_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              h3k4me3_1, ana_2 + "GGk_h3k4me3_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              h3k9me3_1, ana_2 + "GGk_h3k9me3_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              h3k27ac_1, ana_2 + "GGk_h3k27ac_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              h3k36me3_1, ana_2 + "GGk_h3k36me3_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50, hg38, AluGG, fenr=8),
              dnasei_1, ana_2 + "GGk_dnasei_50_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 10, hg38, AluGG, fenr=8),
              mnase_1, ana_2 + "GGk_mnase_10_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50, hg38, AluGG, fenr=8),
              atac_1, ana_2 + "GGk_atac_50_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              rna_3, ana_2 + "GGk_rna_50000_1.csv")
m.read_counts(m.macs_gen(casGG3h_npk1, 50000, hg38, AluGG, fenr=8),
              pol2_1, ana_2 + "GGk_pol2_50000_1.csv")

""" Get epigenetic info for replicate 2 - AluGG """
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              h3k4me1_1, ana_2 + "GGk_h3k4me1_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              h3k4me3_1, ana_2 + "GGk_h3k4me3_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              h3k9me3_1, ana_2 + "GGk_h3k9me3_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              h3k27ac_1, ana_2 + "GGk_h3k27ac_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              h3k36me3_1, ana_2 + "GGk_h3k36me3_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50, hg38, AluGG, fenr=8),
              dnasei_1, ana_2 + "GGk_dnasei_50_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 10, hg38, AluGG, fenr=8),
              mnase_1, ana_2 + "GGk_mnase_10_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50, hg38, AluGG, fenr=8),
              atac_1, ana_2 + "GGk_atac_50_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              rna_3, ana_2 + "GGk_rna_50000_2.csv")
m.read_counts(m.macs_gen(casGG3h_npk2, 50000, hg38, AluGG, fenr=8),
              pol2_1, ana_2 + "GGk_pol2_50000_2.csv")

""" ############################################################################################ """
""" Merge ChIP-seq and epigenetic enrichment at all dCas9 binding sites """

""" Get merged dataset for replicate 1 """
B01 = m.load_nparray(ana_2 + "GGk_h3k4me1_50000_1.csv")
B02 = m.load_nparray(ana_2 + "GGk_h3k4me3_50000_1.csv")
B03 = m.load_nparray(ana_2 + "GGk_h3k9me3_50000_1.csv")
B04 = m.load_nparray(ana_2 + "GGk_h3k27ac_50000_1.csv")
B05 = m.load_nparray(ana_2 + "GGk_h3k36me3_50000_1.csv")
B06 = m.load_nparray(ana_2 + "GGk_dnasei_50_1.csv")
B07 = m.load_nparray(ana_2 + "GGk_mnase_10_1.csv")
B08 = m.load_nparray(ana_2 + "GGk_atac_50_1.csv")
B09 = m.load_nparray(ana_2 + "GGk_rna_50000_1.csv")
B10 = m.load_nparray(ana_2 + "GGk_pol2_50000_1.csv")
B11 = m.load_nparray(ana_2 + "GGk_chromhmm_1.csv")
B12 = m.load_nparray(ana_2 + "GGk_mismatch_1.csv")
B_GG_final_1 = [B01, B02, B03, B04, B05, B06, B07, B08, B09, B10, B11, B12]
num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, pol2, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"
""" generate merged datasets for MRE11 and Cas9 """
head = m.load_npheader(ana_1 + "GGk_cas_rc_kin_1_Ctotal.csv") + rc_head
GG_cas = m.load_nparray(ana_1 + "GGk_cas_rc_kin_1_Ctotal.csv")
m.mergesubsetcounts(GG_cas, B_GG_final_1, num_cols, ana_3 + "GGk-cas-merged_1_Ctotal.csv", head)
head = m.load_npheader(ana_1 + "GGk_mre_rc_kin_1_Ctotal.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GGk_mre_rc_kin_1_Ctotal.csv")
m.mergesubsetcounts(GG_mre, B_GG_final_1, num_cols, ana_3 + "GGk-mre-merged_1_Ctotal.csv", head)
head = m.load_npheader(ana_1 + "GGk_cas_rc_kin_1_Cend.csv") + rc_head
GG_cas = m.load_nparray(ana_1 + "GGk_cas_rc_kin_1_Cend.csv")
m.mergesubsetcounts(GG_cas, B_GG_final_1, num_cols, ana_3 + "GGk-cas-merged_1_Cend.csv", head)
head = m.load_npheader(ana_1 + "GGk_mre_rc_kin_1_Cend.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GGk_mre_rc_kin_1_Cend.csv")
m.mergesubsetcounts(GG_mre, B_GG_final_1, num_cols, ana_3 + "GGk-mre-merged_1_Cend.csv", head)

""" Get merged dataset for replicate 2 """
B01 = m.load_nparray(ana_2 + "GGk_h3k4me1_50000_2.csv")
B02 = m.load_nparray(ana_2 + "GGk_h3k4me3_50000_2.csv")
B03 = m.load_nparray(ana_2 + "GGk_h3k9me3_50000_2.csv")
B04 = m.load_nparray(ana_2 + "GGk_h3k27ac_50000_2.csv")
B05 = m.load_nparray(ana_2 + "GGk_h3k36me3_50000_2.csv")
B06 = m.load_nparray(ana_2 + "GGk_dnasei_50_2.csv")
B07 = m.load_nparray(ana_2 + "GGk_mnase_10_2.csv")
B08 = m.load_nparray(ana_2 + "GGk_atac_50_2.csv")
B09 = m.load_nparray(ana_2 + "GGk_rna_50000_2.csv")
B10 = m.load_nparray(ana_2 + "GGk_pol2_50000_2.csv")
B11 = m.load_nparray(ana_2 + "GGk_chromhmm_2.csv")
B12 = m.load_nparray(ana_2 + "GGk_mismatch_2.csv")
B_GG_final_2 = [B01, B02, B03, B04, B05, B06, B07, B08, B09, B10, B11, B12]
num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, pol2, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"
""" generate merged datasets for MRE11 and Cas9 """
head = m.load_npheader(ana_1 + "GGk_cas_rc_kin_2_Ctotal.csv") + rc_head
GG_cas = m.load_nparray(ana_1 + "GGk_cas_rc_kin_2_Ctotal.csv")
m.mergesubsetcounts(GG_cas, B_GG_final_2, num_cols, ana_3 + "GGk-cas-merged_2_Ctotal.csv", head)
head = m.load_npheader(ana_1 + "GGk_mre_rc_kin_2_Ctotal.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GGk_mre_rc_kin_2_Ctotal.csv")
m.mergesubsetcounts(GG_mre, B_GG_final_2, num_cols, ana_3 + "GGk-mre-merged_2_Ctotal.csv", head)
head = m.load_npheader(ana_1 + "GGk_cas_rc_kin_2_Cend.csv") + rc_head
GG_cas = m.load_nparray(ana_1 + "GGk_cas_rc_kin_2_Cend.csv")
m.mergesubsetcounts(GG_cas, B_GG_final_2, num_cols, ana_3 + "GGk-cas-merged_2_Cend.csv", head)
head = m.load_npheader(ana_1 + "GGk_mre_rc_kin_2_Cend.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GGk_mre_rc_kin_2_Cend.csv")
m.mergesubsetcounts(GG_mre, B_GG_final_2, num_cols, ana_3 + "GGk-mre-merged_2_Cend.csv", head)

""" Split each time point into separate csv files, with each column corresponding to a different
    value of under column 'mmtype' """
collist = ['timepoint_1', 'timepoint_2', 'timepoint_3', 'timepoint_4']
m.sort_mmtype_column(ana_3 + "GGk-cas-merged_1_Ctotal.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-cas-merged_2_Ctotal.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-mre-merged_1_Ctotal.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-mre-merged_2_Ctotal.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-cas-merged_1_Cend.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-cas-merged_2_Cend.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-mre-merged_1_Cend.csv", collist)
m.sort_mmtype_column(ana_3 + "GGk-mre-merged_2_Cend.csv", collist)

""" ############################################################################################ """
""" #################        MACHINE LEARNING MRE11 AND CAS9 ENRICHMENT        ################# """
mT1 = ana_3 + "GGk-mre-merged_1_Ctotal.csv"
mT2 = ana_3 + "GGk-mre-merged_2_Ctotal.csv"
cT1 = ana_3 + "GGk-cas-merged_1_Ctotal.csv"
cT2 = ana_3 + "GGk-cas-merged_2_Ctotal.csv"
t1 = 'timepoint_1'
t2 = 'timepoint_2'
t4 = 'timepoint_4'
ml_list = [(ml.getXy_2orLess, True, 2, mT1, mT2, t2, ana_4 + "gRF_GG_30T_mre_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mT1, mT2, t2, ana_4 + "gRF_GG_30T_mre_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mT1, mT2, t2, ana_4 + "gRF_GG_30T_mre_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mT1, mT2, t2, ana_4 + "gRF_GG_30T_mre_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mT1, mT2, t2, ana_4 + "gRF_GG_30T_mre_2orless_F1.sav"),
           (ml.getXy_2orLess, True, 2, mT1, mT2, t4, ana_4 + "gRF_GG_3hT_mre_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mT1, mT2, t4, ana_4 + "gRF_GG_3hT_mre_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mT1, mT2, t4, ana_4 + "gRF_GG_3hT_mre_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mT1, mT2, t4, ana_4 + "gRF_GG_3hT_mre_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mT1, mT2, t4, ana_4 + "gRF_GG_3hT_mre_2orless_F1.sav"),
           (ml.getXy_all, True, 2, cT1, cT2, t2, ana_4 + "gRF_GG_30T_cas_all_T2.sav"),
           (ml.getXy_all, True, 1, cT1, cT2, t2, ana_4 + "gRF_GG_30T_cas_all_T1.sav"),
           (ml.getXy_all, True, 0, cT1, cT2, t2, ana_4 + "gRF_GG_30T_cas_all_T0.sav"),
           (ml.getXy_all, False, 2, cT1, cT2, t2, ana_4 + "gRF_GG_30T_cas_all_F2.sav"),
           (ml.getXy_all, False, 1, cT1, cT2, t2, ana_4 + "gRF_GG_30T_cas_all_F1.sav"),
           (ml.getXy_all, True, 2, cT1, cT2, t4, ana_4 + "gRF_GG_3hT_cas_all_T2.sav"),
           (ml.getXy_all, True, 1, cT1, cT2, t4, ana_4 + "gRF_GG_3hT_cas_all_T1.sav"),
           (ml.getXy_all, True, 0, cT1, cT2, t4, ana_4 + "gRF_GG_3hT_cas_all_T0.sav"),
           (ml.getXy_all, False, 2, cT1, cT2, t4, ana_4 + "gRF_GG_3hT_cas_all_F2.sav"),
           (ml.getXy_all, False, 1, cT1, cT2, t4, ana_4 + "gRF_GG_3hT_cas_all_F1.sav"),
           ]
train = True
for fun, epi, mm, inpath1, inpath2, t24, outpath in ml_list:
    X1, y1, labels = fun(inpath1, t24, t1, epi, mm)
    X2, y2, labels2 = fun(inpath2, t24, t1, epi, mm)
    X = np.concatenate((X1, X2))
    y = np.concatenate((y1, y2))
    X_train, X_test, y_train, y_test = ml.data_split(X, y)
    if train:
        ml.RandomForestTrainGridCV(X_train, y_train, outpath)
    ml.ModelTest(X_test, y_test, outpath)
