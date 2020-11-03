"""
Script for:
(1) Determination of time-resolved MRE11, 53BP1, and gH2AX enrichment after activation of cgRNA.
"""

import src.mtss as m
import src.chipseq as c
import src.hic as hic
import src.ml as ml
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu - Lab computer)
    desktop = "/mnt/c/Users/Roger/Desktop/"
    hg38 = "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
h2WTin = labhome + "200212_chipseq_WT1/A18_WT1H_final.bam"
bpWTin = labhome + "200212_chipseq_WT1/A15_WT1B_final.bam"
mreGG00m_cg = labhome + "201012_chipseq/A01_final.bam"
mreGG10m_cg = labhome + "201012_chipseq/A02_final.bam"
mreGG30m_cg = labhome + "201012_chipseq/A03_final.bam"
bpGG00m_cg = labhome + "201012_chipseq/A04_final.bam"
bpGG10m_cg = labhome + "201012_chipseq/A14_final.bam"
bpGG30m_cg = labhome + "201012_chipseq/A15_final.bam"
h2GG00m_cg = labhome + "201012_chipseq/A16_final.bam"
h2GG10m_cg = labhome + "201012_chipseq/A17_final.bam"
h2GG30m_cg = labhome + "201012_chipseq/A18_final.bam"
epi = labhome + "Alu_ana_5_timeresolved/2_epigenetics/"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/AluGG-dCas9_3h_1_new_peaks.narrowPeak"
GG_mre11 = labhome + "200206_chipseq/macs/16_AluGG-MRE11_peaks.narrowPeak"

""" Set analysis path """
ana = labhome + "Alu_ana_6_cgRNA/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_merged/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_ml/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_wiggle_windows/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
ana_5 = ana + "5_span_profiles/"
os.makedirs(ana_5) if not os.path.exists(ana_5) else None


""" ############################################################################################ """
""" Get read subsets for vfCRISPR time-resolved MRE11 ChIP-seq for AluGG """
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG00m_cg,
               ana_1 + "GG_cgMRE_00m_rs_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG10m_cg,
               ana_1 + "GG_cgMRE_10m_rs_1")
m.read_subsets(m.macs_gen(GG3h_npk1, 1250, hg38, AluGG, fenr=8), mreGG30m_cg,
               ana_1 + "GG_cgMRE_30m_rs_1")
subset_mre11_1 = [ana_1 + "GG_cgMRE_00m_rs_1.csv",
                  ana_1 + "GG_cgMRE_10m_rs_1.csv",
                  ana_1 + "GG_cgMRE_30m_rs_1.csv"]
m.read_kinetics(subset_mre11_1, ana_1 + "GG_cgMRE_kin_1", 'Ctotal')
m.read_kinetics(subset_mre11_1, ana_1 + "GG_cgMRE_kin_1", 'Cend')


""" ############################################################################################ """
""" Get merged datasets for machine learning """
""" Set arrays for final epigenetic data """
B1 = m.load_nparray(epi +"GGk_h3k4me1_50000_1.csv")
B2 = m.load_nparray(epi +"GGk_h3k4me3_50000_1.csv")
B3 = m.load_nparray(epi +"GGk_h3k9me3_50000_1.csv")
B4 = m.load_nparray(epi +"GGk_h3k27ac_50000_1.csv")
B5 = m.load_nparray(epi +"GGk_h3k36me3_50000_1.csv")
B6 = m.load_nparray(epi +"GGk_dnasei_50_1.csv")
B7 = m.load_nparray(epi +"GGk_mnase_50_1.csv")
B8 = m.load_nparray(epi +"GGk_atac_50000_1.csv")
B9 = m.load_nparray(epi +"GGk_rna_50000_1.csv")
B10 = m.load_nparray(epi +"GGk_chromhmm_1.csv")
B11 = m.load_nparray(epi +"GGk_mismatch_1.csv")
B_GG_1 = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
num_cols = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3]
rc_head = ", h3k4me1, h3k4me3, h3k9me3, h3k27ac, h3k36me3, dnasei, mnase, atac, rna, " \
          "chromhmm, chromhmm_active, mm_pos1, mm_pos2, mm_type"
""" generate merged datasets for MRE11 and Cas9 """
head = m.load_npheader(ana_1 + "GG_cgMRE_kin_1_Ctotal.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GG_cgMRE_kin_1_Ctotal.csv")
m.mergesubsetcounts(GG_mre, B_GG_1, num_cols, ana_2 + "GG_cgMRE_merged_1_Ctotal.csv", head)
head = m.load_npheader(ana_1 + "GG_cgMRE_kin_1_Cend.csv") + rc_head
GG_mre = m.load_nparray(ana_1 + "GG_cgMRE_kin_1_Cend.csv")
m.mergesubsetcounts(GG_mre, B_GG_1, num_cols, ana_2 + "GG_cgMRE_merged_1_Cend.csv", head)

""" Split each time point into separate csv files, with each column corresponding to a different
    value of under column 'mmtype' """
collist = ['timepoint_1', 'timepoint_2', 'timepoint_3']
m.sort_mmtype_column(ana_2 + "GG_cgMRE_merged_1_Ctotal.csv", collist)
m.sort_mmtype_column(ana_2 + "GG_cgMRE_merged_1_Cend.csv", collist)


""" #################                 RANDOM FOREST MRE11                   ################# """
mT = ana_2 + "GG_cgMRE_merged_1_Ctotal.csv"
mE = ana_2 + "GG_cgMRE_merged_1_Cend.csv"
t1 = 'timepoint_1'
t3 = 'timepoint_3'
ml_list = [(ml.getXy_2orLess, True, 2, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_2orless_F1.sav"),
           (ml.getXy_2orLess, True, 2, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_2orless_F1.sav"),
           (ml.getXy_2orLess, True, 2, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_2orless_F1.sav"),
           (ml.getXy_2orLess, True, 2, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_2orless_T2.sav"),
           (ml.getXy_2orLess, True, 1, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_2orless_T1.sav"),
           (ml.getXy_2orLess, True, 0, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_2orless_T0.sav"),
           (ml.getXy_2orLess, False, 2, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_2orless_F2.sav"),
           (ml.getXy_2orLess, False, 1, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_2orless_F1.sav"),
           (ml.getXy_noMM, True, 2, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_noMM_T2.sav"),
           (ml.getXy_noMM, True, 1, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_noMM_T1.sav"),
           (ml.getXy_noMM, True, 0, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_noMM_T0.sav"),
           (ml.getXy_noMM, False, 2, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_noMM_F2.sav"),
           (ml.getXy_noMM, False, 1, mT, t1, ana_3 + "gRF_GG_cgMRE_10mT_noMM_F1.sav"),
           (ml.getXy_noMM, True, 2, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_noMM_T2.sav"),
           (ml.getXy_noMM, True, 1, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_noMM_T1.sav"),
           (ml.getXy_noMM, True, 0, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_noMM_T0.sav"),
           (ml.getXy_noMM, False, 2, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_noMM_F2.sav"),
           (ml.getXy_noMM, False, 1, mE, t1, ana_3 + "gRF_GG_cgMRE_10mE_noMM_F1.sav"),
           (ml.getXy_noMM, True, 2, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_noMM_T2.sav"),
           (ml.getXy_noMM, True, 1, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_noMM_T1.sav"),
           (ml.getXy_noMM, True, 0, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_noMM_T0.sav"),
           (ml.getXy_noMM, False, 2, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_noMM_F2.sav"),
           (ml.getXy_noMM, False, 1, mT, t3, ana_3 + "gRF_GG_cgMRE_30mT_noMM_F1.sav"),
           (ml.getXy_noMM, True, 2, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_noMM_T2.sav"),
           (ml.getXy_noMM, True, 1, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_noMM_T1.sav"),
           (ml.getXy_noMM, True, 0, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_noMM_T0.sav"),
           (ml.getXy_noMM, False, 2, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_noMM_F2.sav"),
           (ml.getXy_noMM, False, 1, mE, t3, ana_3 + "gRF_GG_cgMRE_30mE_noMM_F1.sav")

           ]
train = False
for fun, epi, mm, inpath, t13, outpath in ml_list:
    X, y, labels = fun(inpath, t13, epi, mm)
    X_train, X_test, y_train, y_test = ml.data_split(X, y)
    if train:
        ml.RandomForestTrainGridCV(X_train, y_train, outpath)
    ml.ModelTest(X_test, y_test, outpath)
    # ml.FeatureImportance(X_test, y_test, outpath, labels, count=20)


""" ############################################################################################ """
""" Get vfCRISPR AluGG 53BP1 and gH2AX ChIP-seq averaged enrichment as wiggle files. """
c.to_wiggle_windows(h2GG00m_cg, ana_4 + "GG_cgH2_00m_1", 500)
c.to_wiggle_windows(h2GG10m_cg, ana_4 + "GG_cgH2_10m_1", 500)
c.to_wiggle_windows(h2GG30m_cg, ana_4 + "GG_cgH2_30m_1", 500)
c.to_wiggle_windows(bpGG00m_cg, ana_4 + "GG_cgBP_00m_1", 500)
c.to_wiggle_windows(bpGG10m_cg, ana_4 + "GG_cgBP_10m_1", 500)
c.to_wiggle_windows(bpGG30m_cg, ana_4 + "GG_cgBP_30m_1", 500)


""" ############################################################################################ """
""" Get 53BP1 and gH2AX span and profiles """
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, h2GG10m_cg, h2GG00m_cg, ana_5 + "GG_cgH2_10m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, h2GG30m_cg, h2GG00m_cg, ana_5 + "GG_cgH2_30m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, bpGG10m_cg, bpGG00m_cg, ana_5 + "GG_cgBP_10m_1_span")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg38, AluGG, fenr=8), 2000000)
hic.get_span_width(gen, bpGG30m_cg, bpGG00m_cg, ana_5 + "GG_cgBP_30m_1_span")
