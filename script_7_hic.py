"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.mtss as m
import src.ml as ml
import src.hic as hic
import src.lstm as lstm
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu - Lab computer)
    desktop = "/mnt/c/Users/Roger/Desktop/"
    hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
enc, enc_a = labhome + "public/", labhome + "public/analysis/"
h2WThg19 = labhome + "200212_chipseq_WT1/A18_gh2ax_hg19_final.bam"
bpWThg19 = labhome + "200212_chipseq_WT1/A15_53bp1_hg19_final.bam"
h2GGhg19 = labhome + "200206_chipseq/AluGG-gH2AX_hg19_final.bam"
bpGGhg19 = labhome + "200206_chipseq/AluGG-53BP1_hg19_final.bam"
h2CThg19 = labhome + "200316_chipseq/AluCT-gh2ax-rep1_hg19_final.bam"
bpCThg19 = labhome + "200316_chipseq/AluCT-53bp1-rep1_hg19_final.bam"
h2TAhg19 = labhome + "200316_chipseq/AluTA-gh2ax-rep1_hg19_final.bam"
bpTAhg19 = labhome + "200316_chipseq/AluTA-53bp1-rep1_hg19_final.bam"
bpGG00m_cg = labhome + "201012_chipseq/A04_hg19_final.bam"
bpGG10m_cg = labhome + "201012_chipseq/A14_hg19_final.bam"
bpGG30m_cg = labhome + "201012_chipseq/A15_hg19_final.bam"
h2GG00m_cg = labhome + "201012_chipseq/A16_hg19_final.bam"
h2GG10m_cg = labhome + "201012_chipseq/A17_hg19_final.bam"
h2GG30m_cg = labhome + "201012_chipseq/A18_hg19_final.bam"
mreGGbam = labhome + "200206_chipseq/AluGG-MRE11_hg19_final.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_hg19_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_hg19_final.bam"
iscore_gm12878 = enc + "insulation_scores/GM12878/all_chr_5kb_GM12878.wig"
iscore_huvec = enc + "insulation_scores/HUVEC/all_chr_5kb_HUVEC.wig"
iscore_hmec = enc + "insulation_scores/HMEC/all_chr_5kb_HMEC.wig"
iscore_imr90 = enc + "insulation_scores/IMR90/all_chr_5kb_IMR90.wig"

""" macs2 output """
GG_mre11 = labhome + "200206_chipseq/macs/AluGG-MRE11_hg19_final_peaks.narrowPeak"
CT_mre11 = labhome + "200316_chipseq/macs/AluCT-mre11-rep1_hg19_final_peaks.narrowPeak"
TA_mre11 = labhome + "200316_chipseq/macs/AluTA-mre11-rep1_hg19_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = labhome + "Alu_ana_7_HiC/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_wiggle_windows/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_delta/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_rao_4Cseq/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_insulation_analysis/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
ana_5 = ana + "5_ml_features/"
os.makedirs(ana_5) if not os.path.exists(ana_5) else None


""" ############################################################################################ """
""" Obtain wiggle windows of gH2AX and 53BP1 data """
c.to_wiggle_windows(hg19, h2GGhg19, ana_1 + "GG_gh2ax_hg19", 500)
c.to_wiggle_windows(hg19, bpGGhg19, ana_1 + "GG_53bp1_hg19", 500)
c.to_wiggle_windows(hg19, h2CThg19, ana_1 + "CT_gh2ax_hg19", 500)
c.to_wiggle_windows(hg19, bpCThg19, ana_1 + "CT_53bp1_hg19", 500)
c.to_wiggle_windows(hg19, h2TAhg19, ana_1 + "TA_gh2ax_hg19", 500)
c.to_wiggle_windows(hg19, bpTAhg19, ana_1 + "TA_53bp1_hg19", 500)
c.to_wiggle_windows(hg19, h2WThg19, ana_1 + "WT_gh2ax_hg19", 500)
c.to_wiggle_windows(hg19, bpWThg19, ana_1 + "WT_53bp1_hg19", 500)
c.to_wiggle_windows(hg19, h2GG00m_cg, ana_1 + "GG_cgH2_00m_hg19", 500)
c.to_wiggle_windows(hg19, h2GG10m_cg, ana_1 + "GG_cgH2_10m_hg19", 500)
c.to_wiggle_windows(hg19, h2GG30m_cg, ana_1 + "GG_cgH2_30m_hg19", 500)
c.to_wiggle_windows(hg19, bpGG00m_cg, ana_1 + "GG_cgBP_00m_hg19", 500)
c.to_wiggle_windows(hg19, bpGG10m_cg, ana_1 + "GG_cgBP_10m_hg19", 500)
c.to_wiggle_windows(hg19, bpGG30m_cg, ana_1 + "GG_cgBP_30m_hg19", 500)


""" ############################################################################################ """
""" Obtain absolute changes of gH2AX and 53BP1 wig files """
c.absolutechange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "GG_53bp1_hg19.wig",
                 ana_2 + "GG_53bp1_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "GG_gh2ax_hg19.wig",
                 ana_2 + "GG_gh2ax_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "CT_53bp1_hg19.wig",
                 ana_2 + "CT_53bp1_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "CT_gh2ax_hg19.wig",
                 ana_2 + "CT_gh2ax_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "TA_53bp1_hg19.wig",
                 ana_2 + "TA_53bp1_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "TA_gh2ax_hg19.wig",
                 ana_2 + "TA_gh2ax_00m-3h_hg19_achange")
c.absolutechange(ana_1 + "GG_cgH2_00m_hg19.wig", ana_1 + "GG_cgH2_10m_hg19.wig",
                 ana_2 + "GG_cgH2_00m-10m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgH2_00m_hg19.wig", ana_1 + "GG_cgH2_30m_hg19.wig",
                 ana_2 + "GG_cgH2_00m-30m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgH2_10m_hg19.wig", ana_1 + "GG_cgH2_30m_hg19.wig",
                 ana_2 + "GG_cgH2_10m-30m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgBP_00m_hg19.wig", ana_1 + "GG_cgBP_10m_hg19.wig",
                 ana_2 + "GG_cgBP_00m-10m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgBP_00m_hg19.wig", ana_1 + "GG_cgBP_30m_hg19.wig",
                 ana_2 + "GG_cgBP_00m-30m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgBP_10m_hg19.wig", ana_1 + "GG_cgBP_30m_hg19.wig",
                 ana_2 + "GG_cgBP_10m-30m_hg19_achange")
c.absolutechange(ana_1 + "GG_cgBP_30m_hg19.wig", ana_1 + "GG_53bp1_hg19.wig",
                 ana_2 + "GG_cgBP_30m-3h_hg19_achange")
c.absolutechange(ana_1 + "GG_cgH2_30m_hg19.wig", ana_1 + "GG_gh2ax_hg19.wig",
                 ana_2 + "GG_cgH2_30m-3h_hg19_achange")


""" ############################################################################################ """
""" Obtain 4C-seq profiles from Hi-C data """
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)

gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)


""" ############################################################################################ """
""" Measure relationship between 53BP1/gH2AX enrichment, and insulation scores around cut sites """
# Calculate absolute change in gH2AX (over negative control) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2GGhg19, h2WThg19, ana_4 + "GG_gh2ax_WT-3h_hg19")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2CThg19, h2WThg19, ana_4 + "CT_gh2ax_WT-3h_hg19")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2TAhg19, h2WThg19, ana_4 + "TA_gh2ax_WT-3h_hg19")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "GG_gh2ax_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "CT_gh2ax_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "TA_gh2ax_WT-3h_hg19_achange")
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_merged.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_merged.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_merged.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_merged.csv", head)

# Calculate absolute change in 53BP1 (over negative control) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpGGhg19, bpWThg19, ana_4 + "GG_53bp1_WT-3h_hg19")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpCThg19, bpWThg19, ana_4 + "CT_53bp1_WT-3h_hg19")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpTAhg19, bpWThg19, ana_4 + "TA_53bp1_WT-3h_hg19")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "GG_53bp1_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "CT_53bp1_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "TA_53bp1_WT-3h_hg19_achange")
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_merged.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_merged.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_merged.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_merged.csv", head)

# Calculate derivative of gH2AX relative to each cut site
hic.derivative_from_cutsite(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_delta")
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos_delta.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos_delta.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos_delta.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos_delta.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_pos_delta.csv", head)
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg_delta.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg_delta.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_neg_delta.csv", head)
head = m.load_npheader(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_delta_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_delta_merged.csv"),
             m.load_nparray(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_delta_merged.csv"),
             m.load_nparray(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_delta_merged.csv")],
            ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv", head)

# Calculate derivative of 53BP1 relative to each cut site
hic.derivative_from_cutsite(ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos.csv",
                            ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "GG_53bp1_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "CT_53bp1_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "TA_53bp1_WT-3h_hg19_achange_delta")
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos_delta.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos_delta.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos_delta.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos_delta.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_pos_delta.csv", head)
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg_delta.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg_delta.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg_delta.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg_delta.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_neg_delta.csv", head)
head = m.load_npheader(ana_4 + "GG_53bp1_WT-3h_hg19_achange_delta_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_53bp1_WT-3h_hg19_achange_delta_merged.csv"),
             m.load_nparray(ana_4 + "CT_53bp1_WT-3h_hg19_achange_delta_merged.csv"),
             m.load_nparray(ana_4 + "TA_53bp1_WT-3h_hg19_achange_delta_merged.csv")],
            ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv", head)

# Calculate insulation scores (GM12878) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_4 + "GG_iscore_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_4 + "CT_iscore_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_4 + "TA_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "GG_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "GG_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "CT_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "CT_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "TA_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "TA_iscore_gm12878")
head = m.load_npheader(ana_4 + "GG_iscore_gm12878_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_wigvals_pos.csv")],
            ana_4 + "ALL_iscore_gm12878_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_iscore_gm12878_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_wigvals_neg.csv")],
            ana_4 + "ALL_iscore_gm12878_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_iscore_gm12878_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_merged.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_merged.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_merged.csv")],
            ana_4 + "ALL_iscore_gm12878_merged.csv", head)

# Calculate insulation scores (IMR90) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_4 + "GG_iscore_imr90")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_4 + "CT_iscore_imr90")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_4 + "TA_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "GG_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "GG_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "CT_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "CT_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "TA_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "TA_iscore_imr90")
head = m.load_npheader(ana_4 + "GG_iscore_imr90_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_iscore_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_iscore_imr90_wigvals_pos.csv")],
            ana_4 + "ALL_iscore_imr90_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_iscore_imr90_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_iscore_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_iscore_imr90_wigvals_neg.csv")],
            ana_4 + "ALL_iscore_imr90_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_iscore_imr90_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_imr90_merged.csv"),
             m.load_nparray(ana_4 + "CT_iscore_imr90_merged.csv"),
             m.load_nparray(ana_4 + "TA_iscore_imr90_merged.csv")],
            ana_4 + "ALL_iscore_imr90_merged.csv", head)

# Calculate 4C-seq profiles (GM12878) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "GG_4Cseq_hg19_GM12878.wig", ana_4 + "GG_4Cseq_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "CT_4Cseq_hg19_GM12878.wig", ana_4 + "CT_4Cseq_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "TA_4Cseq_hg19_GM12878.wig", ana_4 + "TA_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "GG_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "GG_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "CT_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "CT_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_4 + "TA_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_4 + "TA_4Cseq_gm12878")
head = m.load_npheader(ana_4 + "GG_4Cseq_gm12878_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_wigvals_pos.csv")],
            ana_4 + "ALL_4Cseq_gm12878_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_4Cseq_gm12878_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_wigvals_neg.csv")],
            ana_4 + "ALL_4Cseq_gm12878_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_4Cseq_gm12878_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_merged.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_merged.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_merged.csv")],
            ana_4 + "ALL_4Cseq_gm12878_merged.csv", head)

# Calculate 4C-seq profiles (IMR90) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "GG_4Cseq_hg19_IMR90.wig", ana_4 + "GG_4Cseq_imr90")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "CT_4Cseq_hg19_IMR90.wig", ana_4 + "CT_4Cseq_imr90")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "TA_4Cseq_hg19_IMR90.wig", ana_4 + "TA_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "GG_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "GG_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "CT_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "CT_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_4 + "TA_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_4 + "TA_4Cseq_imr90")
head = m.load_npheader(ana_4 + "GG_4Cseq_imr90_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_imr90_wigvals_pos.csv")],
            ana_4 + "ALL_4Cseq_imr90_pos.csv", head)
head = m.load_npheader(ana_4 + "GG_4Cseq_imr90_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_imr90_wigvals_neg.csv")],
            ana_4 + "ALL_4Cseq_imr90_neg.csv", head)
head = m.load_npheader(ana_4 + "GG_4Cseq_imr90_merged.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_imr90_merged.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_imr90_merged.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_imr90_merged.csv")],
            ana_4 + "ALL_4Cseq_imr90_merged.csv", head)

# Measure similarity between absolute enrichment change and insulation scores
hic.dtw_randomize(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                  ana_4 + "GG_iscore_gm12878_wigvals_neg.csv",
                  ana_4 + "GG_dtw_gh2ax_gm12878_neg_delta")
hic.categorize_by_insulation_randomize(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                                       ana_4 + "GG_iscore_gm12878_wigvals_neg.csv",
                                       ana_4 + "GG_dtw_gh2ax_gm12878_neg_delta")


""" ############################################################################################ """
""" Obtain features for machine learning model to predict gH2AX span """
# Calculate span width of gH2AX at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.get_span_width(gen, hg19, h2GGhg19, h2WThg19, ana_5 + "GG-gh2ax_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.get_span_width(gen, hg19, h2CThg19, h2WThg19, ana_5 + "CT-gh2ax_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.get_span_width(gen, hg19, h2TAhg19, h2WThg19, ana_5 + "TA-gh2ax_hg19_width")
head = m.load_npheader(ana_5 + "GG-gh2ax_hg19_width.csv")
m.mergerows([m.load_nparray(ana_5 + "GG-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_5 + "CT-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_5 + "TA-gh2ax_hg19_width.csv")],
            ana_5 + "ALL_gh2ax_hg19_width.csv", head)

# Calculate span width of 53BP1 at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.get_span_width(gen, hg19, bpGGhg19, bpWThg19, ana_5 + "GG-53bp1_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.get_span_width(gen, hg19, bpCThg19, bpWThg19, ana_5 + "CT-53bp1_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.get_span_width(gen, hg19, bpTAhg19, bpWThg19, ana_5 + "TA-53bp1_hg19_width")
head = m.load_npheader(ana_5 + "GG-53bp1_hg19_width.csv")
m.mergerows([m.load_nparray(ana_5 + "GG-53bp1_hg19_width.csv"),
             m.load_nparray(ana_5 + "CT-53bp1_hg19_width.csv"),
             m.load_nparray(ana_5 + "TA-53bp1_hg19_width.csv")],
            ana_5 + "ALL_53bp1_hg19_width.csv", head)

# Calculate MRE11 enrichment at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
m.read_counts(gen, mreGGbam, ana_5 + "GG_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
m.read_counts(gen, mreCTbam, ana_5 + "CT_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
m.read_counts(gen, mreTAbam, ana_5 + "TA_mre11_hg19_1250_rc.csv")
head = m.load_npheader(ana_5 + "GG_mre11_hg19_1250_rc.csv")
m.mergerows([m.load_nparray(ana_5 + "GG_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_5 + "CT_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_5 + "TA_mre11_hg19_1250_rc.csv")],
            ana_5 + "ALL_mre11_hg19_1250_rc.csv", head)

# Generate machine learning models for gH2AX
X, y = hic.getXy_insulation(ana_4 + "ALL_iscore_imr90_pos.csv",
                            ana_4 + "ALL_iscore_imr90_neg.csv",
                            ana_5 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_5 + "ALL_gh2ax_hg19_width.csv",
                            ana_5 + "ALL", True)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_5 + "gRF_GG_3h_gh2ax_hg19_1.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gRF_GG_3h_gh2ax_hg19_1.sav")

X, y = hic.getXy_insulation(ana_4 + "ALL_iscore_imr90_pos.csv",
                            ana_4 + "ALL_iscore_imr90_neg.csv",
                            ana_5 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_5 + "ALL_gh2ax_hg19_width.csv",
                            ana_5 + "ALL", False)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_5 + "gRF_GG_3h_gh2ax_hg19_2.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gRF_GG_3h_gh2ax_hg19_2.sav")

lstm.series_to_supervised(ana_4 + 'ALL_gh2ax_WT-3h_hg19_achange_pos.csv',
                          ana_4 + 'ALL_iscore_imr90_pos.csv',
                          ana_4 + "ALL_iscore_imr90_pos.csv")
