"""
Script for:
(1) Analysis of H2AX and 53BP1 data in the context of 3D genome architecture

"""

import src.chipseq as c
import src.mtss as m
import src.ml as ml
import src.hic as hic
import src.lstm as lstm
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
enc, enc_a = datadir + "public/", datadir + "public/analysis/"
h2WThg19 = datadir + "200212_chipseq_WT1/A18_gh2ax_hg19_final.bam"
bpWThg19 = datadir + "200212_chipseq_WT1/A15_53bp1_hg19_final.bam"
h2GGhg19 = datadir + "200206_chipseq/AluGG-gH2AX_hg19_final.bam"
bpGGhg19 = datadir + "200206_chipseq/AluGG-53BP1_hg19_final.bam"
h2CThg19 = datadir + "200316_chipseq/AluCT-gh2ax-rep1_hg19_final.bam"
bpCThg19 = datadir + "200316_chipseq/AluCT-53bp1-rep1_hg19_final.bam"
h2TAhg19 = datadir + "200316_chipseq/AluTA-gh2ax-rep1_hg19_final.bam"
bpTAhg19 = datadir + "200316_chipseq/AluTA-53bp1-rep1_hg19_final.bam"
bpGG00m_cg = datadir + "201012_chipseq/A04_hg19_final.bam"
bpGG10m_cg = datadir + "201012_chipseq/A14_hg19_final.bam"
bpGG30m_cg = datadir + "201012_chipseq/A15_hg19_final.bam"
h2GG00m_cg = datadir + "201012_chipseq/A16_hg19_final.bam"
h2GG10m_cg = datadir + "201012_chipseq/A17_hg19_final.bam"
h2GG30m_cg = datadir + "201012_chipseq/A18_hg19_final.bam"
mreGGbam = datadir + "200206_chipseq/AluGG-MRE11_hg19_final.bam"
mreCTbam = datadir + "200316_chipseq/AluCT-mre11-rep1_hg19_final.bam"
mreTAbam = datadir + "200316_chipseq/AluTA-mre11-rep1_hg19_final.bam"
iscore_gm12878 = enc + "insulation_scores/GM12878/all_chr_5kb_GM12878.wig"
iscore_hmec = enc + "insulation_scores/HMEC/all_chr_5kb_HMEC.wig"
iscore_huvec = enc + "insulation_scores/HUVEC/all_chr_5kb_HUVEC.wig"
iscore_imr90 = enc + "insulation_scores/IMR90/all_chr_5kb_IMR90.wig"
iscore_k562 = enc + "insulation_scores/K562/all_chr_5kb_K562.wig"
iscore_kbm7 = enc + "insulation_scores/KBM7/all_chr_5kb_KBM7.wig"
iscore_nhek = enc + "insulation_scores/NHEK/all_chr_5kb_NHEK.wig"

""" macs2 output """
GG_mre11 = datadir + "200206_chipseq/macs/AluGG-MRE11_hg19_final_peaks.narrowPeak"
CT_mre11 = datadir + "200316_chipseq/macs/AluCT-mre11-rep1_hg19_final_peaks.narrowPeak"
TA_mre11 = datadir + "200316_chipseq/macs/AluTA-mre11-rep1_hg19_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = datadir + "Alu_ana_4_HiC/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_wiggle_windows/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_delta/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_rao_4Cseq/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_chipseq_from_cutsite_wtRNA/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
ana_5 = ana + "5_chipseq_from_cutsite_cgRNA/"
os.makedirs(ana_5) if not os.path.exists(ana_5) else None
ana_6 = ana + "6_insulation_from_cutsite/"
os.makedirs(ana_6) if not os.path.exists(ana_6) else None
ana_7 = ana + "7_4Cseq_from_cutsite/"
os.makedirs(ana_7) if not os.path.exists(ana_7) else None
ana_8 = ana + "8_matplotlib/"
os.makedirs(ana_8) if not os.path.exists(ana_8) else None


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
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_GM12878", datadir + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_GM12878", datadir + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_GM12878", datadir + "public_HiC/GM12878", 5, 2E6)

gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_IMR90", datadir + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_IMR90", datadir + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_IMR90", datadir + "public_HiC/IMR90", 5, 2E6)


""" ############################################################################################ """
""" Measure 53BP1/gH2AX enrichment relative to cut sites for GG, CG, and TA datasets """
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


""" ############################################################################################ """
""" Measure 53BP1/gH2AX enrichment relative to cut sites for cgRNA datasets """
# Calculate absolute change in gH2AX (AluGG cgRNA over no-light) (100k window)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2GG10m_cg, h2GG00m_cg, ana_5 + "GG_gh2ax_00m-10m_hg19")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2GG30m_cg, h2GG10m_cg, ana_5 + "GG_gh2ax_10m-30m_hg19")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2GGhg19, h2GG30m_cg, ana_5 + "GG_gh2ax_30m-3h_hg19")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_gh2ax_00m-10m_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_gh2ax_10m-30m_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_gh2ax_30m-3h_hg19_achange")

# Calculate absolute change in 53BP1 (AluGG cgRNA over no-light) (100k window)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpGG10m_cg, bpGG00m_cg, ana_5 + "GG_53bp1_00m-10m_hg19")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpGG30m_cg, bpGG10m_cg, ana_5 + "GG_53bp1_10m-30m_hg19")
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpGGhg19, bpGG30m_cg, ana_5 + "GG_53bp1_30m-3h_hg19")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_53bp1_00m-10m_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_53bp1_00m-10m_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_53bp1_00m-10m_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_53bp1_10m-30m_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_53bp1_10m-30m_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_53bp1_10m-30m_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_53bp1_30m-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_5 + "GG_53bp1_30m-3h_hg19_achange_neg.csv",
                        outpath=ana_5 + "GG_53bp1_30m-3h_hg19_achange")

# Calculate derivative of gH2AX cgRNA relative to each cut site
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_delta")

# Calculate derivative of 53BP1 cgRNA relative to each cut site
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos.csv",
                            ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg")
hic.derivative_from_cutsite(ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg.csv",
                            ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_5 + "GG_gh2ax_30m-3h_hg19_achange_delta")


""" ############################################################################################ """
""" Obtain insulation scores for different cell types relative to cut sites """
# Calculate insulation scores (GM12878) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_6 + "GG_iscore_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_6 + "CT_iscore_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_6 + "TA_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_gm12878_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_gm12878")
head = m.load_npheader(ana_6 + "GG_iscore_gm12878_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_gm12878_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_gm12878_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_gm12878_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_gm12878_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_gm12878_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_gm12878_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_gm12878_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_gm12878_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_gm12878_merged.csv")],
            ana_6 + "ALL_iscore_gm12878_merged.csv", head)

# Calculate insulation scores (HMEC) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_hmec, ana_6 + "GG_iscore_hmec")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_hmec, ana_6 + "CT_iscore_hmec")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_hmec, ana_6 + "TA_iscore_hmec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_hmec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_hmec_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_hmec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_hmec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_hmec_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_hmec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_hmec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_hmec_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_hmec")
head = m.load_npheader(ana_6 + "GG_iscore_hmec_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_hmec_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_hmec_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_hmec_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_hmec_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_hmec_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_hmec_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_hmec_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_hmec_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_hmec_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_hmec_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_hmec_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_hmec_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_hmec_merged.csv")],
            ana_6 + "ALL_iscore_hmec_merged.csv", head)

# Calculate insulation scores (HUVEC) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_huvec, ana_6 + "GG_iscore_huvec")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_huvec, ana_6 + "CT_iscore_huvec")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_huvec, ana_6 + "TA_iscore_huvec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_huvec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_huvec_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_huvec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_huvec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_huvec_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_huvec")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_huvec_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_huvec_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_huvec")
head = m.load_npheader(ana_6 + "GG_iscore_huvec_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_huvec_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_huvec_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_huvec_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_huvec_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_huvec_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_huvec_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_huvec_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_huvec_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_huvec_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_huvec_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_huvec_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_huvec_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_huvec_merged.csv")],
            ana_6 + "ALL_iscore_huvec_merged.csv", head)

# Calculate insulation scores (IMR90) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_6 + "GG_iscore_imr90")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_6 + "CT_iscore_imr90")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_imr90, ana_6 + "TA_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_imr90")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_imr90_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_imr90_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_imr90")
head = m.load_npheader(ana_6 + "GG_iscore_imr90_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_imr90_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_imr90_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_imr90_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_imr90_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_imr90_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_imr90_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_imr90_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_imr90_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_imr90_merged.csv")],
            ana_6 + "ALL_iscore_imr90_merged.csv", head)

# Calculate insulation scores (K562) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_k562, ana_6 + "GG_iscore_k562")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_k562, ana_6 + "CT_iscore_k562")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_k562, ana_6 + "TA_iscore_k562")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_k562_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_k562_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_k562")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_k562_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_k562_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_k562")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_k562_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_k562_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_k562")
head = m.load_npheader(ana_6 + "GG_iscore_k562_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_k562_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_k562_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_k562_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_k562_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_k562_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_k562_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_k562_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_k562_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_k562_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_k562_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_k562_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_k562_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_k562_merged.csv")],
            ana_6 + "ALL_iscore_k562_merged.csv", head)

# Calculate insulation scores (NHEK) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_kbm7, ana_6 + "GG_iscore_kbm7")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_kbm7, ana_6 + "CT_iscore_kbm7")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_kbm7, ana_6 + "TA_iscore_kbm7")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_kbm7_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_kbm7_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_kbm7")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_kbm7_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_kbm7_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_kbm7")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_kbm7_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_kbm7_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_kbm7")
head = m.load_npheader(ana_6 + "GG_iscore_kbm7_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_kbm7_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_kbm7_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_kbm7_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_kbm7_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_kbm7_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_kbm7_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_kbm7_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_kbm7_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_kbm7_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_kbm7_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_kbm7_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_kbm7_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_kbm7_merged.csv")],
            ana_6 + "ALL_iscore_kbm7_merged.csv", head)

# Calculate insulation scores (NHEK) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_nhek, ana_6 + "GG_iscore_nhek")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_nhek, ana_6 + "CT_iscore_nhek")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_nhek, ana_6 + "TA_iscore_nhek")
hic.merged_from_cutsite(f_data_pos=ana_6 + "GG_iscore_nhek_wigvals_pos.csv",
                        f_data_neg=ana_6 + "GG_iscore_nhek_wigvals_neg.csv",
                        outpath=ana_6 + "GG_iscore_nhek")
hic.merged_from_cutsite(f_data_pos=ana_6 + "CT_iscore_nhek_wigvals_pos.csv",
                        f_data_neg=ana_6 + "CT_iscore_nhek_wigvals_neg.csv",
                        outpath=ana_6 + "CT_iscore_nhek")
hic.merged_from_cutsite(f_data_pos=ana_6 + "TA_iscore_nhek_wigvals_pos.csv",
                        f_data_neg=ana_6 + "TA_iscore_nhek_wigvals_neg.csv",
                        outpath=ana_6 + "TA_iscore_nhek")
head = m.load_npheader(ana_6 + "GG_iscore_nhek_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_nhek_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "CT_iscore_nhek_wigvals_pos.csv"),
             m.load_nparray(ana_6 + "TA_iscore_nhek_wigvals_pos.csv")],
            ana_6 + "ALL_iscore_nhek_pos.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_nhek_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_nhek_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "CT_iscore_nhek_wigvals_neg.csv"),
             m.load_nparray(ana_6 + "TA_iscore_nhek_wigvals_neg.csv")],
            ana_6 + "ALL_iscore_nhek_neg.csv", head)
head = m.load_npheader(ana_6 + "GG_iscore_nhek_merged.csv")
m.mergerows([m.load_nparray(ana_6 + "GG_iscore_nhek_merged.csv"),
             m.load_nparray(ana_6 + "CT_iscore_nhek_merged.csv"),
             m.load_nparray(ana_6 + "TA_iscore_nhek_merged.csv")],
            ana_6 + "ALL_iscore_nhek_merged.csv", head)

""" Determine enrichment values at either >0 or <=0 insulation scores for WT vs 3h gRNAs"""
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_gm12878_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_gm12878_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_hmec_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_hmec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_huvec_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_huvec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_imr90_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_imr90_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_k562_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_k562_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_kbm7_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_kbm7_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_nhek_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_WT-3h_nhek_delta_merged")

""" Determine enrichment values at either >0 or <=0 insulation scores for 00m, 10m, 30m cgRNAs"""
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_gm12878_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_gm12878_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_hmec_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_hmec_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_huvec_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_huvec_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_imr90_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_imr90_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_k562_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_k562_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_kbm7_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_kbm7_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_00m-10m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_nhek_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_00m-10m_nhek_delta_merged")

""" Determine enrichment values at either >0 or <=0 insulation scores for 00m, 10m, 30m cgRNAs"""
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_gm12878_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_gm12878_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_hmec_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_hmec_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_huvec_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_huvec_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_imr90_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_imr90_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_k562_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_k562_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_kbm7_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_kbm7_delta_merged")
hic.categorize_by_insulation_randomize(ana_5 + "GG_gh2ax_10m-30m_hg19_achange_delta_merged.csv",
                                       ana_6 + "GG_iscore_nhek_merged.csv",
                                       ana_6 + "GG_categorize_gh2ax_10m-30m_nhek_delta_merged")


""" ############################################################################################ """
""" Calculate 4C-seq profiles relative to cut sites """
# Calculate 4C-seq profiles (GM12878) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "GG_4Cseq_hg19_GM12878.wig", ana_7 + "GG_4Cseq_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "CT_4Cseq_hg19_GM12878.wig", ana_7 + "CT_4Cseq_gm12878")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "TA_4Cseq_hg19_GM12878.wig", ana_7 + "TA_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_7 + "GG_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_7 + "GG_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_7 + "GG_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_7 + "CT_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_7 + "CT_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_7 + "CT_4Cseq_gm12878")
hic.merged_from_cutsite(f_data_pos=ana_7 + "TA_4Cseq_gm12878_wigvals_pos.csv",
                        f_data_neg=ana_7 + "TA_4Cseq_gm12878_wigvals_neg.csv",
                        outpath=ana_7 + "TA_4Cseq_gm12878")
head = m.load_npheader(ana_7 + "GG_4Cseq_gm12878_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_gm12878_wigvals_pos.csv")],
            ana_7 + "ALL_4Cseq_gm12878_pos.csv", head)
head = m.load_npheader(ana_7 + "GG_4Cseq_gm12878_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_gm12878_wigvals_neg.csv")],
            ana_7 + "ALL_4Cseq_gm12878_neg.csv", head)
head = m.load_npheader(ana_7 + "GG_4Cseq_gm12878_merged.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_gm12878_merged.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_gm12878_merged.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_gm12878_merged.csv")],
            ana_7 + "ALL_4Cseq_gm12878_merged.csv", head)

# Calculate 4C-seq profiles (IMR90) relative to each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "GG_4Cseq_hg19_IMR90.wig", ana_7 + "GG_4Cseq_imr90")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "CT_4Cseq_hg19_IMR90.wig", ana_7 + "CT_4Cseq_imr90")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "TA_4Cseq_hg19_IMR90.wig", ana_7 + "TA_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_7 + "GG_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_7 + "GG_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_7 + "GG_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_7 + "CT_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_7 + "CT_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_7 + "CT_4Cseq_imr90")
hic.merged_from_cutsite(f_data_pos=ana_7 + "TA_4Cseq_imr90_wigvals_pos.csv",
                        f_data_neg=ana_7 + "TA_4Cseq_imr90_wigvals_neg.csv",
                        outpath=ana_7 + "TA_4Cseq_imr90")
head = m.load_npheader(ana_7 + "GG_4Cseq_imr90_wigvals_pos.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_imr90_wigvals_pos.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_imr90_wigvals_pos.csv")],
            ana_7 + "ALL_4Cseq_imr90_pos.csv", head)
head = m.load_npheader(ana_7 + "GG_4Cseq_imr90_wigvals_neg.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_imr90_wigvals_neg.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_imr90_wigvals_neg.csv")],
            ana_7 + "ALL_4Cseq_imr90_neg.csv", head)
head = m.load_npheader(ana_7 + "GG_4Cseq_imr90_merged.csv")
m.mergerows([m.load_nparray(ana_7 + "GG_4Cseq_imr90_merged.csv"),
             m.load_nparray(ana_7 + "CT_4Cseq_imr90_merged.csv"),
             m.load_nparray(ana_7 + "TA_4Cseq_imr90_merged.csv")],
            ana_7 + "ALL_4Cseq_imr90_merged.csv", head)


""" ############################################################################################ """
""" Save plots of enrichment and insulation values relative to cut site using matplotlib """
hic.print_merged_from_cutsite(f_data1=ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_merged.csv",
                              f_data2=ana_6 + "ALL_iscore_gm12878_merged.csv",
                              outpath=ana_8 + "ALL_print_WT-3h_gh2ax_gm12878_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_merged.csv",
                              f_data2=ana_6 + "GG_iscore_gm12878_merged.csv",
                              outpath=ana_8 + "GG_print_00m-10m_gh2ax_gm12878_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_merged.csv",
                              f_data2=ana_6 + "GG_iscore_gm12878_merged.csv",
                              outpath=ana_8 + "GG_print_10m-30m_gh2ax_gm12878_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_merged.csv",
                              f_data2=ana_4 + "ALL_53bp1_WT-3h_hg19_achange_merged.csv",
                              outpath=ana_8 + "ALL_print_WT-3h_gh2ax_53bp1_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_5 + "GG_gh2ax_00m-10m_hg19_achange_merged.csv",
                              f_data2=ana_5 + "GG_53bp1_00m-10m_hg19_achange_merged.csv",
                              outpath=ana_8 + "GG_print_00m-10m_gh2ax_53bp1_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_5 + "GG_gh2ax_10m-30m_hg19_achange_merged.csv",
                              f_data2=ana_5 + "GG_53bp1_10m-30m_hg19_achange_merged.csv",
                              outpath=ana_8 + "GG_print_10m-30m_gh2ax_53bp1_merged", res=5000)
hic.print_merged_from_cutsite(f_data1=ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_merged.csv",
                              f_data2=ana_5 + 'ALL_h3k9me3_bamvals_merged.csv',
                              outpath=ana_8 + "ALL_print_WT-3h_gh2ax_h3k9me3_merged", res=5000)
