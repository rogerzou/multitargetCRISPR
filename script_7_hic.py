"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.msa as msa
import src.mtss as m
import src.ml as ml
import src.hic as hic
import sys
import os
import numpy as np

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
alnpath = labhome + "Alu_ana_1_putative/1_protosearch/psearch_align.csv"

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
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)

gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)

gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)

gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)


""" ############################################################################################ """
""" Measure relationship between 53BP1/gH2AX enrichment, and insulation scores around cut sites """
# Calculate absolute change in enrichment (over negative control) relative to each cut site
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2GGhg19, h2WThg19, ana_4 + "GG_gh2ax_WT-3h_hg19")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpGGhg19, bpWThg19, ana_4 + "GG_53bp1_WT-3h_hg19")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2CThg19, h2WThg19, ana_4 + "CT_gh2ax_WT-3h_hg19")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpCThg19, bpWThg19, ana_4 + "CT_53bp1_WT-3h_hg19")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, h2TAhg19, h2WThg19, ana_4 + "TA_gh2ax_WT-3h_hg19")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.absolute_change_from_cutsite(gen, hg19, bpTAhg19, bpWThg19, ana_4 + "TA_53bp1_WT-3h_hg19")

hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "GG_gh2ax_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "GG_53bp1_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "CT_gh2ax_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "CT_53bp1_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "TA_gh2ax_WT-3h_hg19_achange")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos.csv",
                        f_data_neg=ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg.csv",
                        outpath=ana_4 + "TA_53bp1_WT-3h_hg19_achange")
hic.derivative_from_cutsite(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg")
hic.derivative_from_cutsite(ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg.csv",
                            ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "GG_gh2ax_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "GG_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "GG_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "GG_53bp1_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "CT_gh2ax_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "CT_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "CT_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "CT_53bp1_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "TA_gh2ax_WT-3h_hg19_achange_delta")
hic.merged_from_cutsite(f_data_pos=ana_4 + "TA_53bp1_WT-3h_hg19_achange_pos_delta.csv",
                        f_data_neg=ana_4 + "TA_53bp1_WT-3h_hg19_achange_neg_delta.csv",
                        outpath=ana_4 + "TA_53bp1_WT-3h_hg19_achange_delta")

# Calculate insulation scores relative to each cut site
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_4 + "GG_iscore_gm12878")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, iscore_gm12878, ana_4 + "CT_iscore_gm12878")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
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
head = m.load_npheader(ana_4 + "GG_iscore_gm12878_insulation_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_wigvals_pos.csv")],
            ana_4 + "ALL_iscore_gm12878_pos.csv", head)
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_wigvals_neg.csv")],
            ana_4 + "ALL_iscore_gm12878_neg.csv", head)
m.mergerows([m.load_nparray(ana_4 + "GG_iscore_gm12878_wigvals_merged.csv"),
             m.load_nparray(ana_4 + "CT_iscore_gm12878_wigvals_merged.csv"),
             m.load_nparray(ana_4 + "TA_iscore_gm12878_wigvals_merged.csv")],
            ana_4 + "ALL_iscore_gm12878_merged.csv", head)

# Calculate 4C-seq profiles relative to each cut site
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "GG_4Cseq_hg19_GM12878.wig", ana_4 + "GG_4Cseq_gm12878")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.wigvals_from_cutsite(gen, hg19, ana_3 + "CT_4Cseq_hg19_GM12878.wig", ana_4 + "CT_4Cseq_gm12878")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
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
head = m.load_npheader(ana_4 + "GG_4Cseq_gm12878_pos.csv")
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_wigvals_pos.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_wigvals_pos.csv")],
            ana_4 + "ALL_4Cseq_gm12878_pos.csv", head)
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_wigvals_neg.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_wigvals_neg.csv")],
            ana_4 + "ALL_4Cseq_gm12878_neg.csv", head)
m.mergerows([m.load_nparray(ana_4 + "GG_4Cseq_gm12878_wigvals_merged.csv"),
             m.load_nparray(ana_4 + "CT_4Cseq_gm12878_wigvals_merged.csv"),
             m.load_nparray(ana_4 + "TA_4Cseq_gm12878_wigvals_merged.csv")],
            ana_4 + "ALL_4Cseq_gm12878_merged.csv", head)

# Measure similarity between absolute enrichment change and insulation scores
hic.dtw_randomize(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                  ana_4 + "GG_iscore_gm12878_insulation_neg.csv",
                  ana_4 + "GG_dtw_gh2ax_gm12878_neg_delta")
hic.categorize_by_insulation_randomize(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                                       ana_4 + "GG_iscore_gm12878_insulation_neg.csv",
                                       ana_4 + "GG_dtw_gh2ax_gm12878_neg_delta")


""" ############################################################################################ """
""" Obtain features for machine learning model to predict gH2AX span """
# Calculate span width at each cut site
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.get_span_width(gen, hg19, h2GGhg19, h2WThg19, ana_5 + "GG-gh2ax_hg19_width")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.get_span_width(gen, hg19, h2CThg19, h2WThg19, ana_5 + "CT-gh2ax_hg19_width")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.get_span_width(gen, hg19, h2TAhg19, h2WThg19, ana_5 + "TA-gh2ax_hg19_width")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
hic.get_span_width(gen, hg19, bpGGhg19, bpWThg19, ana_5 + "GG-53bp1_hg19_width")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
hic.get_span_width(gen, hg19, bpCThg19, bpWThg19, ana_5 + "CT-53bp1_hg19_width")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
hic.get_span_width(gen, hg19, bpTAhg19, bpWThg19, ana_5 + "TA-53bp1_hg19_width")
head = m.load_npheader(ana_5 + "GG-gh2ax_hg19_width.csv")
m.mergerows([m.load_nparray(ana_5 + "GG-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_5 + "CT-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_5 + "TA-gh2ax_hg19_width.csv")],
            ana_5 + "ALL_gh2ax_hg19_width.csv", head)
head = m.load_npheader(ana_5 + "GG-53bp1_hg19_width.csv")
m.mergerows([m.load_nparray(ana_5 + "GG-53bp1_hg19_width.csv"),
             m.load_nparray(ana_5 + "CT-53bp1_hg19_width.csv"),
             m.load_nparray(ana_5 + "TA-53bp1_hg19_width.csv")],
            ana_5 + "ALL_53bp1_hg19_width.csv", head)

# Calculate MRE11 enrichment at each cut site
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluGG), 2E6)
m.read_counts(gen, mreGGbam, ana_5 + "GG_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluCT), 2E6)
m.read_counts(gen, mreCTbam, ana_5 + "CT_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 750, AluTA), 2E6)
m.read_counts(gen, mreTAbam, ana_5 + "TA_mre11_hg19_1250_rc.csv")
head = m.load_npheader(ana_5 + "GG_mre11_hg19_1250_rc.csv")
m.mergerows([m.load_nparray(ana_5 + "GG_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_5 + "CT_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_5 + "TA_mre11_hg19_1250_rc.csv")],
            ana_5 + "ALL_mre11_hg19_1250_rc.csv", head)

X, y = hic.getXy_insulation(ana_4 + "ALL_iscore_gm12878_insulation_pos.csv",
                            ana_4 + "ALL_iscore_gm12878_insulation_neg.csv",
                            ana_5 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_5 + "ALL_gh2ax_hg19_width.csv",
                            ana_5 + "ALL_features", True)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_5 + "gRF_GG_3h_gh2ax_hg19_1.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gRF_GG_3h_gh2ax_hg19_1.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, ana_5 + "gNN_GG_3h_gh2ax_hg19_1.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gNN_GG_3h_gh2ax_hg19_1.sav")

X, y = hic.getXy_insulation(ana_4 + "ALL_iscore_gm12878_insulation_pos.csv",
                            ana_4 + "ALL_iscore_gm12878_insulation_neg.csv",
                            ana_5 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_5 + "ALL_gh2ax_hg19_width.csv",
                            ana_5 + "ALL_features", False)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_5 + "gRF_GG_3h_gh2ax_hg19_2.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gRF_GG_3h_gh2ax_hg19_2.sav")
ml.NeuralNetworkTrainGridCV(X_train, y_train, ana_5 + "gNN_GG_3h_gh2ax_hg19_2.sav")
ml.ModelTest(X_test, y_test, ana_5 + "gNN_GG_3h_gh2ax_hg19_2.sav")
