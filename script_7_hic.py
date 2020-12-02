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
iscore_hmec = enc + "insulation_scores/HMEC/all_chr_5kb_HMEC.wig"
iscore_huvec = enc + "insulation_scores/HUVEC/all_chr_5kb_HUVEC.wig"
iscore_imr90 = enc + "insulation_scores/IMR90/all_chr_5kb_IMR90.wig"
iscore_k562 = enc + "insulation_scores/K562/all_chr_5kb_K562.wig"
iscore_kbm7 = enc + "insulation_scores/KBM7/all_chr_5kb_KBM7.wig"
iscore_nhek = enc + "insulation_scores/NHEK/all_chr_5kb_NHEK.wig"
h3k4me1_1 = enc + "H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "ATACseq_HEK293_SRR6418075.bam"              # ATAC-seq (medium deep)
# atac_2 = enc + "ATACseq_HEK293_SRR5627157.bam"              # ATAC-seq (bad)
mnase_1 = enc + "MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3
ctcf_1 = enc + "CTCF_HEK293_ENCFF354GWS.bam"                # CTCF
smc3_1 = enc + "SMC3_GM12878_ENCFF375OQR.bam"               # SMC3

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
ana_9 = ana + "9_ml_features/"
os.makedirs(ana_9) if not os.path.exists(ana_9) else None


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
""" Measure epigenetic enrichment relative to cut sites for GG, CG, and TA datasets """
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me1_1, ana_5 + "GG_h3k4me1", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me3_1, ana_5 + "GG_h3k4me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k9me3_1, ana_5 + "GG_h3k9me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k27ac_1, ana_5 + "GG_h3k27ac", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k36me3_1, ana_5 + "GG_h3k36me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, atac_1, ana_5 + "GG_atacseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, rna_3, ana_5 + "GG_rnaseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, ctcf_1, ana_5 + "GG_ctcf", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.bamvals_from_cutsite(gen, hg19, smc3_1, ana_5 + "GG_smc3", span_rad=10000)
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_h3k4me1_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_h3k4me1_bamvals_neg.csv",
                        outpath=ana_5 + "GG_h3k4me1_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_h3k4me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_h3k4me3_bamvals_neg.csv",
                        outpath=ana_5 + "GG_h3k4me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_h3k9me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_h3k9me3_bamvals_neg.csv",
                        outpath=ana_5 + "GG_h3k9me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_h3k27ac_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_h3k27ac_bamvals_neg.csv",
                        outpath=ana_5 + "GG_h3k27ac_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_h3k36me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_h3k36me3_bamvals_neg.csv",
                        outpath=ana_5 + "GG_h3k36me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_atacseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_atacseq_bamvals_neg.csv",
                        outpath=ana_5 + "GG_atacseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_rnaseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_rnaseq_bamvals_neg.csv",
                        outpath=ana_5 + "GG_rnaseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_ctcf_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_ctcf_bamvals_neg.csv",
                        outpath=ana_5 + "GG_ctcf_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "GG_smc3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "GG_smc3_bamvals_neg.csv",
                        outpath=ana_5 + "GG_smc3_bamvals")

gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me1_1, ana_5 + "CT_h3k4me1", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me3_1, ana_5 + "CT_h3k4me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k9me3_1, ana_5 + "CT_h3k9me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k27ac_1, ana_5 + "CT_h3k27ac", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k36me3_1, ana_5 + "CT_h3k36me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, atac_1, ana_5 + "CT_atacseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, rna_3, ana_5 + "CT_rnaseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, ctcf_1, ana_5 + "CT_ctcf", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.bamvals_from_cutsite(gen, hg19, smc3_1, ana_5 + "CT_smc3", span_rad=10000)
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_h3k4me1_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_h3k4me1_bamvals_neg.csv",
                        outpath=ana_5 + "CT_h3k4me1_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_h3k4me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_h3k4me3_bamvals_neg.csv",
                        outpath=ana_5 + "CT_h3k4me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_h3k9me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_h3k9me3_bamvals_neg.csv",
                        outpath=ana_5 + "CT_h3k9me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_h3k27ac_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_h3k27ac_bamvals_neg.csv",
                        outpath=ana_5 + "CT_h3k27ac_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_h3k36me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_h3k36me3_bamvals_neg.csv",
                        outpath=ana_5 + "CT_h3k36me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_atacseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_atacseq_bamvals_neg.csv",
                        outpath=ana_5 + "CT_atacseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_rnaseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_rnaseq_bamvals_neg.csv",
                        outpath=ana_5 + "CT_rnaseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_ctcf_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_ctcf_bamvals_neg.csv",
                        outpath=ana_5 + "CT_ctcf_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "CT_smc3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "CT_smc3_bamvals_neg.csv",
                        outpath=ana_5 + "CT_smc3_bamvals")

gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me1_1, ana_5 + "TA_h3k4me1", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k4me3_1, ana_5 + "TA_h3k4me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k9me3_1, ana_5 + "TA_h3k9me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k27ac_1, ana_5 + "TA_h3k27ac", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, h3k36me3_1, ana_5 + "TA_h3k36me3", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, atac_1, ana_5 + "TA_atacseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, rna_3, ana_5 + "TA_rnaseq", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, ctcf_1, ana_5 + "TA_ctcf", span_rad=10000)
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.bamvals_from_cutsite(gen, hg19, smc3_1, ana_5 + "TA_smc3", span_rad=10000)
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_h3k4me1_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_h3k4me1_bamvals_neg.csv",
                        outpath=ana_5 + "TA_h3k4me1_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_h3k4me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_h3k4me3_bamvals_neg.csv",
                        outpath=ana_5 + "TA_h3k4me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_h3k9me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_h3k9me3_bamvals_neg.csv",
                        outpath=ana_5 + "TA_h3k9me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_h3k27ac_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_h3k27ac_bamvals_neg.csv",
                        outpath=ana_5 + "TA_h3k27ac_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_h3k36me3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_h3k36me3_bamvals_neg.csv",
                        outpath=ana_5 + "TA_h3k36me3_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_atacseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_atacseq_bamvals_neg.csv",
                        outpath=ana_5 + "TA_atacseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_rnaseq_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_rnaseq_bamvals_neg.csv",
                        outpath=ana_5 + "TA_rnaseq_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_ctcf_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_ctcf_bamvals_neg.csv",
                        outpath=ana_5 + "TA_ctcf_bamvals")
hic.merged_from_cutsite(f_data_pos=ana_5 + "TA_smc3_bamvals_pos.csv",
                        f_data_neg=ana_5 + "TA_smc3_bamvals_neg.csv",
                        outpath=ana_5 + "TA_smc3_bamvals")

head = m.load_npheader(ana_5 + "GG_h3k4me1_bamvals_pos.csv")
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me1_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me1_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me1_bamvals_pos.csv")],
            ana_5 + "ALL_h3k4me1_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me3_bamvals_pos.csv")],
            ana_5 + "ALL_h3k4me3_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k9me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_h3k9me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_h3k9me3_bamvals_pos.csv")],
            ana_5 + "ALL_h3k9me3_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k27ac_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_h3k27ac_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_h3k27ac_bamvals_pos.csv")],
            ana_5 + "ALL_h3k27ac_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k36me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_h3k36me3_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_h3k36me3_bamvals_pos.csv")],
            ana_5 + "ALL_h3k36me3_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_atacseq_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_atacseq_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_atacseq_bamvals_pos.csv")],
            ana_5 + "ALL_atacseq_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_rnaseq_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_rnaseq_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_rnaseq_bamvals_pos.csv")],
            ana_5 + "ALL_rnaseq_bamvals_pos.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_ctcf_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "CT_ctcf_bamvals_pos.csv"),
             m.load_nparray(ana_5 + "TA_ctcf_bamvals_pos.csv")],
            ana_5 + "ALL_ctcf_bamvals_pos.csv", head)

head = m.load_npheader(ana_5 + "GG_h3k4me1_bamvals_neg.csv")
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me1_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me1_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me1_bamvals_neg.csv")],
            ana_5 + "ALL_h3k4me1_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me3_bamvals_neg.csv")],
            ana_5 + "ALL_h3k4me3_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k9me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_h3k9me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_h3k9me3_bamvals_neg.csv")],
            ana_5 + "ALL_h3k9me3_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k27ac_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_h3k27ac_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_h3k27ac_bamvals_neg.csv")],
            ana_5 + "ALL_h3k27ac_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k36me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_h3k36me3_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_h3k36me3_bamvals_neg.csv")],
            ana_5 + "ALL_h3k36me3_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_atacseq_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_atacseq_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_atacseq_bamvals_neg.csv")],
            ana_5 + "ALL_atacseq_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_rnaseq_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_rnaseq_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_rnaseq_bamvals_neg.csv")],
            ana_5 + "ALL_rnaseq_bamvals_neg.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_ctcf_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "CT_ctcf_bamvals_neg.csv"),
             m.load_nparray(ana_5 + "TA_ctcf_bamvals_neg.csv")],
            ana_5 + "ALL_ctcf_bamvals_neg.csv", head)

head = m.load_npheader(ana_5 + "GG_h3k4me1_bamvals_merged.csv")
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me1_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me1_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me1_bamvals_merged.csv")],
            ana_5 + "ALL_h3k4me1_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k4me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_h3k4me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_h3k4me3_bamvals_merged.csv")],
            ana_5 + "ALL_h3k4me3_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k9me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_h3k9me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_h3k9me3_bamvals_merged.csv")],
            ana_5 + "ALL_h3k9me3_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k27ac_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_h3k27ac_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_h3k27ac_bamvals_merged.csv")],
            ana_5 + "ALL_h3k27ac_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_h3k36me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_h3k36me3_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_h3k36me3_bamvals_merged.csv")],
            ana_5 + "ALL_h3k36me3_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_atacseq_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_atacseq_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_atacseq_bamvals_merged.csv")],
            ana_5 + "ALL_atacseq_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_rnaseq_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_rnaseq_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_rnaseq_bamvals_merged.csv")],
            ana_5 + "ALL_rnaseq_bamvals_merged.csv", head)
m.mergerows([m.load_nparray(ana_5 + "GG_ctcf_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "CT_ctcf_bamvals_merged.csv"),
             m.load_nparray(ana_5 + "TA_ctcf_bamvals_merged.csv")],
            ana_5 + "ALL_ctcf_bamvals_merged.csv", head)


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

""" Determine enrichment values at either >0 or <=0 insulation scores """
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_gm12878_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_gm12878_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_hmec_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_hmec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_huvec_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_huvec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_imr90_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_imr90_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_k562_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_k562_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_kbm7_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_kbm7_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_nhek_merged.csv",
                                       ana_6 + "ALL_categorize_gh2ax_nhek_delta_merged")

hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_gm12878_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_gm12878_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_hmec_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_hmec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_huvec_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_huvec_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_imr90_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_imr90_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_k562_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_k562_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_kbm7_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_kbm7_delta_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_53bp1_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_nhek_merged.csv",
                                       ana_6 + "ALL_categorize_53bp1_nhek_delta_merged")

# Measure similarity between absolute enrichment change and insulation scores
hic.dtw_randomize(ana_4 + "GG_gh2ax_WT-3h_hg19_achange_neg_delta.csv",
                  ana_6 + "GG_iscore_gm12878_wigvals_neg.csv",
                  ana_6 + "GG_dtw_gh2ax_gm12878_neg_delta")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_merged.csv",
                                       ana_6 + "ALL_iscore_gm12878_merged.csv",
                                       ana_6 + "GG_catins_gh2ax_gm12878_merged")
hic.categorize_by_insulation_randomize(ana_4 + "ALL_gh2ax_WT-3h_hg19_achange_delta_merged.csv",
                                       ana_6 + "ALL_iscore_gm12878_merged.csv",
                                       ana_6 + "GG_catins_gh2ax_gm12878_delta_merged")


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
                              f_data2=ana_5 + 'ALL_h3k4me1_bamvals_merged.csv',
                              outpath=ana_8 + "ALL_print_WT-3h_gh2ax_h3k4me1_merged", res=5000)

""" ############################################################################################ """
""" Obtain features for machine learning model to predict gH2AX span """
# Calculate span width of gH2AX at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.get_span_width(gen, hg19, h2GGhg19, h2WThg19, ana_9 + "GG-gh2ax_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.get_span_width(gen, hg19, h2CThg19, h2WThg19, ana_9 + "CT-gh2ax_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.get_span_width(gen, hg19, h2TAhg19, h2WThg19, ana_9 + "TA-gh2ax_hg19_width")
head = m.load_npheader(ana_9 + "GG-gh2ax_hg19_width.csv")
m.mergerows([m.load_nparray(ana_9 + "GG-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_9 + "CT-gh2ax_hg19_width.csv"),
             m.load_nparray(ana_9 + "TA-gh2ax_hg19_width.csv")],
            ana_9 + "ALL_gh2ax_hg19_width.csv", head)

# Calculate span width of 53BP1 at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
hic.get_span_width(gen, hg19, bpGGhg19, bpWThg19, ana_9 + "GG-53bp1_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
hic.get_span_width(gen, hg19, bpCThg19, bpWThg19, ana_9 + "CT-53bp1_hg19_width")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
hic.get_span_width(gen, hg19, bpTAhg19, bpWThg19, ana_9 + "TA-53bp1_hg19_width")
head = m.load_npheader(ana_9 + "GG-53bp1_hg19_width.csv")
m.mergerows([m.load_nparray(ana_9 + "GG-53bp1_hg19_width.csv"),
             m.load_nparray(ana_9 + "CT-53bp1_hg19_width.csv"),
             m.load_nparray(ana_9 + "TA-53bp1_hg19_width.csv")],
            ana_9 + "ALL_53bp1_hg19_width.csv", head)

# Calculate MRE11 enrichment at each cut site
gen = hic.gen_filter_dist(m.macs_gen(GG_mre11, 1250, hg19, AluGG), 2E6)
m.read_counts(gen, mreGGbam, ana_9 + "GG_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(m.macs_gen(CT_mre11, 1250, hg19, AluCT), 2E6)
m.read_counts(gen, mreCTbam, ana_9 + "CT_mre11_hg19_1250_rc.csv")
gen = hic.gen_filter_dist(m.macs_gen(TA_mre11, 1250, hg19, AluTA), 2E6)
m.read_counts(gen, mreTAbam, ana_9 + "TA_mre11_hg19_1250_rc.csv")
head = m.load_npheader(ana_9 + "GG_mre11_hg19_1250_rc.csv")
m.mergerows([m.load_nparray(ana_9 + "GG_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_9 + "CT_mre11_hg19_1250_rc.csv"),
             m.load_nparray(ana_9 + "TA_mre11_hg19_1250_rc.csv")],
            ana_9 + "ALL_mre11_hg19_1250_rc.csv", head)

# Generate machine learning models for gH2AX
X, y = hic.getXy_insulation(ana_6 + "ALL_iscore_imr90_pos.csv",
                            ana_6 + "ALL_iscore_imr90_neg.csv",
                            ana_9 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_9 + "ALL_gh2ax_hg19_width.csv",
                            ana_9 + "ALL", True)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_9 + "gRF_GG_3h_gh2ax_hg19_1.sav")
ml.ModelTest(X_test, y_test, ana_9 + "gRF_GG_3h_gh2ax_hg19_1.sav")

X, y = hic.getXy_insulation(ana_6 + "ALL_iscore_imr90_pos.csv",
                            ana_6 + "ALL_iscore_imr90_neg.csv",
                            ana_9 + "ALL_mre11_hg19_1250_rc.csv",
                            ana_9 + "ALL_gh2ax_hg19_width.csv",
                            ana_9 + "ALL", False)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.RandomForestTrainGridCV(X_train, y_train, ana_9 + "gRF_GG_3h_gh2ax_hg19_2.sav")
ml.ModelTest(X_test, y_test, ana_9 + "gRF_GG_3h_gh2ax_hg19_2.sav")


epi = [ana_5 + 'ALL_h3k4me1_bamvals_pos.csv',
       ana_5 + 'ALL_h3k4me3_bamvals_pos.csv',
       ana_5 + 'ALL_h3k9me3_bamvals_pos.csv',
       ana_5 + 'ALL_h3k27ac_bamvals_pos.csv',
       ana_5 + 'ALL_h3k36me3_bamvals_pos.csv',
       ana_5 + 'ALL_atacseq_bamvals_pos.csv',
       ana_5 + 'ALL_rnaseq_bamvals_pos.csv',
       ana_5 + 'ALL_ctcf_bamvals_pos.csv',
       ana_5 + 'ALL_smc3_bamvals_pos.csv',
       ana_6 + 'ALL_iscore_gm12878_pos.csv']
lstm.save_Xy_matrix(y_file=ana_4 + 'ALL_gh2ax_WT-3h_hg19_achange_pos.csv', X_files=epi,
                    outfile=ana_9 + 'ALL_Xy_features_pos.csv')

X, y = lstm.load_Xy_matrix(outfile=ana_9 + 'ALL_Xy_features_pos.csv')
X, y = lstm.remove_outliers(X, y, outfile=ana_9 + 'ALL_Xy_features2_pos.csv')
X, y = lstm.modify_matrix(X, y, classifier=False, normalize=False)
X_train, X_test, y_train, y_test = ml.data_split(X, y)
ml.NeuralNetworkTrainGridCV(X_train, y_train, 'hic_nn.sav', classifier=False)
ml.ModelTest(X_test, y_test, 'hic_nn.sav', classifier=False)
