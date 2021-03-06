"""
Script for analysis of ATAC-seq results
"""

import src.mtss as m
import src.msa as msa
import src.hic as hic
import sys
import os

if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
enc, enc_a = datadir + "public/", datadir + "public/analysis/"
newWTr1 = datadir + "210225_atac/N16_hg38_merged.bam"
newGGr1 = datadir + "210225_atac/N18_hg38_merged.bam"
newWTr2 = datadir + "210225_atac/N19_hg38_merged.bam"
newGGr2 = datadir + "210225_atac/N01_hg38_merged.bam"
newGG00r1 = datadir + "210225_atac/N08_hg38_merged.bam"
newGG10r1 = datadir + "210225_atac/N09_hg38_merged.bam"
newGG30r1 = datadir + "210225_atac/N10_hg38_merged.bam"
newGG00r2 = datadir + "210225_atac/N04_hg38_merged.bam"
newGG10r2 = datadir + "210225_atac/N05_hg38_merged.bam"
newGG30r2 = datadir + "210225_atac/N06_hg38_merged.bam"
mreWTbam = datadir + "200212_chipseq_WT1/A17_mre11_hg38_final.bam"
mreGGbam = datadir + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
alnpath = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = datadir + "Alu_ana_7_atac/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_nucleosomes/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_spanning/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None

""" ############################################################################################ """
""" Obtain read counts in kb scale to mb scale for ATAC-seq """
m.read_counts(msa.target_gen(alnpath, hg38, 250, AluGG),
              newWTr1, ana_1 + "newWTr1_250_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 250, AluGG),
              newGGr1, ana_1 + "newGGr1_250_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 500, AluGG),
              newWTr1, ana_1 + "newWTr1_500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 500, AluGG),
              newGGr1, ana_1 + "newGGr1_500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1000, AluGG),
              newWTr1, ana_1 + "newWTr1_1000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1000, AluGG),
              newGGr1, ana_1 + "newGGr1_1000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              newWTr1, ana_1 + "newWTr1_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              newGGr1, ana_1 + "newGGr1_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 2500, AluGG),
              newWTr1, ana_1 + "newWTr1_2500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 2500, AluGG),
              newGGr1, ana_1 + "newGGr1_2500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              newWTr1, ana_1 + "newWTr1_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              newGGr1, ana_1 + "newGGr1_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 10000, AluGG),
              newWTr1, ana_1 + "newWTr1_10000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 10000, AluGG),
              newGGr1, ana_1 + "newGGr1_10000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 25000, AluGG),
              newWTr1, ana_1 + "newWTr1_25000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 25000, AluGG),
              newGGr1, ana_1 + "newGGr1_25000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 35000, AluGG),
              newWTr1, ana_1 + "newWTr1_35000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 35000, AluGG),
              newGGr1, ana_1 + "newGGr1_35000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 50000, AluGG),
              newWTr1, ana_1 + "newWTr1_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 50000, AluGG),
              newGGr1, ana_1 + "newGGr1_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 75000, AluGG),
              newWTr1, ana_1 + "newWTr1_75000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 75000, AluGG),
              newGGr1, ana_1 + "newGGr1_75000_rc.csv")
atac_3h = [ana_1 + "newGGr1_250_rc.csv",
           ana_1 + "newGGr1_500_rc.csv",
           ana_1 + "newGGr1_1000_rc.csv",
           ana_1 + "newGGr1_1500_rc.csv",
           ana_1 + "newGGr1_2500_rc.csv",
           ana_1 + "newGGr1_5000_rc.csv",
           ana_1 + "newGGr1_10000_rc.csv",
           ana_1 + "newGGr1_25000_rc.csv",
           ana_1 + "newGGr1_35000_rc.csv",
           ana_1 + "newGGr1_50000_rc.csv",
           ana_1 + "newGGr1_75000_rc.csv"]
m.aggregate_values(atac_3h, ana_1 + "newGGr1_merged_rc.csv", col_index=5)
atac_WT = [ana_1 + "newWTr1_250_rc.csv",
           ana_1 + "newWTr1_500_rc.csv",
           ana_1 + "newWTr1_1000_rc.csv",
           ana_1 + "newWTr1_1500_rc.csv",
           ana_1 + "newWTr1_2500_rc.csv",
           ana_1 + "newWTr1_5000_rc.csv",
           ana_1 + "newWTr1_10000_rc.csv",
           ana_1 + "newWTr1_25000_rc.csv",
           ana_1 + "newWTr1_35000_rc.csv",
           ana_1 + "newWTr1_50000_rc.csv",
           ana_1 + "newWTr1_75000_rc.csv"]
m.aggregate_values(atac_WT, ana_1 + "newWTr1_merged_rc.csv", col_index=5)


# """ ############################################################################################ """
# """ Generate peak profiles centered at the cut site for all putative on-target sites """
# # ATAC-seq replicate 1
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newWTr1,
#                     ana_2 + "newWTr1_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGGr1,
#                     ana_2 + "newGGr1_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG00r1,
#                     ana_2 + "newGG00r1_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG10r1,
#                     ana_2 + "newGG10r1_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG30r1,
#                     ana_2 + "newGG30r1_ppw", span_rad=1500, res=1, wind_rad=2)
# # ATAC-seq replicate 2
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newWTr2,
#                     ana_2 + "newWTr2_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGGr2,
#                     ana_2 + "newGGr2_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG00r2,
#                     ana_2 + "newGG00r2_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG10r2,
#                     ana_2 + "newGG10r2_ppw", span_rad=1500, res=1, wind_rad=2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, newGG30r2,
#                     ana_2 + "newGG30r2_ppw", span_rad=1500, res=1, wind_rad=2)
#
#
# """ ############################################################################################ """
# """ Get span widths and read counts for ATAC-seq vs MRE11 ChIP-seq. """
# # ATAC-seq replicate 1
# hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, newGGr1, newWTr1,
#                    ana_2 + "newGGWTr1_width", w_rad=50, skip=5, false_ct=10)
# hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, mreGGbam, mreWTbam,
#                    ana_2 + "mreGGWTr1_width", w_rad=50, skip=5, false_ct=10)
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               newWTr1, ana_2 + "newWTr1_1500_rc.csv")
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               newGGr1, ana_2 + "newGGr1_1500_rc.csv")
# # ATAC-seq replicate 2
# hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, newGGr2, newWTr2,
#                    ana_2 + "newGGWTr2_width", w_rad=50, skip=5, false_ct=10)
# hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, mreGGbam, mreWTbam,
#                    ana_2 + "mreGGWTr2_width", w_rad=50, skip=5, false_ct=10)
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               newWTr2, ana_2 + "newWTr2_1500_rc.csv")
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               newGGr2, ana_2 + "newGGr2_1500_rc.csv")
# # MRE11 ChIP-seq
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               mreWTbam, ana_2 + "mreWTr1_1500_rc.csv")
# m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
#               mreGGbam, ana_2 + "mreGGr1_1500_rc.csv")
#
#
# """ ############################################################################################ """
# """ Determine ATAC-seq nucleosomal vs nucleosome-free fragments """
# # ATAC-seq replicate 1
# m.read_atac_nucleosomes(msa.target_gen(alnpath, hg38, 1500, AluGG), newWTr1, ana_3 + "newWTr1")
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newWTr1_bamnucl.bam",
#                     ana_3 + "newWTr1_nucl", span_rad=1500, res=1, wind_rad=2, norm_type=newWTr1)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newWTr1_bamfree.bam",
#                     ana_3 + "newWTr1_free", span_rad=1500, res=1, wind_rad=2, norm_type=newWTr1)
# m.read_atac_nucleosomes(msa.target_gen(alnpath, hg38, 1500, AluGG), newGGr1, ana_3 + "newGGr1")
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newGGr1_bamnucl.bam",
#                     ana_3 + "newGGr1_nucl", span_rad=1500, res=1, wind_rad=2, norm_type=newGGr1)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newGGr1_bamfree.bam",
#                     ana_3 + "newGGr1_free", span_rad=1500, res=1, wind_rad=2, norm_type=newGGr1)
# # ATAC-seq replicate 2
# m.read_atac_nucleosomes(msa.target_gen(alnpath, hg38, 1500, AluGG), newWTr2, ana_3 + "newWTr2")
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newWTr2_bamnucl.bam",
#                     ana_3 + "newWTr2_nucl", span_rad=1500, res=1, wind_rad=2, norm_type=newWTr2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newWTr2_bamfree.bam",
#                     ana_3 + "newWTr2_free", span_rad=1500, res=1, wind_rad=2, norm_type=newWTr2)
# m.read_atac_nucleosomes(msa.target_gen(alnpath, hg38, 1500, AluGG), newGGr2, ana_3 + "newGGr2")
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newGGr2_bamnucl.bam",
#                     ana_3 + "newGGr2_nucl", span_rad=1500, res=1, wind_rad=2, norm_type=newGGr2)
# m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, ana_3 + "newGGr2_bamfree.bam",
#                     ana_3 + "newGGr2_free", span_rad=1500, res=1, wind_rad=2, norm_type=newGGr2)
#
#
# """ ############################################################################################ """
# """ Determine ATAC-seq reads that span cut site to demonstrate increased accessibility after DNA
#     ligation/repair. """
# # ATAC-seq replicate 1
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newWTr1, ana_4 + "newWTr1_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGGr1, ana_4 + "newGGr1_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG00r1, ana_4 + "newGG00r1_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG10r1, ana_4 + "newGG10r1_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG30r1, ana_4 + "newGG30r1_rs")
# # ATAC-seq replicate 2
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newWTr2, ana_4 + "newWTr2_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGGr2, ana_4 + "newGGr2_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG00r2, ana_4 + "newGG00r2_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG10r2, ana_4 + "newGG10r2_rs")
# m.read_subsets(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38,
#                newGG30r2, ana_4 + "newGG30r2_rs")
#
# # ATAC-seq replicate 1
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newWTr1_rs_span.bam", ana_4 + "newWTr1_span",
#                              norm_type=newWTr1)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGGr1_rs_span.bam", ana_4 + "newGGr1_span",
#                              norm_type=newGGr1)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG00r1_rs_span.bam", ana_4 + "newGG00r1_span",
#                              norm_type=newGG00r1)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG10r1_rs_span.bam", ana_4 + "newGG10r1_span",
#                              norm_type=newGG10r1)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG30r1_rs_span.bam", ana_4 + "newGG30r1_span",
#                              norm_type=newGG30r1)
# # ATAC-seq replicate 2
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newWTr2_rs_span.bam", ana_4 + "newWTr2_span",
#                              norm_type=newWTr2)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGGr2_rs_span.bam", ana_4 + "newGGr2_span",
#                              norm_type=newGGr2)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG00r2_rs_span.bam", ana_4 + "newGG00r2_span",
#                              norm_type=newGG00r2)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG10r2_rs_span.bam", ana_4 + "newGG10r2_span",
#                              norm_type=newGG10r2)
# m.peak_profile_bp_resolution(msa.target_gen(alnpath, hg38, 1500, AluGG),
#                              ana_4 + "newGG30r2_rs_span.bam", ana_4 + "newGG30r2_span",
#                              norm_type=newGG30r2)
