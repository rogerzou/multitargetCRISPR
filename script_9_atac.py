"""
Script for analysis of ATAC-seq results
"""

import src.mtss as m
import src.msa as msa
import src.hic as hic
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
atacWTse = labhome + "201110_atac/N03_sorted.bam"
atacGGse3h = labhome + "201110_atac/N02_sorted.bam"
atacWTpe = labhome + "201207_atac/N701_hg38_final.bam"
atacGGpe3h = labhome + "201207_atac/N703_hg38_final.bam"
atacGGpe00m = labhome + "201207_atac/N706_hg38_final.bam"
atacGGpe10m = labhome + "201207_atac/N704_hg38_final.bam"
atacGGpe30m = labhome + "201207_atac/N705_hg38_final.bam"

mreWTbam = labhome + "200212_chipseq_WT1/A17_mre11_hg38_final.bam"
mreGGbam = labhome + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
alnpath = labhome + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"
h3k4me1_1 = enc + "hg38/H3K4me1_HEK293_rep1_ENCFF909ESY.bam"     # enhancers
h3k4me3_1 = enc + "hg38/H3K4me3_HEK293_rep1_ENCFF912BYL.bam"     # transcription activation
h3k9me3_1 = enc + "hg38/H3K9me3_HEK293_rep1_ENCFF141ZEQ.bam"     # heterochromatin
h3k27ac_1 = enc + "hg38/H3K27ac_HEK293_rep1_ENCFF588KSR.bam"     # enhancer
h3k36me3_1 = enc + "hg38/H3K36me3_HEK293_rep1_ENCFF593SUW.bam"   # gene bodies
dnasei_1 = enc + "hg38/DNaseI_HEK293T_ENCFF120XFB.bam"           # DNase I hypersensitivity
atac_1 = enc + "hg38/ATACseq_HEK293_SRR6418075.bam"              # ATAC-seq (medium deep)
mnase_1 = enc + "hg38/MNaseseq_HEK293_ERR2403161.bam"            # MNase-seq
rna_3 = enc + "hg38/RNAseq_HEK293_SRR5627161.bam"                # RNA-seq #3

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
GG3h_npk1 = labhome + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"
GG_mre11 = labhome + "200206_chipseq/macs/AluGG-MRE11_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = labhome + "Alu_ana_9_atac/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_subsets/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Obtain read counts in kb scale to mb scale for ATAC-seq """
m.read_counts(msa.target_gen(alnpath, hg38, 250, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_250_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 250, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_250_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 500, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 500, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_1000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_1000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 2500, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_2500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 2500, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_2500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 5000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_5000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 10000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_10000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 10000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_10000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 25000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_25000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 25000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_25000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 35000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_35000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 35000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_35000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 50000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 50000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_50000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 75000, AluGG),
              atacWTse, ana_1 + "WT-ON_atac_75000_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 75000, AluGG),
              atacGGse3h, ana_1 + "GG-ON_atac_3h_75000_rc.csv")


atac_3h = [ana_1 + "GG-ON_atac_3h_250_rc.csv",
           ana_1 + "GG-ON_atac_3h_500_rc.csv",
           ana_1 + "GG-ON_atac_3h_1000_rc.csv",
           ana_1 + "GG-ON_atac_3h_1500_rc.csv",
           ana_1 + "GG-ON_atac_3h_2500_rc.csv",
           ana_1 + "GG-ON_atac_3h_5000_rc.csv",
           ana_1 + "GG-ON_atac_3h_10000_rc.csv",
           ana_1 + "GG-ON_atac_3h_25000_rc.csv",
           ana_1 + "GG-ON_atac_3h_35000_rc.csv",
           ana_1 + "GG-ON_atac_3h_50000_rc.csv",
           ana_1 + "GG-ON_atac_3h_75000_rc.csv"]
m.aggregate_values(atac_3h, ana_1 + "GG-ON_atac_3h_merged_rc.csv", col_index=5)
atac_WT = [ana_1 + "WT-ON_atac_250_rc.csv",
           ana_1 + "WT-ON_atac_500_rc.csv",
           ana_1 + "WT-ON_atac_1000_rc.csv",
           ana_1 + "WT-ON_atac_1500_rc.csv",
           ana_1 + "WT-ON_atac_2500_rc.csv",
           ana_1 + "WT-ON_atac_5000_rc.csv",
           ana_1 + "WT-ON_atac_10000_rc.csv",
           ana_1 + "WT-ON_atac_25000_rc.csv",
           ana_1 + "WT-ON_atac_35000_rc.csv",
           ana_1 + "WT-ON_atac_50000_rc.csv",
           ana_1 + "WT-ON_atac_75000_rc.csv"]
m.aggregate_values(atac_WT, ana_1 + "WT-ON_atac_merged_rc.csv", col_index=5)


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites """
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacWTse,
                    ana_2 + "WT_atac", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacGGse3h,
                    ana_2 + "GG-3h_atac", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacWTpe,
                    ana_2 + "atac_pe_WT_1w", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacGGpe3h,
                    ana_2 + "atac_pe_GG_3h_1w", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacGGpe00m,
                    ana_2 + "atac_pe_cgGG_00m_1w", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacGGpe10m,
                    ana_2 + "atac_pe_cgGG_10m_1w", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath, hg38, 1500, AluGG), hg38, atacGGpe30m,
                    ana_2 + "atac_pe_cgGG_30m_1w", span_rad=1500, res=1, wind_rad=2)

hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, atacGGse3h, atacWTse,
                   ana_2 + "GG-atac_hg38_width", w_rad=50, skip=5, false_ct=10)
hic.get_span_width(msa.target_gen(alnpath, hg38, 100, AluGG), hg38, mreGGbam, mreWTbam,
                   ana_2 + "GG-mre11_hg38_width", w_rad=50, skip=5, false_ct=10)
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacWTse, ana_2 + "WT-ON_atac_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG),
              atacGGse3h, ana_2 + "GG-ON_atac_3h_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG), mreWTbam,
              ana_2 + "WT-mre11_hg38_1500_rc.csv")
m.read_counts(msa.target_gen(alnpath, hg38, 1500, AluGG), mreGGbam,
              ana_2 + "GG-mre11_hg38_1500_rc.csv")
