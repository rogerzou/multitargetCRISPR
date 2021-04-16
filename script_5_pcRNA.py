"""
Script for:
(1) Determination of time-resolved MRE11 enrichment after activation of cgRNA.
"""

import src.mtss as m
import src.msa as msa
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
atacWTr1 = datadir + "210225_atac/N16_hg38_merged.bam"
atacWTr2 = datadir + "210225_atac/N19_hg38_merged.bam"
mreWTbam = datadir + "200212_chipseq_WT1/A17_mre11_hg38_final.bam"
mreGG0hPCr1 = datadir + "210325_chipseq/A01_hg38_merged.bam"
mreGG1hPCr1 = datadir + "210325_chipseq/A02_hg38_merged.bam"
mreGG2hPCr1 = datadir + "210325_chipseq/A03_hg38_merged.bam"
mreGG4hPCr1 = datadir + "210325_chipseq/A04_hg38_merged.bam"
mreGG0hPCr2 = datadir + "210325_chipseq/A05_hg38_merged.bam"
mreGG1hPCr2 = datadir + "210325_chipseq/A06_hg38_merged.bam"
mreGG2hPCr2 = datadir + "210325_chipseq/A08_hg38_merged.bam"
mreGG4hPCr2 = datadir + "210325_chipseq/A09_hg38_merged.bam"
atacGG0hPCr1 = datadir + "210325_atac/N702_hg38_merged.bam"
atacGG1hPCr1 = datadir + "210325_atac/N703_hg38_merged.bam"
atacGG2hPCr1 = datadir + "210325_atac/N704_hg38_merged.bam"
atacGG4hPCr1 = datadir + "210325_atac/N705_hg38_merged.bam"
atacGG0hPCr2 = datadir + "210325_atac/N706_hg38_merged.bam"
atacGG1hPCr2 = datadir + "210325_atac/N707_hg38_merged.bam"
atacGG2hPCr2 = datadir + "210325_atac/N708_hg38_merged.bam"
atacGG4hPCr2 = datadir + "210325_atac/N709_hg38_merged.bam"
alnpath_hg38 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"

""" macs2 output """
casGG3h_npk1 = datadir + "200804_chipseq/macs/A03_hg38_final_peaks.narrowPeak"

""" Set analysis path """
ana = datadir + "Alu_ana_5_pcRNA/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_nucleosomes/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_peak_profiles/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Determine ATAC-seq nucleosomal vs nucleosome-free fragments """
# ATAC-seq cgRNA
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG0hPCr1, ana_1 + "atacGG0hPCr1")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG1hPCr1, ana_1 + "atacGG1hPCr1")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG2hPCr1, ana_1 + "atacGG2hPCr1")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG4hPCr1, ana_1 + "atacGG4hPCr1")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG0hPCr2, ana_1 + "atacGG0hPCr2")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG1hPCr2, ana_1 + "atacGG1hPCr2")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG2hPCr2, ana_1 + "atacGG2hPCr2")
m.read_atac_nucleosomes(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                        atacGG4hPCr2, ana_1 + "atacGG4hPCr2")


""" ############################################################################################ """
"""  """
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacWTr1,
                    ana_2 + "atacWTr1_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG0hPCr1,
                    ana_2 + "atacGG0hPCr1_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG1hPCr1,
                    ana_2 + "atacGG1hPCr1_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG2hPCr1,
                    ana_2 + "atacGG2hPCr1_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG4hPCr1,
                    ana_2 + "atacGG4hPCr1_ppw", span_rad=1500, res=1, wind_rad=2)

m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacWTr2,
                    ana_2 + "atacWTr2_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG0hPCr2,
                    ana_2 + "atacGG0hPCr2_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG1hPCr2,
                    ana_2 + "atacGG1hPCr2_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG2hPCr2,
                    ana_2 + "atacGG2hPCr2_ppw", span_rad=1500, res=1, wind_rad=2)
m.peak_profile_wide(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG), hg38, atacGG4hPCr2,
                    ana_2 + "atacGG4hPCr2_ppw", span_rad=1500, res=1, wind_rad=2)


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq. """
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreWTbam, ana_2 + "mreWTbam")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG0hPCr1, ana_2 + "mreGG0hPCr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG1hPCr1, ana_2 + "mreGG1hPCr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG2hPCr1, ana_2 + "mreGG2hPCr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG4hPCr1, ana_2 + "mreGG4hPCr1")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG0hPCr2, ana_2 + "mreGG0hPCr2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG1hPCr2, ana_2 + "mreGG1hPCr2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG2hPCr2, ana_2 + "mreGG2hPCr2")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGG4hPCr2, ana_2 + "mreGG4hPCr2")
