"""
Script for:
(1) Determination of on-target paired-end read uniqueness for Cas9 and MRE11 ChIP-seq.
(2) Determination of enrichment at on-target sites for Cas9, MRE11, gH2AX, and 53BP1 ChIP-seq.
(3) Generate peak profiles centered at each cut site for  Cas9, MRE11, gH2AX, and 53BP1 ChIP-seq.
"""

import src.mtss as m
import src.msa as msa
import src.hic as hic
import sys
import os

""" Determine run paths based on operating system """
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
mreGGbam = labhome + "200206_chipseq/16_AluGG-MRE11_final.bam"
casGGbam = labhome + "200206_chipseq/15_AluGG-Cas9_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_rmdup.bam"
casTAbam = labhome + "200316_chipseq/AluTA-cas9-rep1_rmdup.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_rmdup.bam"
casCTbam = labhome + "200316_chipseq/AluCT-cas9-rep1_rmdup.bam"
mreGGbam_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_rmdup.bam"
mreGGbam_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_rmdup.bam"
h2GGbam = labhome + "200206_chipseq/13_AluGG-gH2AX_final.bam"
bpGGbam = labhome + "200206_chipseq/14_AluGG-53BP1_final.bam"
h2TAbam = labhome + "200316_chipseq/AluTA-gh2ax-rep1_rmdup.bam"
bpTAbam = labhome + "200316_chipseq/AluTA-53bp1-rep1_rmdup.bam"
h2CTbam = labhome + "200316_chipseq/AluCT-gh2ax-rep1_rmdup.bam"
bpCTbam = labhome + "200316_chipseq/AluCT-53bp1-rep1_rmdup.bam"
alnpath = labhome + "Alu_ana_1_putative/protosearch/psearch_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = labhome + "Alu_ana_2_ontarget/"
os.makedirs(ana) if not os.path.exists(ana) else None


""" ############################################################################################ """
""" Determine multiple sequence alignments (up to 300) for each ChIP-seq paired-end reads around
    putative on-target sites (Figure 1G) """
ana_1 = ana + "1_msa_ontarget/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluGG), casGGbam, ana_1 + "GG-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluCT), casCTbam, ana_1 + "CT-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluTA), casTAbam, ana_1 + "TA-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluGG), mreGGbam, ana_1 + "GG-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluCT), mreCTbam, ana_1 + "CT-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluTA), mreTAbam, ana_1 + "TA-ON_mre11")
msa.bowtie2_msa_paired(ana_1 + "GG-ON_cas9", hg38)
msa.bowtie2_msa_paired(ana_1 + "CT-ON_cas9", hg38)
msa.bowtie2_msa_paired(ana_1 + "TA-ON_cas9", hg38)
msa.bowtie2_msa_paired(ana_1 + "GG-ON_mre11", hg38)
msa.bowtie2_msa_paired(ana_1 + "CT-ON_mre11", hg38)
msa.bowtie2_msa_paired(ana_1 + "TA-ON_mre11", hg38)
msa.parse_msa_sam_paired(ana_1 + "GG-ON_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-ON_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-ON_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "GG-ON_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-ON_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-ON_mre11_msa")
msa.get_msa_stats(ana_1 + "GG-ON_cas9_msa")
msa.get_msa_stats(ana_1 + "CT-ON_cas9_msa")
msa.get_msa_stats(ana_1 + "TA-ON_cas9_msa")
msa.get_msa_stats(ana_1 + "GG-ON_mre11_msa")
msa.get_msa_stats(ana_1 + "CT-ON_mre11_msa")
msa.get_msa_stats(ana_1 + "TA-ON_mre11_msa")


""" ############################################################################################ """
""" Determine paired-end read subsets for all putative on-target sites (Figure 1) """
ana_2 = ana + "2_subsets_counts/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam, ana_2 + "GG-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluCT), mreCTbam, ana_2 + "CT-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluTA), mreTAbam, ana_2 + "TA-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluGG), casGGbam, ana_2 + "GG-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluCT), casCTbam, ana_2 + "CT-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluTA), casTAbam, ana_2 + "TA-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam_nD, ana_2 + "GG-ON-noD_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam_PK, ana_2 + "GG-ON-PKi_mre11_rs")
GG_mre = m.load_nparray(ana_2 + "GG-ON_mre11_rs.csv")
CT_mre = m.load_nparray(ana_2 + "CT-ON_mre11_rs.csv")
TA_mre = m.load_nparray(ana_2 + "TA-ON_mre11_rs.csv")
GG_cas = m.load_nparray(ana_2 + "GG-ON_cas9_rs.csv")
CT_cas = m.load_nparray(ana_2 + "CT-ON_cas9_rs.csv")
TA_cas = m.load_nparray(ana_2 + "TA-ON_cas9_rs.csv")
m.mergerows([GG_mre, CT_mre, TA_mre], ana_2 + "MERGE-ON_mre_rs.csv")
m.mergerows([GG_cas, CT_cas, TA_cas], ana_2 + "MERGE-ON_cas_rs.csv")

""" Get corresponding enrichment data from gH2AX and 53BP1 (Figure 1) """
m.read_counts(msa.target_gen(alnpath, 100000, AluGG), h2GGbam, ana_2 + "GG-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluGG), bpGGbam, ana_2 + "GG-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluCT), h2CTbam, ana_2 + "CT-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluCT), bpCTbam, ana_2 + "CT-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluTA), h2TAbam, ana_2 + "TA-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluTA), bpTAbam, ana_2 + "TA-ON_53bp1_rc.csv")
GG_h2 = m.load_nparray(ana_2 + "GG-ON_gh2ax_rc.csv")
CT_h2 = m.load_nparray(ana_2 + "CT-ON_gh2ax_rc.csv")
TA_h2 = m.load_nparray(ana_2 + "TA-ON_gh2ax_rc.csv")
GG_bp = m.load_nparray(ana_2 + "GG-ON_53bp1_rc.csv")
CT_bp = m.load_nparray(ana_2 + "CT-ON_53bp1_rc.csv")
TA_bp = m.load_nparray(ana_2 + "TA-ON_53bp1_rc.csv")
m.mergerows([GG_h2, CT_h2, TA_h2], ana_2 + "MERGE-ON_gh2ax_rc.csv")
m.mergerows([GG_bp, CT_bp, TA_bp], ana_2 + "MERGE-ON_53bp1_rc.csv")


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq.
    (Figures 1D to 1F) """
ana_3 = ana + "3_profiles_mre11-cas9/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluGG), mreGGbam, ana_3 + "GG-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluCT), mreCTbam, ana_3 + "CT-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluTA), mreTAbam, ana_3 + "TA-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluGG), casGGbam, ana_3 + "GG-ON_cas9")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluCT), casCTbam, ana_3 + "CT-ON_cas9")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluTA), casTAbam, ana_3 + "TA-ON_cas9")

""" Generate peak profiles centered at cut site for only abutting reads (Figure S4) """
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluTA),
                             ana_2 + "TA-ON_mre11_rs_abut.bam", ana_3 + "TA-ON_mre11_rs_abut")
m.peak_profile_bp_resolution(msa.target_gen(alnpath, 1500, AluTA),
                             ana_2 + "TA-ON_cas9_rs_abut.bam", ana_3 + "TA-ON_cas9_rs_abut")


""" ############################################################################################ """
""" Get 53BP1 and gH2AX span and profiles """
ana_4 = ana + "4_profiles_gh2ax_53bp1/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
hic.get_span_width(gen, h2GGbam, h2WTin, ana_4 + "GG-ON_gh2ax_span")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
hic.get_span_width(gen, h2CTbam, h2WTin, ana_4 + "CT-ON_gh2ax_span")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
hic.get_span_width(gen, h2TAbam, h2WTin, ana_4 + "TA-ON_gh2ax_span")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
hic.get_span_width(gen, bpGGbam, bpWTin, ana_4 + "GG-ON_53bp1_span")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
hic.get_span_width(gen, bpCTbam, bpWTin, ana_4 + "CT-ON_53bp1_span")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
hic.get_span_width(gen, bpTAbam, bpWTin, ana_4 + "TA-ON_53bp1_span")

""" Generate peak profiles centered at the cut site for all putative on-target sites from
    53BP1 and gH2AX ChIP-seq. (Figures ) """
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, h2GGbam, ana_4 + "GG-ON_gh2ax")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, h2CTbam, ana_4 + "CT-ON_gh2ax")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, h2TAbam, ana_4 + "TA-ON_gh2ax")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, bpGGbam, ana_4 + "GG-ON_53bp1")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, bpCTbam, ana_4 + "CT-ON_53bp1")
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, bpTAbam, ana_4 + "TA-ON_53bp1")
