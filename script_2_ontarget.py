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
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" File paths """
h2WTin = labhome + "200212_chipseq_WT1/A18_gh2ax_hg38_final.bam"
bpWTin = labhome + "200212_chipseq_WT1/A15_53bp1_hg38_final.bam"
mreGGbam = labhome + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
casGGbam = labhome + "200206_chipseq/AluGG-Cas9_hg38_final.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_hg38_final.bam"
casCTbam = labhome + "200316_chipseq/AluCT-cas9-rep1_hg38_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_hg38_final.bam"
casTAbam = labhome + "200316_chipseq/AluTA-cas9-rep1_hg38_final.bam"
mreGGbam_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_hg38_final.bam"
mreGGbam_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_hg38_final.bam"
h2GGbam = labhome + "200206_chipseq/AluGG-gH2AX_hg38_final.bam"
bpGGbam = labhome + "200206_chipseq/AluGG-53BP1_hg38_final.bam"
h2TAbam = labhome + "200316_chipseq/AluTA-gh2ax-rep1_hg38_final.bam"
bpTAbam = labhome + "200316_chipseq/AluTA-53bp1-rep1_hg38_final.bam"
h2CTbam = labhome + "200316_chipseq/AluCT-gh2ax-rep1_hg38_final.bam"
bpCTbam = labhome + "200316_chipseq/AluCT-53bp1-rep1_hg38_final.bam"
alnpath = labhome + "Alu_ana_1_putative/1_protosearch/psearch_align.csv"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = labhome + "Alu_ana_2_ontarget/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_msa_ontarget/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_subsets_counts/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_profiles_mre11-cas9/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_profiles_53bp1-gh2ax/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


""" ############################################################################################ """
""" Determine multiple sequence alignments (up to 300) for each ChIP-seq paired-end reads around
    putative on-target sites (Figures 1A-D, S1) """
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluGG), casGGbam, ana_1 + "GG-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluCT), casCTbam, ana_1 + "CT-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 750, AluTA), casTAbam, ana_1 + "TA-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluGG), mreGGbam, ana_1 + "GG-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluCT), mreCTbam, ana_1 + "CT-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath, 1250, AluTA), mreTAbam, ana_1 + "TA-ON_mre11")
msa.bowtie2_msa_paired(ana_1 + "GG-ON_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-ON_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-ON_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "GG-ON_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-ON_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-ON_mre11", hg38[1])
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
""" For all putative on-target sites, determine paired-end read subsets for Cas9 and MRE11 """
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam, ana_2 + "GG-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluCT), mreCTbam, ana_2 + "CT-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluTA), mreTAbam, ana_2 + "TA-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluGG), casGGbam, ana_2 + "GG-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluCT), casCTbam, ana_2 + "CT-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 750, AluTA), casTAbam, ana_2 + "TA-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam_nD, ana_2 + "GG-ON-noD_mre11_rs")
m.read_subsets(msa.target_gen(alnpath, 1250, AluGG), mreGGbam_PK, ana_2 + "GG-ON-PKi_mre11_rs")
""" For all putative on-target sites, determine read counts for 53BP1 and gH2AX """
m.read_counts(msa.target_gen(alnpath, 100000, AluGG), h2GGbam, ana_2 + "GG-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluCT), h2CTbam, ana_2 + "CT-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluTA), h2TAbam, ana_2 + "TA-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluGG), bpGGbam, ana_2 + "GG-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluCT), bpCTbam, ana_2 + "CT-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath, 100000, AluTA), bpTAbam, ana_2 + "TA-ON_53bp1_rc.csv")
""" Merge datasets """
num_cols = [1, 1]
head = m.load_npheader(ana_2 + "GG-ON_mre11_rs.csv") + ", 53bp1, gh2ax"
GG_mre = m.load_nparray(ana_2 + "GG-ON_mre11_rs.csv")
CT_mre = m.load_nparray(ana_2 + "CT-ON_mre11_rs.csv")
TA_mre = m.load_nparray(ana_2 + "TA-ON_mre11_rs.csv")
GG_cas = m.load_nparray(ana_2 + "GG-ON_cas9_rs.csv")
CT_cas = m.load_nparray(ana_2 + "CT-ON_cas9_rs.csv")
TA_cas = m.load_nparray(ana_2 + "TA-ON_cas9_rs.csv")
GG = [m.load_nparray(ana_2 + "GG-ON_53bp1_rc.csv"), m.load_nparray(ana_2 + "GG-ON_gh2ax_rc.csv")]
CT = [m.load_nparray(ana_2 + "CT-ON_53bp1_rc.csv"), m.load_nparray(ana_2 + "CT-ON_gh2ax_rc.csv")]
TA = [m.load_nparray(ana_2 + "TA-ON_53bp1_rc.csv"), m.load_nparray(ana_2 + "TA-ON_gh2ax_rc.csv")]
GG_mre_m = m.mergesubsetcounts(GG_mre, GG, num_cols, ana_2 + "GG-ON_mre_merged.csv", head)
CT_mre_m = m.mergesubsetcounts(CT_mre, CT, num_cols, ana_2 + "CT-ON_mre_merged.csv", head)
TA_mre_m = m.mergesubsetcounts(TA_mre, TA, num_cols, ana_2 + "TA-ON_mre_merged.csv", head)
GG_cas_m = m.mergesubsetcounts(GG_cas, GG, num_cols, ana_2 + "GG-ON_cas_merged.csv", head)
CT_cas_m = m.mergesubsetcounts(CT_cas, CT, num_cols, ana_2 + "CT-ON_cas_merged.csv", head)
TA_cas_m = m.mergesubsetcounts(TA_cas, TA, num_cols, ana_2 + "TA-ON_cas_merged.csv", head)
m.mergerows([GG_mre_m, CT_mre_m, TA_mre_m], ana_2 + "ALL-ON_mre_merged.csv", head)
m.mergerows([GG_cas_m, CT_cas_m, TA_cas_m], ana_2 + "ALL-ON_cas_merged.csv", head)


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites from
    Cas9 and MRE11 ChIP-seq. (Figures 1E-F) """
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
""" Generate peak profiles centered at the cut site for all putative on-target sites separated by
    2MB from53BP1 and gH2AX ChIP-seq. (Figures 1E-F) """
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, h2GGbam, ana_4 + "GG-ON_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, h2CTbam, ana_4 + "CT-ON_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, h2TAbam, ana_4 + "TA-ON_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, bpGGbam, ana_4 + "GG-ON_53bp1", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, bpCTbam, ana_4 + "CT-ON_53bp1", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, bpTAbam, ana_4 + "TA-ON_53bp1", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "GG-WTneg_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "CT-WTneg_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "TA-WTneg_gh2ax", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "GG-WTneg_53bp1", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "CT-WTneg_53bp1", span_rad=1000000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "TA-WTneg_53bp1", span_rad=1000000)
