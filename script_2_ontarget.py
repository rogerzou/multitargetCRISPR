"""
Script for:
(1) Determining percentage of predicted ambiguous reads at on-target sites for GG, CT, TA gRNAs.
(2) Determining percentage of ambiguous reads at all target sites for Cas9 and MRE11 ChIP-seq.
(3) Determining enrichment (RPM) at on-target sites for Cas9, MRE11, gH2AX, and 53BP1 ChIP-seq.
(4) Generating peak profiles centered at each cut site for Cas9, MRE11, gH2AX, and 53BP1 ChIP-seq.
(5) Generating profiles of gH2AX and 53BP1 data genome-wide for easy visualization.
"""

import src.chipseq as c
import src.mtss as m
import src.msa as msa
import src.hic as hic
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" File paths """
h2WTin = datadir + "200212_chipseq_WT1/A18_gh2ax_hg38_final.bam"
bpWTin = datadir + "200212_chipseq_WT1/A15_53bp1_hg38_final.bam"
mreGGbam = datadir + "200206_chipseq/AluGG-MRE11_hg38_final.bam"
casGGbam = datadir + "200206_chipseq/AluGG-Cas9_hg38_final.bam"
mreCTbam = datadir + "200316_chipseq/AluCT-mre11-rep1_hg38_final.bam"
casCTbam = datadir + "200316_chipseq/AluCT-cas9-rep1_hg38_final.bam"
mreTAbam = datadir + "200316_chipseq/AluTA-mre11-rep1_hg38_final.bam"
casTAbam = datadir + "200316_chipseq/AluTA-cas9-rep1_hg38_final.bam"
mreGGbam_nD = datadir + "200316_chipseq/AluGG-mre11-noD-rep1_hg38_final.bam"
mreGGbam_PK = datadir + "200316_chipseq/AluGG-mre11-PKi-rep1_hg38_final.bam"
h2GGbam = datadir + "200206_chipseq/AluGG-gH2AX_hg38_final.bam"
bpGGbam = datadir + "200206_chipseq/AluGG-53BP1_hg38_final.bam"
h2TAbam = datadir + "200316_chipseq/AluTA-gh2ax-rep1_hg38_final.bam"
bpTAbam = datadir + "200316_chipseq/AluTA-53bp1-rep1_hg38_final.bam"
h2CTbam = datadir + "200316_chipseq/AluCT-gh2ax-rep1_hg38_final.bam"
bpCTbam = datadir + "200316_chipseq/AluCT-53bp1-rep1_hg38_final.bam"
alnpath_hg38 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_align.csv"

""" macs2 peak detection """
casGGnpk = datadir + "200206_chipseq/macs/AluGG-Cas9_hg38_final_peaks.narrowPeak"
casCTnpk = datadir + "200316_chipseq/macs/AluCT-cas9-rep1_hg38_final_peaks.narrowPeak"
casTAnpk = datadir + "200316_chipseq/macs/AluTA-cas9-rep1_hg38_final_peaks.narrowPeak"

""" Sequences """
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"

""" Set analysis path """
ana = datadir + "Alu_ana_2_ontarget/"
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
""" For Cas9 and MRE11 ChIP-seq at all on-target sites, determine number of reads that
    (1) align once, (2) align multiple times with one optimal alignment, or
    (3) align multiple times with multiple optimal alignments. """
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 750, AluGG),
                         casGGbam, ana_1 + "GG-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 750, AluCT),
                         casCTbam, ana_1 + "CT-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 750, AluTA),
                         casTAbam, ana_1 + "TA-ON_cas9")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG),
                         mreGGbam, ana_1 + "GG-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT),
                         mreCTbam, ana_1 + "CT-ON_mre11")
msa.get_bamfile_pe_reads(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA),
                         mreTAbam, ana_1 + "TA-ON_mre11")
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
""" For Cas9 and MRE11 ChIP-seq at all on- and off-target sites, determine number of reads that
    (1) align once, (2) align multiple times with one optimal alignment, or
    (3) align multiple times with multiple optimal alignments. """
msa.get_bamfile_pe_reads(m.macs_gen(casGGnpk, 750, hg38, AluGG), casGGbam, ana_1 + "GG-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casCTnpk, 750, hg38, AluCT), casCTbam, ana_1 + "CT-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casTAnpk, 750, hg38, AluTA), casTAbam, ana_1 + "TA-C9_cas9")
msa.get_bamfile_pe_reads(m.macs_gen(casGGnpk, 1250, hg38, AluGG), mreGGbam, ana_1 + "GG-C9_mre11")
msa.get_bamfile_pe_reads(m.macs_gen(casCTnpk, 1250, hg38, AluCT), mreCTbam, ana_1 + "CT-C9_mre11")
msa.get_bamfile_pe_reads(m.macs_gen(casTAnpk, 1250, hg38, AluTA), mreTAbam, ana_1 + "TA-C9_mre11")
msa.bowtie2_msa_paired(ana_1 + "GG-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-C9_cas9", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "GG-C9_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "CT-C9_mre11", hg38[1])
msa.bowtie2_msa_paired(ana_1 + "TA-C9_mre11", hg38[1])
msa.parse_msa_sam_paired(ana_1 + "GG-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-C9_cas9_msa")
msa.parse_msa_sam_paired(ana_1 + "GG-C9_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "CT-C9_mre11_msa")
msa.parse_msa_sam_paired(ana_1 + "TA-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "GG-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "CT-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "TA-C9_cas9_msa")
msa.get_msa_stats(ana_1 + "GG-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "CT-C9_mre11_msa")
msa.get_msa_stats(ana_1 + "TA-C9_mre11_msa")


""" ############################################################################################ """
""" For all putative on-target sites, determine paired-end read subsets for Cas9 and MRE11. """
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               mreGGbam, ana_2 + "GG-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT), hg38,
               mreCTbam, ana_2 + "CT-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA), hg38,
               mreTAbam, ana_2 + "TA-ON_mre11_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 750, AluGG), hg38,
               casGGbam, ana_2 + "GG-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 750, AluCT), hg38,
               casCTbam, ana_2 + "CT-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 750, AluTA), hg38,
               casTAbam, ana_2 + "TA-ON_cas9_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               mreGGbam_nD, ana_2 + "GG-ON-noD_mre11_rs")
m.read_subsets(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), hg38,
               mreGGbam_PK, ana_2 + "GG-ON-PKi_mre11_rs")
""" For all putative on-target sites, determine read counts for 53BP1 and gH2AX """
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluGG),
              h2GGbam, ana_2 + "GG-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluCT),
              h2CTbam, ana_2 + "CT-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluTA),
              h2TAbam, ana_2 + "TA-ON_gh2ax_rc.csv")
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluGG),
              bpGGbam, ana_2 + "GG-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluCT),
              bpCTbam, ana_2 + "CT-ON_53bp1_rc.csv")
m.read_counts(msa.target_gen(alnpath_hg38, hg38, 100000, AluTA),
              bpTAbam, ana_2 + "TA-ON_53bp1_rc.csv")
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
    Cas9 and MRE11 ChIP-seq. """
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             mreGGbam, ana_3 + "GG-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluCT),
                             mreCTbam, ana_3 + "CT-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluTA),
                             mreTAbam, ana_3 + "TA-ON_mre11")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluGG),
                             casGGbam, ana_3 + "GG-ON_cas9")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluCT),
                             casCTbam, ana_3 + "CT-ON_cas9")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluTA),
                             casTAbam, ana_3 + "TA-ON_cas9")

""" Generate peak profiles centered at cut site for only abutting reads. """
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluTA),
                             ana_2 + "TA-ON_mre11_rs_abut.bam", ana_3 + "TA-ON_mre11_rs_abut")
m.peak_profile_bp_resolution(msa.target_gen(alnpath_hg38, hg38, 1500, AluTA),
                             ana_2 + "TA-ON_cas9_rs_abut.bam", ana_3 + "TA-ON_cas9_rs_abut")


""" ############################################################################################ """
""" Generate peak profiles centered at the cut site for all putative on-target sites separated by
    2MB from 53BP1 and gH2AX ChIP-seq. """
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, h2GGbam, ana_4 + "GG-ON_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, h2CTbam, ana_4 + "CT-ON_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, h2TAbam, ana_4 + "TA-ON_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "GG-WT_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "CT-WT_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, h2WTin, ana_4 + "TA-WT_gh2ax",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, bpGGbam, ana_4 + "GG-ON_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, bpCTbam, ana_4 + "CT-ON_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, bpTAbam, ana_4 + "TA-ON_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluGG), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "GG-WT_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluCT), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "CT-WT_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)
gen = hic.gen_filter_dist(msa.target_gen(alnpath_hg38, hg38, 1250, AluTA), distance=2000000)
m.peak_profile_wide(gen, hg38, bpWTin, ana_4 + "TA-WT_53bp1",
                    span_rad=2000000, res=10000, wind_rad=5000)


""" ############################################################################################ """
""" Generate profiles of gH2AX and 53BP1 data genome-wide for easy visualization using IGV. """
c.to_wiggle_windows(hg38, h2GGbam, ana_4 + "GG_gh2ax", 500)
c.to_wiggle_windows(hg38, bpGGbam, ana_4 + "GG_53bp1", 500)
c.to_wiggle_windows(hg38, h2CTbam, ana_4 + "CT_gh2ax", 500)
c.to_wiggle_windows(hg38, bpCTbam, ana_4 + "CT_53bp1", 500)
c.to_wiggle_windows(hg38, h2TAbam, ana_4 + "TA_gh2ax", 500)
c.to_wiggle_windows(hg38, bpTAbam, ana_4 + "TA_53bp1", 500)
c.to_wiggle_windows(hg38, h2WTin, ana_4 + "WT_gh2ax", 500)
c.to_wiggle_windows(hg38, bpWTin, ana_4 + "WT_53bp1", 500)
