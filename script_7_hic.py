"""
Script for:
(1) Analysis of H2AX and 53BP1 data

"""

import src.chipseq as c
import src.mtss as m
import src.hic as hic
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
    
""" macs2 peak detection """
mreGGnpk = labhome + "200206_chipseq/macs/AluGG-MRE11_hg19_final_peaks.narrowPeak"
mreCTnpk = labhome + "200316_chipseq/macs/AluCT-mre11-rep1_hg19_final_peaks.narrowPeak"
mreTAnpk = labhome + "200316_chipseq/macs/AluTA-mre11-rep1_hg19_final_peaks.narrowPeak"

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
ana_4 = ana + "4_insu_scores/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None


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
""" Obtain percent and absolute changes of gH2AX and 53BP1 wig files """
c.percentchange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "GG_53bp1_hg19.wig",
                ana_2 + "GG_53bp1_00m-3h_hg19_pchange")
c.percentchange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "GG_gh2ax_hg19.wig",
                ana_2 + "GG_gh2ax_00m-3h_hg19_pchange")
c.percentchange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "CT_53bp1_hg19.wig",
                ana_2 + "CT_53bp1_00m-3h_hg19_pchange")
c.percentchange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "CT_gh2ax_hg19.wig",
                ana_2 + "CT_gh2ax_00m-3h_hg19_pchange")
c.percentchange(ana_1 + "WT_53bp1_hg19.wig", ana_1 + "TA_53bp1_hg19.wig",
                ana_2 + "TA_53bp1_00m-3h_hg19_pchange")
c.percentchange(ana_1 + "WT_gh2ax_hg19.wig", ana_1 + "TA_gh2ax_hg19.wig",
                ana_2 + "TA_gh2ax_00m-3h_hg19_pchange")

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

c.percentchange(ana_1 + "GG_cgH2_00m_hg19.wig", ana_1 + "GG_cgH2_10m_hg19.wig",
                ana_2 + "GG_cgH2_00m-10m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgH2_00m_hg19.wig", ana_1 + "GG_cgH2_30m_hg19.wig",
                ana_2 + "GG_cgH2_00m-30m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgH2_10m_hg19.wig", ana_1 + "GG_cgH2_30m_hg19.wig",
                ana_2 + "GG_cgH2_10m-30m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgBP_00m_hg19.wig", ana_1 + "GG_cgBP_10m_hg19.wig",
                ana_2 + "GG_cgBP_00m-10m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgBP_00m_hg19.wig", ana_1 + "GG_cgBP_30m_hg19.wig",
                ana_2 + "GG_cgBP_00m-30m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgBP_10m_hg19.wig", ana_1 + "GG_cgBP_30m_hg19.wig",
                ana_2 + "GG_cgBP_10m-30m_hg19_pchange")
c.percentchange(ana_1 + "GG_cgBP_30m_hg19.wig", ana_1 + "GG_53bp1_hg19.wig",
                ana_2 + "GG_cgBP_30m-3h_hg19_pchange")
c.percentchange(ana_1 + "GG_cgH2_30m_hg19.wig", ana_1 + "GG_gh2ax_hg19.wig",
                ana_2 + "GG_cgH2_30m-3h_hg19_pchange")

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

# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, h2GGhg19, ana_1 + "GG_gh2ax_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, bpGGhg19, ana_1 + "GG_53bp1_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
# m.peak_profile_wide(gen, hg19, h2CThg19, ana_1 + "CT_gh2ax_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
# m.peak_profile_wide(gen, hg19, bpCThg19, ana_1 + "CT_53bp1_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
# m.peak_profile_wide(gen, hg19, h2TAhg19, ana_1 + "TA_gh2ax_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
# m.peak_profile_wide(gen, hg19, bpTAhg19, ana_1 + "TA_53bp1_hg19", res=5000, wind_rad=20000)

# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, h2GG00m_cg, ana_2 + "GG_cgH2_00m_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, h2GG10m_cg, ana_2 + "GG_cgH2_10m_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, h2GG30m_cg, ana_2 + "GG_cgH2_30m_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, bpGG00m_cg, ana_2 + "GG_cgBP_00m_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, bpGG10m_cg, ana_2 + "GG_cgBP_10m_hg19", res=5000, wind_rad=20000)
# gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
# m.peak_profile_wide(gen, hg19, bpGG30m_cg, ana_2 + "GG_cgBP_30m_hg19", res=5000, wind_rad=20000)


""" ############################################################################################ """
""" Obtain 4C-seq profiles from Hi-C data """
gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_GM12878", labhome + "public_HiC/GM12878", 5, 2E6)

gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_HUVEC", labhome + "public_HiC/HUVEC", 5, 2E6)

gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_IMR90", labhome + "public_HiC/IMR90", 5, 2E6)

gen = hic.gen_filter_dist(m.macs_gen(mreGGnpk, 1250, hg19, AluGG), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "GG_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreCTnpk, 1250, hg19, AluCT), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "CT_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)
gen = hic.gen_filter_dist(m.macs_gen(mreTAnpk, 1250, hg19, AluTA), 2E6)
hic.rao_fourCseq_gen(gen, ana_3 + "TA_4Cseq_hg19_NHEK", labhome + "public_HiC/NHEK", 5, 2E6)

""" ############################################################################################ """
""" Obtain insulation scores from raw Hi-C matrices """
wig_filename = "all_chr_5kb_GM12878" # modify as needed
hic.gen_insu_scores(raw_matrix_path, labhome, ana_4, wig_filename, 250000, 5000)

    

