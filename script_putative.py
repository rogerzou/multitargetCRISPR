"""
Script for:
(1) Computational identification of potential gRNAs from Alu repetitive sequence.
(2) Determination of putative genome-wide on-target sites for each gRNA.
(3) Evaluate uniqueness of each paired-end ChIP-seq read.
"""

import src.mtss as m


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana = labhome + "Alu_analysis1/"
mreGGbam = labhome + "200206_chipseq/16_AluGG-MRE11_final.bam"
casGGbam = labhome + "200206_chipseq/15_AluGG-Cas9_final.bam"
mreTAbam = labhome + "200316_chipseq/AluTA-mre11-rep1_rmdup.bam"
casTAbam = labhome + "200316_chipseq/AluTA-cas9-rep1_rmdup.bam"
mreCTbam = labhome + "200316_chipseq/AluCT-mre11-rep1_rmdup.bam"
casCTbam = labhome + "200316_chipseq/AluCT-cas9-rep1_rmdup.bam"
mreGGbam_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_rmdup.bam"
mreGGbam_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_rmdup.bam"
alnpath = ana + 'basemut_align.csv'

""" Sequences """
Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTG" \
      "GCCAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC" \
      "TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG" \
      "ACTCCGTCTC"
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"


""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, calculate distances between adjacent targets
    (Figures 1A to 1C) """
basemut = ana + 'basemut'
m.get_targets_fasta(basemut, Alu, 9)   # get list of all protospacer sequences as FASTA file
m.get_targets_bowtie2(basemut, hg38)   # from FASTA, MSA up to 1000 locations in hg38 as SAM file
m.get_targets_stats(basemut)           # from SAM, summarize MSA as well as location in genes
m.get_targets_dist(alnpath, basemut)   # from MSA, get distance between each putative target site

""" Determine paired-end read subsets for all putative on-target sites (Figure 1) """
m.read_subsets(m.target_gen(alnpath, 1250, AluGG), mreGGbam, ana + "GG-ON_mre11_rs")
m.read_subsets(m.target_gen(alnpath, 1250, AluCT), mreCTbam, ana + "CT-ON_mre11_rs")
m.read_subsets(m.target_gen(alnpath, 1250, AluTA), mreTAbam, ana + "TA-ON_mre11_rs")
m.read_subsets(m.target_gen(alnpath, 750, AluGG), casGGbam, ana + "GG-ON_cas9_rs")
m.read_subsets(m.target_gen(alnpath, 750, AluCT), casCTbam, ana + "CT-ON_cas9_rs")
m.read_subsets(m.target_gen(alnpath, 750, AluTA), casTAbam, ana + "TA-ON_cas9_rs")
m.read_subsets(m.target_gen(alnpath, 1250, AluGG), mreGGbam_nD, ana + "GG-ON-noD_mre11_rs")
m.read_subsets(m.target_gen(alnpath, 1250, AluGG), mreGGbam_PK, ana + "GG-ON-PKi_mre11_rs")
GG_mre = m.load_nparray(ana + "GG-ON_mre11_rs.csv")
CT_mre = m.load_nparray(ana + "CT-ON_mre11_rs.csv")
TA_mre = m.load_nparray(ana + "TA-ON_mre11_rs.csv")
GG_cas = m.load_nparray(ana + "GG-ON_cas9_rs.csv")
CT_cas = m.load_nparray(ana + "CT-ON_cas9_rs.csv")
TA_cas = m.load_nparray(ana + "TA-ON_cas9_rs.csv")
m.mergerows([GG_mre, CT_mre, TA_mre], ana + "MERGE-ON_mre_rs.csv")
m.mergerows([GG_cas, CT_cas, TA_cas], ana + "MERGE-ON_cas_rs.csv")

""" Generate peak profiles centered at cut site for only abutting reads (Figure S4) """
m.peak_profile(m.target_gen(alnpath, 1500, AluTA),
               ana + "TA-ON_mre11_rs_abut.bam", ana + "TA-ON_mre11_rs_abut")
m.peak_profile(m.target_gen(alnpath, 1500, AluTA),
               ana + "TA-ON_cas9_rs_abut.bam", ana + "TA-ON_cas9_rs_abut")

""" Generate peak profiles centered at the cut site for all putative on-target sites
    (Figures 1D to 1F) """
m.peak_profile(m.target_gen(alnpath, 1500, AluGG), mreGGbam, ana + "GG-ON_mre11")
m.peak_profile(m.target_gen(alnpath, 1500, AluCT), mreCTbam, ana + "CT-ON_mre11")
m.peak_profile(m.target_gen(alnpath, 1500, AluTA), mreTAbam, ana + "TA-ON_mre11")
m.peak_profile(m.target_gen(alnpath, 1500, AluGG), casGGbam, ana + "GG-ON_cas9")
m.peak_profile(m.target_gen(alnpath, 1500, AluCT), casCTbam, ana + "CT-ON_cas9")
m.peak_profile(m.target_gen(alnpath, 1500, AluTA), casTAbam, ana + "TA-ON_cas9")
GG_mre = m.load_nparray(ana + "GG-ON_mre11_bpeaks.csv")
CT_mre = m.load_nparray(ana + "CT-ON_mre11_bpeaks.csv")
TA_mre = m.load_nparray(ana + "TA-ON_mre11_bpeaks.csv")
GG_cas = m.load_nparray(ana + "GG-ON_cas9_bpeaks.csv")
CT_cas = m.load_nparray(ana + "CT-ON_cas9_bpeaks.csv")
TA_cas = m.load_nparray(ana + "TA-ON_cas9_bpeaks.csv")
m.mergerows([GG_mre, CT_mre, TA_mre], ana + "MERGED-ON_mre11_bpeaks.csv")
m.mergerows([GG_cas, CT_cas, TA_cas], ana + "MERGED-ON_cas9_bpeaks.csv")

""" Determine multiple sequence alignments (up to 300) for each ChIP-seq paired-end reads around
    putative on-target sites (Figure 1G) """
m.find_msa(m.target_gen(alnpath, 750, AluGG), casGGbam, ana + "GG-ON_cas9_msa", hg38)
m.find_msa(m.target_gen(alnpath, 750, AluCT), casCTbam, ana + "CT-ON_cas9_msa", hg38)
m.find_msa(m.target_gen(alnpath, 750, AluTA), casTAbam, ana + "TA-ON_cas9_msa", hg38)
m.find_msa(m.target_gen(alnpath, 1250, AluGG), mreGGbam, ana + "GG-ON_mre11_msa", hg38)
m.find_msa(m.target_gen(alnpath, 1250, AluCT), mreCTbam, ana + "CT-ON_mre11_msa", hg38)
m.find_msa(m.target_gen(alnpath, 1250, AluTA), mreTAbam, ana + "TA-ON_mre11_msa", hg38)
