"""
Script for:
(1) Computational identification of potential gRNAs from Alu repetitive sequence.
(2) Determination of putative genome-wide on-target sites for each gRNA.
(3) Determination of epigenetic context and on-target ChIP-seq read uniqueness for each gRNA.
"""

import src.msa as msa
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu - Lab computer)
    desktop = "/mnt/c/Users/Roger/Desktop/"
    hg38path = "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"
    hg19path = "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"
    labhome = "/mnt/z/rzou4/NGS_data/4_damage/"
elif sys.platform == "darwin":                              # File paths (Mac - Personal computer)
    desktop = "/Users/rogerzou/Desktop/"
    hg38path = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
    hg19path = "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"
    labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
else:
    sys.exit()

""" Sequences """
Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTG" \
      "GCCAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC" \
      "TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG" \
      "ACTCCGTCTC"

""" Set analysis path """
ana = labhome + "Alu_ana_1_putative/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_protosearch/"
psearch_hg38 = ana_1 + "psearch_hg38"
psearch_hg19 = ana_1 + "psearch_hg19"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_artificial/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" ############################################################################################ """
""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, the epigenetic characteristics of each site,
    and distances between adjacent sites.
    (Figures 1A to 1C) """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg38, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg38 as SAM file
msa.get_targets_bowtie2(psearch_hg38, hg38path)
# from SAM, summarize MSA (including gene + epigenetic status)
msa.get_targets_stats(msa.gen_putative(psearch_hg38 + ".sam"), psearch_hg38)
# from MSA, get distance between each putative target site
msa.get_targets_dist(psearch_hg38 + "_align.csv", psearch_hg38)

""" ############################################################################################ """
""" For each potential protospacer sequence, generate artificial paired-end ChIP-seq reads.
    Determine number of reads that (1) align once, (2) align multiple times with one optimal
    alignment, and (3) align multiple times with multiple optimal alignments.
    Compare alignments using paired-end information vs single-end information alone.
    (Figures 1D-1E)"""
# generate artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + ".sam"), ana_2 + "psearch_PE", desktop)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_PE", hg38path)
msa.parse_msa_sam_paired(ana_2 + "psearch_PE_msa")      # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_PE_msa")             # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_PE_1", hg38path)
msa.bowtie2_msa_single(ana_2 + "psearch_PE_2", hg38path)
msa.parse_msa_sam_single(ana_2 + "psearch_PE_1_msa")    # output CSV of SE alignments for r1 and r2
msa.parse_msa_sam_single(ana_2 + "psearch_PE_2_msa")
msa.get_msa_stats(ana_2 + "psearch_PE_1_msa")           # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_PE_2_msa")


""" ############################################################################################ """
""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites for HG19! """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg19, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg19 as SAM file
msa.get_targets_bowtie2(psearch_hg19, hg19path)
