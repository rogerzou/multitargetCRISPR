"""
Script for:
(1) Computational identification of potential gRNAs from Alu repetitive sequence.
(2) Determination of putative genome-wide on-target sites for each gRNA.
(3) Determination of epigenetic context and predicted percentage of ambiguous reads for each gRNA.
"""

import src.msa as msa
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    genome_savepath = "/mnt/c/Users/Roger/Desktop/"         # Directory to store indexed genome
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    genome_savepath = "/Users/rogerzou/Desktop/"            # Directory to store indexed genome
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" Sequences """
Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTG" \
      "GCCAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC" \
      "TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG" \
      "ACTCCGTCTC"

""" Set analysis path """
ana = datadir + "Alu_ana_1_putative/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_protosearch/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_artificial/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
psearch_hg38 = ana_1 + "psearch_hg38"
psearch_hg19 = ana_1 + "psearch_hg19"


""" ############################################################################################ """
""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, the epigenetic characteristics of each site,
    and distances between adjacent sites.
    (Fig. 1A-1C, S1A) """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg38, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg38 as SAM file
msa.get_targets_bowtie2(psearch_hg38, hg38[1])
# from SAM, summarize MSA (including gene + epigenetic status)
msa.get_targets_stats(msa.gen_putative(psearch_hg38 + ".sam"), 'hg38', psearch_hg38)
# from MSA, get distance between each putative target site
msa.get_targets_dist(psearch_hg38 + "_align.csv", psearch_hg38)


""" ############################################################################################ """
""" For each potential protospacer sequence, generate artificial paired-end ChIP-seq reads.
    Determine number of reads that (1) align once, (2) align multiple times with one optimal
    alignment, and (3) align multiple times with multiple optimal alignments.
    Compare alignments using paired-end information vs single-end information alone.
    (Fig. 1D, S1B-D) """
# generate 2x36bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads_hg38(msa.gen_putative(psearch_hg38 + ".sam"),
                                ana_2 + "psearch_PE_36bp", genome_savepath, rlen=36)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_PE_36bp", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_PE_36bp_msa")      # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_PE_36bp_msa")             # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_PE_36bp_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_PE_36bp_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_PE_36bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_PE_36bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_PE_36bp_1_msa")           # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_PE_36bp_2_msa")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads_hg38(msa.gen_putative(psearch_hg38 + ".sam"),
                                ana_2 + "psearch_PE_75bp", genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_PE_75bp", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_PE_75bp_msa")      # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_PE_75bp_msa")             # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_PE_75bp_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_PE_75bp_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_PE_75bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_PE_75bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_PE_75bp_1_msa")           # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_PE_75bp_2_msa")


""" ############################################################################################ """
""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site for 
    hg19. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg19, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg19 as SAM file
msa.get_targets_bowtie2(psearch_hg19, hg19[1])
# from SAM, summarize MSA (including gene + epigenetic status)
AluGG, AluCT, AluTA = "CCTGTAGTCCCAGCTACTGG", "CCTGTAGTCCCAGCTACTCT", "CCTGTAGTCCCAGCTACTTA"
gen = msa.gen_putative(psearch_hg19 + ".sam", subset=[AluGG, AluCT, AluTA])
msa.get_targets_stats(gen, 'hg19', psearch_hg19 + "_ontargets")
