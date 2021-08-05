"""
Script for:
(1) Identifying potential gRNAs from Alu repetitive sequence.
(2) Determining putative genome-wide on-target sites for each gRNA.
(3) Determining the epigenetic context and predicted percentage of ambiguous reads for each gRNA.
(4) Determining the nucleotide composition around on-target sites for GG, CT, and TA gRNAs.
"""

import src.msa as msa
import sys
import os

""" Determine run paths based on operating system """
if sys.platform == "linux" or sys.platform == "linux2":     # File paths (Ubuntu)
    genome_savepath = "/mnt/c/Users/Roger/Desktop/"         # Directory to store indexed genome
    hg38 = ['hg38', "/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/mnt/c/Users/Roger/bioinformatics/mm10_bowtie2/mm10.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    genome_savepath = "/Users/rogerzou/Desktop/"            # Directory to store indexed genome
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/Users/rogerzou/bioinformatics/mm10_bowtie2/mm10.fa"]
    datadir = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"  # Directory for input and output data
else:
    sys.exit()

""" Sequences """
Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTG" \
      "GCCAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC" \
      "TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG" \
      "ACTCCGTCTC"
B4 = "GGGCTGGGGAGGATGGCTCAGTCGGTAAAGTGCTTGCTGTGCAAGCATGAGGACCTGAGTTCAGATCCCCAGAACCCACATATAAAAAAGC" \
     "CAGGCATGGTGGTATGTGCTTGTAATCCCAGYGCTGGGGAGGCAGAGACAGGAGGATCCCTGGGGCTCGCTGGCCAGCCAGCCTAGCCTAA" \
     "TTGGTGAGCTCCAGGTTCAGTGAGAGACCCTGTCTCAAAAAATAAGGTGGAGAGTAACTGAGGAAGACACCTGAGGTTGACCTCTGGCCTC" \
     "CACATACACACACACACACACACACA"

""" Set analysis path """
ana = datadir + "Alu_ana_1_putative/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_protosearch/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_artificial/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_sequences/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
psearch_hg38 = ana_1 + "psearch_hg38"
psearch_hg19 = ana_1 + "psearch_hg19"
psearch_mm10 = ana_1 + "psearch_mm10"


""" ############################################################################################ """
""" Starting from Alu, find all sequences in hg38 with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, the epigenetic characteristics of each site,
    and distances between adjacent sites. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg38, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg38 as SAM file
msa.get_targets_bowtie2(psearch_hg38, hg38[1])
# from SAM, summarize MSA (including gene + epigenetic status)
msa.get_targets_stats(msa.gen_putative(psearch_hg38 + ".sam"), hg38[0], psearch_hg38, chromhmm=True)
# from MSA, get distance between each putative target site
msa.get_targets_dist(psearch_hg38 + "_align.csv", psearch_hg38)

""" ############################################################################################ """
""" For each potential protospacer sequence from Alu, generate artificial paired-end ChIP-seq reads.
    Determine number of reads that (1) have an optimal alignment or (2) have multiple optimal
    alignments. Compare alignments using paired-end information vs single-end information alone. """
# generate 2x36bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_36bp", hg38[0], genome_savepath, rlen=36)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_36bp", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_msa")  # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_msa")         # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_1_msa")       # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_2_msa")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_75bp", hg38[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_75bp", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_msa")  # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_msa")         # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_1_msa")       # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_2_msa")


""" ############################################################################################ """
""" Starting from Alu, find all sequences in hg19 with at most 3 mismatches within 9 bases from PAM.
    Determine putative on-target genomic sites and epigenetic characteristics of a subset of sites:
    AluGG, AluTA, and AluCT. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg19, seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg19 as SAM file
msa.get_targets_bowtie2(psearch_hg19, hg19[1])
# from SAM, summarize MSA (including gene + epigenetic status)
AluGG, AluCT, AluTA = "CCTGTAGTCCCAGCTACTGG", "CCTGTAGTCCCAGCTACTCT", "CCTGTAGTCCCAGCTACTTA"
gen = msa.gen_putative(psearch_hg19 + ".sam", subset=[AluGG, AluCT, AluTA])
msa.get_targets_stats(gen, hg19[0], psearch_hg19 + "_ontargets")


""" ############################################################################################ """
""" Determine the local genomic sequence for all putative on-targets for a subset of gRNAs:
    AluGG, AluTA, and AluCT. """
AluGG, AluCT, AluTA = "CCTGTAGTCCCAGCTACTGG", "CCTGTAGTCCCAGCTACTCT", "CCTGTAGTCCCAGCTACTTA"
gen = msa.gen_putative(psearch_hg19 + ".sam", subset=[AluGG, AluCT, AluTA])
msa.get_target_sequences(gen, ana_3 + "psearch_hg19_Alu_GG-CT-TA_seq",
                         hg19[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "psearch_hg19_Alu_GG-CT-TA_seq.csv",
                           AluGG, ana_3 + "parse_GG")
msa.parse_target_sequences(ana_3 + "psearch_hg19_Alu_GG-CT-TA_seq.csv",
                           AluCT, ana_3 + "parse_CT")
msa.parse_target_sequences(ana_3 + "psearch_hg19_Alu_GG-CT-TA_seq.csv",
                           AluTA, ana_3 + "parse_TA")
msa.parse_imshow(ana_3 + "parse_GG_num")
msa.parse_imshow(ana_3 + "parse_CT_num")
msa.parse_imshow(ana_3 + "parse_TA_num")


""" ############################################################################################ """
""" Starting from Alu, find all sequences in mm10 with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site.
    For each potential protospacer sequence from Alu, generate artificial paired-end ChIP-seq reads.
    Determine number of optimal alignments with 2x36bp PE reads. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_mm10 + "_Alu", seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in mm10 as SAM file
msa.get_targets_bowtie2(psearch_mm10 + "_Alu", mm10[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_mm10 + "_Alu" + ".sam")
msa.get_targets_stats(gen, mm10[0], psearch_mm10 + "_Alu" + "_ontargets")

# generate 2x36bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_mm10 + "_Alu" + ".sam"),
                           ana_2 + "psearch_mm10_Alu_PE_36bp", mm10[0], genome_savepath, rlen=36)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_mm10_Alu_PE_36bp", mm10[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_mm10_Alu_PE_36bp_msa")  # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_mm10_Alu_PE_36bp_msa")         # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_Alu_PE_36bp_1", mm10[1])
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_Alu_PE_36bp_2", mm10[1])
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_Alu_PE_36bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_Alu_PE_36bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_mm10_Alu_PE_36bp_1_msa")       # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_mm10_Alu_PE_36bp_2_msa")


""" ############################################################################################ """
""" Starting from B4, find all sequences in mm10 with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site.
    For each potential protospacer sequence from Alu, generate artificial paired-end ChIP-seq reads.
    Determine number of optimal alignments with 2x36bp PE reads. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_mm10 + "_B4", seqstr=B4, numbases=9)
# from FASTA, MSA up to 1000 locations in mm10 as SAM file
msa.get_targets_bowtie2(psearch_mm10 + "_B4", mm10[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_mm10 + "_B4" + ".sam")
msa.get_targets_stats(gen, mm10[0], psearch_mm10 + "_B4" + "_ontargets")

# generate 2x36bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_mm10 + "_B4" + ".sam"),
                           ana_2 + "psearch_mm10_B4_PE_36bp", mm10[0], genome_savepath, rlen=36)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_mm10_B4_PE_36bp", mm10[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_mm10_B4_PE_36bp_msa")   # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_36bp_msa")          # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_B4_PE_36bp_1", mm10[1])
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_B4_PE_36bp_2", mm10[1])
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_B4_PE_36bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_B4_PE_36bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_36bp_1_msa")        # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_36bp_2_msa")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_mm10 + "_B4" + ".sam"),
                           ana_2 + "psearch_mm10_B4_PE_75bp", mm10[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_mm10_B4_PE_75bp", mm10[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_mm10_B4_PE_75bp_msa")   # output CSV file of PE alignments
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_75bp_msa")          # get statistics for PE alignment
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_B4_PE_75bp_1", mm10[1])
msa.bowtie2_msa_single(ana_2 + "psearch_mm10_B4_PE_75bp_2", mm10[1])
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_B4_PE_75bp_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_mm10_B4_PE_75bp_2_msa")
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_75bp_1_msa")        # get statistics for SE alignment
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_75bp_2_msa")
