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
    genome_savepath = "/mnt/c/Users/rzou4/Desktop/"         # Directory to store indexed genome
    hg38 = ['hg38', "/mnt/c/Users/rzou4/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/mnt/c/Users/rzou4/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/mnt/c/Users/rzou4/bioinformatics/mm10_bowtie2/mm10.fa"]
    dr11 = ['dr11', "/mnt/c/Users/rzou4/bioinformatics/dr11_bowtie2/dr11.fa"]
    datadir = "/mnt/z/rzou4/NGS_data/4_damage/"             # Directory for input and output data
elif sys.platform == "darwin":                              # File paths (macOS)
    genome_savepath = "/Users/rogerzou/Desktop/"            # Directory to store indexed genome
    hg38 = ['hg38', "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"]
    hg19 = ['hg19', "/Users/rogerzou/bioinformatics/hg19_bowtie2/hg19.fa"]
    mm10 = ['mm10', "/Users/rogerzou/bioinformatics/mm10_bowtie2/mm10.fa"]
    dr11 = ['dr11', "/Users/rogerzou/bioinformatics/dr11_bowtie2/dr11.fa"]
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
DR1 = "GTGCATTTATGAAGTGTGCTTCACACAGGTGAGTGGGCTTGACAAACCACCTGTAGAAACACTCTTCTCTCTATAAAAAAAAAAAAAAAC" \
      "ACTCCCCCCTCCTTACTCTAGCACTTAATTCTCTGAGCACTAACAGTTCCTTTGTATAATTAGCACTTCTTGTGTGTATTGCCTCTTCTT" \
      "GTTGAATCGCTGAATGCCTCCTCAATTGTAAGTCGCTTTGGACAAAAGCGTCTGCTAAATGACTAAATGTAAATGT"
DR2 = "GTCGGCGCCAATAGCCTAGTGGTTAGTGCGTCGACACATAGCACCGAGGTGCTCGCAGCGACCCGAGTTCGATTCCCGTCTCGAGGTCCT" \
      "TTGCTGATCCTTCCCCTATCTCTGCTCCCCACACTTTCCTGTCTCTATATCTCCACTGTCCTATCAATAAAGGTGAAAACCCCTAAAAAA" \
      "TAAT"
DAN = "GGCGACGCAGTGGCGCAGTAGGTAGCGCTGTCGCCTCACAGCAAGAAGGTCGCTGGTTCGAGCCTCGGCTGGGTCAGTTGGCGTTTCTGT" \
      "GTGGAGTTTGCATGTTCTCCCTGCGTTCGCGTGGGTTTCCTCCGGGTGCTCCGGTTTCCCCCACAGTCCAAAGACATGCGGTACAGGTGA" \
      "ATTGGGTARGCTAAATTGTCCGTAGTGTATGAGTGTGTGTGAATGAGTGTGTATGGGTGTTTCCCAGTGATGGGTTGCGGCTGGAAGGGC" \
      "ATCCGCTGCGTAAAACATATGCTGGATAAGTTGGCGGTTCATTCCGCTGTGGCGACCCCGGATTAATAAAGGGACTAAGCCGAAAAGAAA" \
      "ATGAATGAATGAATGAATRAATTATATAA"

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
psearch_dr11 = ana_1 + "psearch_dr11"


""" ############################################################################################ """
""" ############# Human ############# """
""" Find all sequences in genome with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site.
    For each protospacer, generate & align artificial 2x36bp or 2x75bp PE vs SE ChIP-seq reads. """
""" Alu """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg38 + "_Alu", seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg38 as SAM file
msa.get_targets_bowtie2(psearch_hg38 + "_Alu", hg38[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_hg38 + "_Alu" + ".sam")
msa.get_targets_stats(gen, hg38[0], psearch_hg38 + "_Alu", chromhmm=True)
# from MSA, get distance between each putative target site
msa.get_targets_dist(psearch_hg38 + "_Alu" + "_align.csv", psearch_hg38 + "_Alu")

# generate 2x36bp artificial PE ChIP-seq reads at all potential protospacers (100-300 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_36bp_100-300", hg38[0],
                           genome_savepath, rlen=36, ct_min=100, ct_max=300)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_100-300_2_msa")
# generate 2x36bp artificial PE ChIP-seq reads at all potential protospacers (050-099 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_36bp_050-099", hg38[0],
                           genome_savepath, rlen=36, ct_min=50, ct_max=99)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_050-099_2_msa")
# generate 2x36bp artificial PE ChIP-seq reads at all potential protospacers (005-049 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_36bp_005-049", hg38[0],
                           genome_savepath, rlen=36, ct_min=5, ct_max=49)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_36bp_005-049_2_msa")


# generate 2x75bp artificial PE ChIP-seq reads at all potential protospacers (100-300 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_75bp_100-300", hg38[0],
                           genome_savepath, rlen=75, ct_min=100, ct_max=300)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_100-300_2_msa")
# generate 2x75bp artificial PE ChIP-seq reads at all potential protospacers (050-099 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_75bp_050-099", hg38[0],
                           genome_savepath, rlen=75, ct_min=50, ct_max=99)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_050-099_2_msa")
# generate 2x75bp artificial PE ChIP-seq reads at all potential protospacers (005-049 targets)
msa.get_artifical_pe_reads(msa.gen_putative(psearch_hg38 + "_Alu" + ".sam"),
                           ana_2 + "psearch_hg38_Alu_PE_75bp_005-049", hg38[0],
                           genome_savepath, rlen=75, ct_min=5, ct_max=49)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049", hg38[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_msa")
# align artificial ChIP-seq reads to genome with SE alignment
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_1", hg38[1])
msa.bowtie2_msa_single(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_2", hg38[1])
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_1_msa")
msa.parse_msa_sam_single(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_2_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_1_msa")
msa.get_msa_stats(ana_2 + "psearch_hg38_Alu_PE_75bp_005-049_2_msa")


""" ############################################################################################ """
""" ############# Mouse ############# """
""" Find all sequences in genome with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site.
    For each protospacer, generate & align artificial 2x75bp PE vs SE ChIP-seq reads. """
""" Alu """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_mm10 + "_Alu", seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in mm10 as SAM file
msa.get_targets_bowtie2(psearch_mm10 + "_Alu", mm10[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_mm10 + "_Alu" + ".sam")
msa.get_targets_stats(gen, mm10[0], psearch_mm10 + "_Alu")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_mm10 + "_Alu" + ".sam"),
                           ana_2 + "psearch_mm10_Alu_PE_75bp", mm10[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_mm10_Alu_PE_75bp", mm10[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_mm10_Alu_PE_75bp_msa")
msa.get_msa_stats(ana_2 + "psearch_mm10_Alu_PE_75bp_msa")

""" B4 """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_mm10 + "_B4", seqstr=B4, numbases=9)
# from FASTA, MSA up to 1000 locations in mm10 as SAM file
msa.get_targets_bowtie2(psearch_mm10 + "_B4", mm10[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_mm10 + "_B4" + ".sam")
msa.get_targets_stats(gen, mm10[0], psearch_mm10 + "_B4")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_mm10 + "_B4" + ".sam"),
                           ana_2 + "psearch_mm10_B4_PE_75bp", mm10[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_mm10_B4_PE_75bp", mm10[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_mm10_B4_PE_75bp_msa")
msa.get_msa_stats(ana_2 + "psearch_mm10_B4_PE_75bp_msa")


""" ############################################################################################ """
""" ############# Zebrafish ############# """
""" Find all sequences in genome with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites and epigenetic characteristics of each site.
    For each protospacer, generate and align artificial 2x75bp PE ChIP-seq reads. """
""" DR1 """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_dr11 + "_DR1", seqstr=DR1, numbases=9)
# from FASTA, MSA up to 1000 locations in dr11 as SAM file
msa.get_targets_bowtie2(psearch_dr11 + "_DR1", dr11[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_dr11 + "_DR1" + ".sam")
msa.get_targets_stats(gen, dr11[0], psearch_dr11 + "_DR1")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_dr11 + "_DR1" + ".sam"),
                           ana_2 + "psearch_dr11_DR1_PE_75bp", dr11[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_dr11_DR1_PE_75bp", dr11[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_dr11_DR1_PE_75bp_msa")
msa.get_msa_stats(ana_2 + "psearch_dr11_DR1_PE_75bp_msa")

""" DR2 """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_dr11 + "_DR2", seqstr=DR2, numbases=9)
# from FASTA, MSA up to 1000 locations in dr11 as SAM file
msa.get_targets_bowtie2(psearch_dr11 + "_DR2", dr11[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_dr11 + "_DR2" + ".sam")
msa.get_targets_stats(gen, dr11[0], psearch_dr11 + "_DR2")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_dr11 + "_DR2" + ".sam"),
                           ana_2 + "psearch_dr11_DR2_PE_75bp", dr11[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_dr11_DR2_PE_75bp", dr11[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_dr11_DR2_PE_75bp_msa")
msa.get_msa_stats(ana_2 + "psearch_dr11_DR2_PE_75bp_msa")

""" DANA """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_dr11 + "_DAN", seqstr=DR2, numbases=9)
# from FASTA, MSA up to 1000 locations in dr11 as SAM file
msa.get_targets_bowtie2(psearch_dr11 + "_DAN", dr11[1])
# from SAM, summarize MSA (including gene + epigenetic status)
gen = msa.gen_putative(psearch_dr11 + "_DAN" + ".sam")
msa.get_targets_stats(gen, dr11[0], psearch_dr11 + "_DAN")

# generate 2x75bp artificial paired-end ChIP-seq reads at all potential protospacer sequences
msa.get_artifical_pe_reads(msa.gen_putative(psearch_dr11 + "_DAN" + ".sam"),
                           ana_2 + "psearch_dr11_DAN_PE_75bp", dr11[0], genome_savepath, rlen=75)
# align artificial ChIP-seq reads to genome with PE alignment
msa.bowtie2_msa_paired(ana_2 + "psearch_dr11_DAN_PE_75bp", dr11[1])
msa.parse_msa_sam_paired(ana_2 + "psearch_dr11_DAN_PE_75bp_msa")
msa.get_msa_stats(ana_2 + "psearch_dr11_DAN_PE_75bp_msa")


""" ############################################################################################ """
""" Starting from Alu, find all sequences in hg19 with at most 3 mismatches within 9 bases from PAM.
    Determine putative on-target genomic sites and epigenetic characteristics of a subset of sites:
    AluGG, AluTA, and AluCT. """
# get list of all protospacer sequences as FASTA file
msa.get_targets_fasta(psearch_hg19 + "_Alu", seqstr=Alu, numbases=9)
# from FASTA, MSA up to 1000 locations in hg19 as SAM file
msa.get_targets_bowtie2(psearch_hg19 + "_Alu", hg19[1])
# from SAM, summarize MSA (including gene + epigenetic status)
AluGG, AluCT, AluTA = "CCTGTAGTCCCAGCTACTGG", "CCTGTAGTCCCAGCTACTCT", "CCTGTAGTCCCAGCTACTTA"
gen = msa.gen_putative(psearch_hg19 + "_Alu" + ".sam", subset=[AluGG, AluCT, AluTA])
msa.get_targets_stats(gen, hg19[0], psearch_hg19 + "_subset")


""" ############################################################################################ """
""" Determine the local genomic sequence for all putative on-targets for a subset of gRNAs. """

AluGG, AluCT, AluTA = "CCTGTAGTCCCAGCTACTGG", "CCTGTAGTCCCAGCTACTCT", "CCTGTAGTCCCAGCTACTTA"
gen = msa.gen_putative(psearch_hg38 + "_Alu" + ".sam", subset=[AluGG, AluCT, AluTA])
msa.get_target_sequences(gen, ana_3 + "psearch_hg38_Alu_subseq", hg38[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "psearch_hg38_Alu_subseq.csv", AluGG, ana_3 + "parse_GG")
msa.parse_target_sequences(ana_3 + "psearch_hg38_Alu_subseq.csv", AluCT, ana_3 + "parse_CT")
msa.parse_target_sequences(ana_3 + "psearch_hg38_Alu_subseq.csv", AluTA, ana_3 + "parse_TA")
msa.parse_imshow(ana_3 + "parse_GG_num")
msa.parse_imshow(ana_3 + "parse_CT_num")
msa.parse_imshow(ana_3 + "parse_TA_num")

mm10_Alu_025t = "CTGTAATCCCAGCACTTTGG"
mm10_Alu_063t = "ACTCGGGAGGCAGAGGAAGG"
mm10_Alu_125t = "GTGCTGGGATTAAAGGTGTG"
ss_mm10_Alu = [mm10_Alu_025t, mm10_Alu_063t, mm10_Alu_125t]
gen = msa.gen_putative(psearch_mm10 + "_Alu" + ".sam", subset=ss_mm10_Alu)
msa.get_target_sequences(gen, ana_3 + "ts_mm10_Alu", mm10[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "ts_mm10_Alu.csv", mm10_Alu_025t, ana_3 + "ss_mm10_Alu_025t")
msa.parse_target_sequences(ana_3 + "ts_mm10_Alu.csv", mm10_Alu_063t, ana_3 + "ss_mm10_Alu_063t")
msa.parse_target_sequences(ana_3 + "ts_mm10_Alu.csv", mm10_Alu_125t, ana_3 + "ss_mm10_Alu_125t")
msa.parse_imshow(ana_3 + "ss_mm10_Alu_025t_num")
msa.parse_imshow(ana_3 + "ss_mm10_Alu_063t_num")
msa.parse_imshow(ana_3 + "ss_mm10_Alu_125t_num")

mm10_B4_035t = "ACAGGGTCTCTCATTGAACC"
mm10_B4_063t = "TGAACTCAGGTCATCAGACT"
mm10_B4_158t = "ACAGGGTCTCTCACTGAACC"
ss_mm10_B4 = [mm10_B4_035t, mm10_B4_063t, mm10_B4_158t]
gen = msa.gen_putative(psearch_mm10 + "_B4" + ".sam", subset=ss_mm10_B4)
msa.get_target_sequences(gen, ana_3 + "ts_mm10_B4", mm10[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "ts_mm10_B4.csv", mm10_B4_035t, ana_3 + "ss_mm10_B4_035t")
msa.parse_target_sequences(ana_3 + "ts_mm10_B4.csv", mm10_B4_063t, ana_3 + "ss_mm10_B4_063t")
msa.parse_target_sequences(ana_3 + "ts_mm10_B4.csv", mm10_B4_158t, ana_3 + "ss_mm10_B4_158t")
msa.parse_imshow(ana_3 + "ss_mm10_B4_035t_num")
msa.parse_imshow(ana_3 + "ss_mm10_B4_063t_num")
msa.parse_imshow(ana_3 + "ss_mm10_B4_158t_num")

dr11_DR1_027t = "TTATGAAGTGTGCCTCACAC"
dr11_DR1_066t = "GTCAAGCCCACTCACCTGCA"
dr11_DR1_132t = "GAGAAGAGTGTTTCTACAGG"
ss_dr11_DR1 = [dr11_DR1_027t, dr11_DR1_066t, dr11_DR1_132t]
gen = msa.gen_putative(psearch_dr11 + "_DR1" + ".sam", subset=ss_dr11_DR1)
msa.get_target_sequences(gen, ana_3 + "ts_dr11_DR1", dr11[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "ts_dr11_DR1.csv", dr11_DR1_027t, ana_3 + "ss_dr11_DR1_027t")
msa.parse_target_sequences(ana_3 + "ts_dr11_DR1.csv", dr11_DR1_066t, ana_3 + "ss_dr11_DR1_066t")
msa.parse_target_sequences(ana_3 + "ts_dr11_DR1.csv", dr11_DR1_132t, ana_3 + "ss_dr11_DR1_132t")
msa.parse_imshow(ana_3 + "ss_dr11_DR1_027t_num")
msa.parse_imshow(ana_3 + "ss_dr11_DR1_066t_num")
msa.parse_imshow(ana_3 + "ss_dr11_DR1_132t_num")

dr11_DR2_024t = "TCAGCAAAGGACCTTGAGCC"
dr11_DR2_044t = "CGAGTTCGATTCCCGGCTCG"
dr11_DR2_125t = "CGCACTAACCACTAGGCTAT"
ss_dr11_DR2 = [dr11_DR2_024t, dr11_DR2_044t, dr11_DR2_125t]
gen = msa.gen_putative(psearch_dr11 + "_DR2" + ".sam", subset=ss_dr11_DR2)
msa.get_target_sequences(gen, ana_3 + "ts_dr11_DR2", dr11[0], genome_savepath, win=500)
msa.parse_target_sequences(ana_3 + "ts_dr11_DR2.csv", dr11_DR2_024t, ana_3 + "ss_dr11_DR2_024t")
msa.parse_target_sequences(ana_3 + "ts_dr11_DR2.csv", dr11_DR2_044t, ana_3 + "ss_dr11_DR2_044t")
msa.parse_target_sequences(ana_3 + "ts_dr11_DR2.csv", dr11_DR2_125t, ana_3 + "ss_dr11_DR2_125t")
msa.parse_imshow(ana_3 + "ss_dr11_DR2_024t_num")
msa.parse_imshow(ana_3 + "ss_dr11_DR2_044t_num")
msa.parse_imshow(ana_3 + "ss_dr11_DR2_125t_num")
