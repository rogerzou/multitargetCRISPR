"""
Script for:
(1) Identifying potential gRNAs from Alu repetitive sequence.
(2) Determining putative genome-wide on-target sites for each gRNA.
(3) Determining the epigenetic context and predicted percentage of ambiguous reads for each gRNA.
(4) Determining the nucleotide composition around on-target sites for GG, CT, and TA gRNAs.
"""

import src.msa as msa
import src.ltr as ltr
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
ana = datadir + "Alu_ana_10_lineage/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_search/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None

""" Set file paths """
psearch_hg38 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38"


""" ############################################################################################ """
""" For a given number (or sets of numbers) of targets, do the following:
1) find all the gRNAs that have that many targets
2) find the alignments for all the targets of each of those gRNAs
3) find nested PCR primers for all those alignments, using primer3
4) filter the primers found according to how uniquely they bind to the genome, using bowtie2.
"""

# Find nested PCR primers for all the alignments of all the mgRNAs that have "num_targets" targets
num_targets = 15
outfile = ana_1 + "targets_ct" + str(num_targets)
gen = msa.gen_putative(psearch_hg38 + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath, ct_values=[num_targets])

# Run wrapper function that runs bowtie to align primer pairs to genome, parses SAM output, then
# obtains statistics on the best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])

# # Run bowtie2 to align pairs of primers to the genome in order to check for uniqueness.
# # Check outer PCR primers (PCR-1)
# ltr.bowtie2_msa_primers(outfile + "num_targets_out", hg38[1])
# # Check inner PCR primers (PCR-2)
# ltr.bowtie2_msa_primers(outfile + "num_targets_inn", hg38[1])
# # parse sam file to output info about outer PCR primers
# msa.parse_msa_sam_paired(outfile + "num_targets_out" + "_msa")
# # parse sam file to output info about inner PCR primers
# msa.parse_msa_sam_paired(outfile + "num_targets_inn" + "_msa")
# # get statistics for each gRNA
# ltr.get_stats_primers(outfile + "num_targets")
