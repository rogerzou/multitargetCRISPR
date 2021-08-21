"""
Script for:
(1) Finding gRNAs suitable for lineage tracing.
"""

import src.msa as msa
import sys
import os
import src.ltr as ltr

""" Define Run Paths """
genome_savepath = "/mnt/d/DATA/MULTITARGET_GRNA_DESIGN/"         # Directory to store indexed genome
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]
hg19 = ['hg19', "/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19.fa"]
mm10 = ['mm10', "/mnt/d/DATABASES/MOUSE_REF_GENOME/MM10_B/mm10.fa"]
datadir = "/mnt/d/DATA/MULTITARGET_GRNA_DESIGN/"             # Directory for input and output data
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
""" ############################################################################################ """
""" Starting from Alu, find all sequences in hg38 with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, the epigenetic characteristics of each site,
    and distances between adjacent sites. """
# get list of all protospacer sequences as FASTA file (only if this hasn't been done before)
msa.get_targets_fasta(psearch_hg38, seqstr=Alu, numbases=9) \
    if not os.path.exists(psearch_hg38 + ".fa") else None
# from FASTA, MSA up to 1000 locations in hg38 as SAM file
msa.get_targets_bowtie2(psearch_hg38, hg38[1]) \
    if not os.path.exists(psearch_hg38 + ".sam") else None
# from SAM, build generator with info about the alignments
gen = msa.gen_putative(psearch_hg38 + ".sam")
""" ############################################################################################ """
""" For a given number (or sets of numbers) of targets, do the following:
1) find all the gRNAs that have that many targets
2) find the alignments for all the targets of each of those gRNAs
3) find nested PCR primers for all those alignments, using primer3
4) filter the primers found according to how uniquely they bind to the genome, using bowtie2.
"""
# Select the number of targets desired for the mgRNA and define the output file
num_targets = 20
outfile = ana_3 + "targets_ct" + str(num_targets)
# Find nested PCR primers for all the alignments of all the mgRNAs that have "num_targets" targets
# ltr.get_primers_nested(gen, outfile + "num_targets", hg38[0], genome_savepath, ct_values=[num_targets])
# Run bowtie2 to align pairs of primers to the genome in order to check for uniqueness.
ltr.bowtie2_msa_primers(outfile + "_out", hg38[1], k_count=10) # Check outer PCR primers (PCR-1)
ltr.bowtie2_msa_primers(outfile + "_in", hg38[1], k_count=10) # Check inner PCR primers (PCR-2)
ltr.parse_sam_primers(outfile + "_out" + "_msa", label=str(0)) # parse sam file to output info about outer PCR primers
ltr.parse_sam_primers(outfile + "_in" + "_msa", label=str(1)) # parse sam file to output info about inner PCR primers
ltr.get_stats_primers(outfile) # get statistics for each gRNA