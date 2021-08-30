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
    dr11 = ['dr11', "/mnt/c/Users/Roger/bioinformatics/dr11_bowtie2/dr11.fa"]
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

""" Set analysis path """
ana = datadir + "Alu_ana_10_lineage/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_search/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_primers/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None


""" Set file paths """
psearch_hg38_Alu = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_Alu"
psearch_dr11_DR1 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR1"
psearch_dr11_DR2 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR2"
psearch_dr11_DAN = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DAN"

""" ############################################################################################ """
""" For a given number (or sets of numbers) of targets, do the following:
1) find all the gRNAs that have that many targets
2) find the alignments for all the targets of each of those gRNAs
3) find nested PCR primers for all those alignments, using primer3
4) filter the primers found according to how uniquely they bind to the genome, using bowtie2.
"""

# Find PCR primers for target sites of all danRer11 DR1 mgRNAs that target between 10-40 sites
outfile = ana_1 + "lineage_dr11_DR1_10-40"
gen = msa.gen_putative(psearch_dr11_DR1 + ".sam")
ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])

# Find PCR primers for target sites of all danRer11 DR2 mgRNAs that target between 10-40 sites
outfile = ana_1 + "lineage_dr11_DR2_10-40"
gen = msa.gen_putative(psearch_dr11_DR2 + ".sam")
ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])

# Find PCR primers for target sites of all danRer11 DANA mgRNAs that target between 10-40 sites
outfile = ana_1 + "lineage_dr11_DAN_10-40"
gen = msa.gen_putative(psearch_dr11_DAN + ".sam")
ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])

# Find PCR primers for target sites of all hg38 Alu mgRNAs that target between 10-40 sites
outfile = ana_1 + "lineage_hg38_Alu_10-40"
gen = msa.gen_putative(psearch_hg38_Alu + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath, ct_values=[*range(10, 41)])
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])


""" ############################################################################################ """
dr11_DR1_1 = "TGGGCTTGACAAACCACCTGNGG"
dr11_DR1_2 = "GTCTGCTAAATGACGTCATGNGG"
inpath = ana_1 + "lineage_dr11_DR1_10-40"
ltr.get_nested_primers(inpath + "_inn_msa", inpath + "_out_msa", ana_2 + "dr11_DR1_1", dr11_DR1_1)
ltr.get_nested_primers(inpath + "_inn_msa", inpath + "_out_msa", ana_2 + "dr11_DR1_2", dr11_DR1_2)
