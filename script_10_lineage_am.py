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

dr11 = ['dr11', "/mnt/d/DATABASES/dr11_bowtie2/dr11.fa"]
datadir = "/mnt/d/DATA/MULTITARGET_GRNA_DESIGN/"             # Directory for input and output data
genome_savepath = "/mnt/d/DATA/MULTITARGET_GRNA_DESIGN/"         # Directory to store indexed genome
hg38 = ['hg38', "/mnt/d/DATABASES/HUMAN_REF_GENOME/GRCh38_noalt_as/GRCh38_noalt_as.fa"]

""" Set analysis path """
ana = datadir + "Alu_ana_10_lineage/"
os.makedirs(ana) if not os.path.exists(ana) else None
ana_1 = ana + "1_search/"
os.makedirs(ana_1) if not os.path.exists(ana_1) else None
ana_2 = ana + "2_primers/"
os.makedirs(ana_2) if not os.path.exists(ana_2) else None
ana_3 = ana + "3_ampngs/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_chip = datadir + "Alu_ana_1_putative/1_protosearch/"
psearch_hg38 = ana_chip + "psearch_hg38"
psearch_out = ana_1 + "psearch_hg38_"

""" Set file paths
psearch_hg38_Alu = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38"
psearch_dr11_DR1 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR1"
psearch_dr11_DR2 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR2"
psearch_dr11_DAN = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DAN"
"""
lineageNGS_210830 = datadir + "210830_lineageNGS/"
lineageNGS_210908 = datadir + "210908_lineageNGS/DATA_ANA/"

# """ ############################################################################################ """
# """ For a given number (or sets of numbers) of targets, do the following:
# 1) find all the gRNAs that have that many targets
# 2) find the alignments for all the targets of each of those gRNAs
# 3) find nested PCR primers for all those alignments, using primer3
# 4) filter the primers found according to how uniquely they bind to the genome, using bowtie2.
# """
# # Find PCR primers for target sites of all danRer11 DR1 mgRNAs that target between 10-40 sites
# outfile = ana_1 + "lineage_dr11_DR1_10-40"
# gen = msa.gen_putative(psearch_dr11_DR1 + ".sam")
# ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# # Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
# ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])
#
# # Find PCR primers for target sites of all danRer11 DR2 mgRNAs that target between 10-40 sites
# outfile = ana_1 + "lineage_dr11_DR2_10-40"
# gen = msa.gen_putative(psearch_dr11_DR2 + ".sam")
# ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# # Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
# ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])
#
# # Find PCR primers for target sites of all danRer11 DANA mgRNAs that target between 10-40 sites
# outfile = ana_1 + "lineage_dr11_DAN_10-40"
# gen = msa.gen_putative(psearch_dr11_DAN + ".sam")
# ltr.get_primers_nested(gen, outfile, dr11[0], genome_savepath, ct_values=[*range(10, 41)])
# # Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
# ltr.bowtie_parse_stats_wrapper(outfile, dr11[1])
"""
# Find PCR primers for target sites of all hg38 Alu mgRNAs that target between 10-40 sites
outfile = ana_1 + "lineage_hg38_Alu"
gen = msa.gen_putative(psearch_hg38 + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath, ct_values=[10])
# # Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])
# ###########################################################################################
# Get targets stats for the best 10-target-mgRNAs
ct10S1 = "CCAGGCTGGAGTGCAGTGCT"
ct10S2 = "AGCTACTCGGGAGGCCAAAG"
ct10S3 = "TGCACTCCAGCCTGGATAAC"
ct10S4 = "GGAGTGCAGTGGCAATATCT"
ct10S5 = "TTTCACCGTGTTGGCCAGCA"
ct10S6 = "AAGTGCTGGGATTATAGGCC"
gen = msa.gen_putative(psearch_hg38 + ".sam", subset=[ct10S3])
msa.get_targets_stats(gen, hg38[0], psearch_out + "ct10S3", chromhmm=True)
gen = msa.gen_putative(psearch_hg38 + ".sam", subset=[ct10S4])
msa.get_targets_stats(gen, hg38[0], psearch_out + "ct10S4", chromhmm=True)
gen = msa.gen_putative(psearch_hg38 + ".sam", subset=[ct10S5])
msa.get_targets_stats(gen, hg38[0], psearch_out + "ct10S5", chromhmm=True)
gen = msa.gen_putative(psearch_hg38 + ".sam", subset=[ct10S6])
msa.get_targets_stats(gen, hg38[0], psearch_out + "ct10S6", chromhmm=True) """
# ###########################################################################################
# Analyze mutation rates in NGS data using three different approaches.
mut_rate_file1 = lineageNGS_210908 + "mutation_rates_reg_dose_1.csv"
mut_rate_file2 = lineageNGS_210908 + "mutation_rates_reg_dose_2.csv"
mut_rate_file3 = lineageNGS_210908 + "mutation_rates_reg_dose_3.csv"
out_rate1 = open(mut_rate_file1, 'w')
out_rate2 = open(mut_rate_file2, 'w')
out_rate3 = open(mut_rate_file3, 'w')
Time_Points = [0, 0.5, 1, 2, 3, 5]
cter = 0
out_rate1.write("Time, T1, T2, T3, T4, T5, T6, T7, T8\n")
out_rate2.write("Time, T1, T2, T3, T4, T5, T6, T7, T8\n")
out_rate3.write("Time, T1, T2, T3, T4, T5, T6, T7, T8\n")
for file in "WT", "Reg-12h", "Reg-24h", "Reg-2d", "Reg-3d", "Reg-5d":
    infile = lineageNGS_210908 + file
    Mut1 = ltr.lineage_ngs(infile, hg38[1])
    Mut2 = ltr.lineage_ngs2(infile)
    Mut3 = ltr.lineage_ngs3(infile)
    out_rate1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (Time_Points[cter], Mut1[0], Mut1[1], Mut1[2], Mut1[3]
                                                  , Mut1[4], Mut1[5], Mut1[6], Mut1[7]))
    out_rate2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (Time_Points[cter], Mut2[0], Mut2[1], Mut2[2], Mut2[3]
                                                  , Mut2[4], Mut2[5], Mut2[6], Mut2[7]))
    out_rate3.write("%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (Time_Points[cter], Mut3[0], Mut3[1], Mut3[2], Mut3[3]
                                                   , Mut3[4], Mut3[5], Mut3[6], Mut3[7]))
    cter = cter + 1
out_rate1.close()
out_rate2.close()
out_rate3.close()