"""
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
ana_3 = ana + "3_ampngs/"
os.makedirs(ana_3) if not os.path.exists(ana_3) else None
ana_4 = ana + "4_ampngs2/"
os.makedirs(ana_4) if not os.path.exists(ana_4) else None
ana_5 = ana + "5_ampngs3/"
os.makedirs(ana_5) if not os.path.exists(ana_5) else None
ana_6 = ana + "6_lineage1/"
os.makedirs(ana_6) if not os.path.exists(ana_6) else None

""" Set file paths """
psearch_hg38_Alu = datadir + "Alu_ana_1_putative/1_protosearch/psearch_hg38_Alu"
psearch_dr11_DR1 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR1"
psearch_dr11_DR2 = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DR2"
psearch_dr11_DAN = datadir + "Alu_ana_1_putative/1_protosearch/psearch_dr11_DAN"
lin_210830 = datadir + "210830_lineageNGS/"     # amplicon NGS sequencing from 210830
lin_210929 = datadir + "210929_lineageNGS/"     # amplicon NGS sequencing from 210929
lin_211019 = datadir + "211019_lineageNGS/"     # amplicon NGS sequencing from 211019
lin_211028 = datadir + "211028_mutationNGS/"    # amplicon NGS sequencing from 211028
lin_211102 = datadir + "211102_lineageNGS/"     # amplicon NGS sequencing from 211102
lin_211119 = datadir + "211119_mutationNGS/"    # amplicon NGS sequencing from 211119

hg38_Alu_05_1 = "CAGGCGTGAGCTACTACGCC"
hg38_Alu_10_1 = "CCAGGCTGGAGTGCAGTGCT"
hg38_Alu_20_1 = "GCACTCCAGCCTGGGTTACA"

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


""" ############################################################################################ """
# Find PCR primers for target sites of hg38 Alu mgRNAs that target 05 sites
outfile = ana_2 + "lineage_hg38_Alu_05"
gen = msa.gen_putative(psearch_hg38_Alu + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath,
                       ct_values=[5], rad1=75, rad2=150, rad3=250)
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])

# Find PCR primers for target sites of hg38 Alu mgRNAs that target 10 sites (different parameters)
outfile = ana_2 + "lineage_hg38_Alu_10"
gen = msa.gen_putative(psearch_hg38_Alu + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath,
                       ct_values=[10], rad1=100, rad2=200, rad3=350)
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])

# Find PCR primers for target sites of hg38 Alu mgRNAs that target 20 sites
outfile = ana_2 + "lineage_hg38_Alu_20"
gen = msa.gen_putative(psearch_hg38_Alu + ".sam")
ltr.get_primers_nested(gen, outfile, hg38[0], genome_savepath,
                       ct_values=[20], rad1=75, rad2=150, rad3=250)
# Run wrapper function that align primer pairs to genome, parses SAM, obtain stats on best gRNAs
ltr.bowtie_parse_stats_wrapper(outfile, hg38[1])

# Retrieve the discovered PCR primers only for selected hg38 Alu mgRNAs
ip_05 = ana_2 + "lineage_hg38_Alu_05"
ip_10 = ana_2 + "lineage_hg38_Alu_10"
ip_20 = ana_2 + "lineage_hg38_Alu_20"
ltr.get_nested_primers_for_protospacer(ip_05 + "_inn_msa", ip_05 + "_out_msa",
                                       ana_2 + "hg38_Alu_05_1", hg38_Alu_05_1)
ltr.get_nested_primers_for_protospacer(ip_10 + "_inn_msa", ip_10 + "_out_msa",
                                       ana_2 + "hg38_Alu_10_1", hg38_Alu_10_1)
ltr.get_nested_primers_for_protospacer(ip_20 + "_inn_msa", ip_20 + "_out_msa",
                                       ana_2 + "hg38_Alu_20_1", hg38_Alu_20_1)


""" ############################################################################################ """
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_10_1])
msa.get_target_sequences(gen, ana_3 + "hg38_Alu_10_1", hg38[0], genome_savepath, win=500)
# align len read1, align len read2, read for mutation analysis, reverse complement?
al1, al2, ri, rc = 75, 75, 2, True

ltr.lineage_ngs_fq2sam(lin_210929 + "T0", hg38[1], ana_3 + "T0", readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T2-Ind", hg38[1], ana_3 + "T2-Ind", readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T2-NoInd", hg38[1], ana_3 + "T2-NoInd",  readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T6-Ind", hg38[1], ana_3 + "T6-Ind", readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T6-NoInd", hg38[1], ana_3 + "T6-NoInd", readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T10-Ind", hg38[1], ana_3 + "T10-Ind", readi=ri)
ltr.lineage_ngs_fq2sam(lin_210929 + "T10-NoInd", hg38[1], ana_3 + "T10-NoInd", readi=ri)

ltr.lineage_ngs_sam2dict(ana_3 + "T0", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T2-Ind", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T2-NoInd", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T6-Ind", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T6-NoInd", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T10-Ind", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_3 + "T10-NoInd", ana_3 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)

ltr.lineage_ngs_dict2csv(ana_3 + "T0", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T2-Ind", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T2-NoInd", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T6-Ind", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T6-NoInd", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T10-Ind", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_3 + "T10-NoInd", hg38_Alu_10_1)

list1 = [ana_3 + "T0", ana_3 + "T2-Ind", ana_3 + "T6-Ind", ana_3 + "T10-Ind"]
list2 = [ana_3 + "T0", ana_3 + "T2-NoInd", ana_3 + "T6-NoInd", ana_3 + "T10-NoInd"]
ltr.lineage_ngs_np2sum(list1, "1")
ltr.lineage_ngs_np2sum(list2, "2")
ltr.lineage_ngs_aggregate(list1, "1", ana_3 + "summary_Ind.csv")
ltr.lineage_ngs_aggregate(list2, "2", ana_3 + "summary_noInd.csv")


""" ############################################################################################ """
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_05_1])
msa.get_target_sequences(gen, ana_4 + "hg38_Alu_05_1", hg38[0], genome_savepath, win=500)
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_10_1])
msa.get_target_sequences(gen, ana_4 + "hg38_Alu_10_1", hg38[0], genome_savepath, win=500)
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_20_1])
msa.get_target_sequences(gen, ana_4 + "hg38_Alu_20_1", hg38[0], genome_savepath, win=500)
# align len read1, align len read2, read for mutation analysis, reverse complement?
al1, al2, ri, rc = 75, 50, 1, False

ltr.lineage_ngs_fq2sam(lin_211019 + "ct5-Set1", hg38[1], ana_4 + "ct5-Set1", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct10-Set1", hg38[1], ana_4 + "ct10-Set1", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct10-Set2", hg38[1], ana_4 + "ct10-Set2", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct10-Set3", hg38[1], ana_4 + "ct10-Set3", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct10-Set4", hg38[1], ana_4 + "ct10-Set4", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct10-Set5", hg38[1], ana_4 + "ct10-Set5", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct20-Set1", hg38[1], ana_4 + "ct20-Set1", al1=al1, al2=al2)
ltr.lineage_ngs_fq2sam(lin_211019 + "ct20-Set2", hg38[1], ana_4 + "ct20-Set2", al1=al1, al2=al2)

ltr.lineage_ngs_sam2dict(ana_4 + "ct5-Set1", ana_4 + "hg38_Alu_05_1.csv", hg38_Alu_05_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct10-Set1", ana_4 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct10-Set2", ana_4 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct10-Set3", ana_4 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct10-Set4", ana_4 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct10-Set5", ana_4 + "hg38_Alu_10_1.csv", hg38_Alu_10_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct20-Set1", ana_4 + "hg38_Alu_20_1.csv", hg38_Alu_20_1)
ltr.lineage_ngs_sam2dict(ana_4 + "ct20-Set2", ana_4 + "hg38_Alu_20_1.csv", hg38_Alu_20_1)

ltr.lineage_ngs_dict2csv(ana_4 + "ct5-Set1", hg38_Alu_05_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct10-Set1", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct10-Set2", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct10-Set3", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct10-Set4", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct10-Set5", hg38_Alu_10_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct20-Set1", hg38_Alu_20_1)
ltr.lineage_ngs_dict2csv(ana_4 + "ct20-Set2", hg38_Alu_20_1)

list1 = [ana_4 + "ct5-Set1"]
list2 = [ana_4 + "ct10-Set1", ana_4 + "ct10-Set2", ana_4 + "ct10-Set3",
         ana_4 + "ct10-Set4", ana_4 + "ct10-Set5"]
list3 = [ana_4 + "ct20-Set1", ana_4 + "ct20-Set2"]
ltr.lineage_ngs_np2sum(list1, "1")
ltr.lineage_ngs_np2sum(list2, "2")
ltr.lineage_ngs_np2sum(list3, "3")
ltr.lineage_ngs_aggregate(list1, "1", ana_4 + "summary_ct05.csv")
ltr.lineage_ngs_aggregate(list2, "2", ana_4 + "summary_ct10.csv")
ltr.lineage_ngs_aggregate(list3, "3", ana_4 + "summary_ct20.csv")


""" ############################################################################################ """
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_05_1])
msa.get_target_sequences(gen, ana_5 + "hg38_Alu_05_1", hg38[0], genome_savepath, win=500)
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_10_1])
msa.get_target_sequences(gen, ana_5 + "hg38_Alu_10_1", hg38[0], genome_savepath, win=500)
# align len read1, align len read2, read for mutation analysis, reverse complement?
al1, al2, ri, rc = 75, 50, 1, False

ct05_R1 = ["ct5-R1-Ind-T2", "ct5-R1-Ind-T6", "ct5-R1-Ind-T10",
           "ct5-R1-NI-T2", "ct5-R1-NI-T6", "ct5-R1-NI-T10", "ct5-R1-T0"]
ct05_R2 = ["ct5-R2-Ind-T2", "ct5-R2-Ind-T6", "ct5-R2-Ind-T10",
           "ct5-R2-NI-T2", "ct5-R2-NI-T6", "ct5-R2-NI-T10", "ct5-R2-T0"]
ct10_R1 = ["ct10-R1-Ind-T2", "ct10-R1-Ind-T6", "ct10-R1-Ind-T10",
           "ct10-R1-NI-T2", "ct10-R1-NI-T6", "ct10-R1-NI-T10", "ct10-R1-T0"]
ct10_R2 = ["ct10-R2-Ind-T2", "ct10-R2-Ind-T6", "ct10-R2-Ind-T10",
           "ct10-R2-NI-T2", "ct10-R2-NI-T6", "ct10-R2-NI-T10", "ct10-R2-T0"]

[ltr.lineage_ngs_fq2sam(lin_211028 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct05_R1]
[ltr.lineage_ngs_fq2sam(lin_211028 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct05_R2]
[ltr.lineage_ngs_fq2sam(lin_211028 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct10_R1]
[ltr.lineage_ngs_fq2sam(lin_211028 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct10_R2]

[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_05_1.csv", hg38_Alu_05_1) for x in ct05_R1]
[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_05_1.csv", hg38_Alu_05_1) for x in ct05_R2]
[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_10_1.csv", hg38_Alu_10_1) for x in ct10_R1]
[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_10_1.csv", hg38_Alu_10_1) for x in ct10_R2]

[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_05_1) for x in ct05_R1]
[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_05_1) for x in ct05_R2]
[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_10_1) for x in ct10_R1]
[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_10_1) for x in ct10_R2]

ct05_R1 = [ana_5 + "ct5-R1-T0", ana_5 + "ct5-R1-Ind-T2", ana_5 + "ct5-R1-Ind-T6",
           ana_5 + "ct5-R1-Ind-T10",
           ana_5 + "ct5-R1-NI-T2", ana_5 + "ct5-R1-NI-T6", ana_5 + "ct5-R1-NI-T10"]
ct05_R2 = [ana_5 + "ct5-R2-T0", ana_5 + "ct5-R2-Ind-T2", ana_5 + "ct5-R2-Ind-T6",
           ana_5 + "ct5-R2-Ind-T10",
           ana_5 + "ct5-R2-NI-T2", ana_5 + "ct5-R2-NI-T6", ana_5 + "ct5-R2-NI-T10"]
ct10_R1 = [ana_5 + "ct10-R1-T0", ana_5 + "ct10-R1-Ind-T2", ana_5 + "ct10-R1-Ind-T6",
           ana_5 + "ct10-R1-Ind-T10",
           ana_5 + "ct10-R1-NI-T2", ana_5 + "ct10-R1-NI-T6", ana_5 + "ct10-R1-NI-T10"]
ct10_R2 = [ana_5 + "ct10-R2-T0", ana_5 + "ct10-R2-Ind-T2", ana_5 + "ct10-R2-Ind-T6",
           ana_5 + "ct10-R2-Ind-T10",
           ana_5 + "ct10-R2-NI-T2", ana_5 + "ct10-R2-NI-T6", ana_5 + "ct10-R2-NI-T10"]
ltr.lineage_ngs_np2sum(ct05_R1, "ct05_R1")
ltr.lineage_ngs_np2sum(ct05_R2, "ct05_R2")
ltr.lineage_ngs_np2sum(ct10_R1, "ct10_R1")
ltr.lineage_ngs_np2sum(ct10_R2, "ct10_R2")
ltr.lineage_ngs_aggregate(ct05_R1, "ct05_R1", ana_5 + "summary_ct05_R1.csv")
ltr.lineage_ngs_aggregate(ct05_R2, "ct05_R2", ana_5 + "summary_ct05_R2.csv")
ltr.lineage_ngs_aggregate(ct10_R1, "ct10_R1", ana_5 + "summary_ct10_R1.csv")
ltr.lineage_ngs_aggregate(ct10_R2, "ct10_R2", ana_5 + "summary_ct10_R2.csv")


""" ############################################################################################ """
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_20_1])
msa.get_target_sequences(gen, ana_5 + "hg38_Alu_20_1", hg38[0], genome_savepath, win=500)
# align len read1, align len read2, read for mutation analysis, reverse complement?
al1, al2, ri, rc = 75, 50, 1, False

ct20_R1 = ["ct20-T0-R1-Set3", "ct20-T0-R1-Set4", "ct20-T0-R1-Set5",
           "ct20-T2-R1-Set3", "ct20-T2-R1-Set4", "ct20-T2-R1-Set5",
           "ct20-T6-R1-Set3", "ct20-T6-R1-Set4", "ct20-T6-R1-Set5",
           "ct20-T10-R1-Set3", "ct20-T10-R1-Set4", "ct20-T10-R1-Set5"]
ct20_R2 = ["ct20-T0-R2-Set3", "ct20-T0-R2-Set4", "ct20-T0-R2-Set5",
           "ct20-T2-R2-Set3", "ct20-T2-R2-Set4", "ct20-T2-R2-Set5",
           "ct20-T6-R2-Set3", "ct20-T6-R2-Set4", "ct20-T6-R2-Set5",
           "ct20-T10-R2-Set3", "ct20-T10-R2-Set4", "ct20-T10-R2-Set5"]

[ltr.lineage_ngs_fq2sam(lin_211119 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct20_R1]
[ltr.lineage_ngs_fq2sam(lin_211119 + x, hg38[1], ana_5 + x, al1=al1, al2=al2) for x in ct20_R2]

[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_20_1.csv", hg38_Alu_20_1) for x in ct20_R1]
[ltr.lineage_ngs_sam2dict(ana_5 + x, ana_5 + "hg38_Alu_20_1.csv", hg38_Alu_20_1) for x in ct20_R2]

[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_20_1) for x in ct20_R1]
[ltr.lineage_ngs_dict2csv(ana_5 + x, hg38_Alu_20_1) for x in ct20_R2]

ct20_S3_R1 = [ana_5 + "ct20-T0-R1-Set3", ana_5 + "ct20-T2-R1-Set3",
              ana_5 + "ct20-T6-R1-Set3", ana_5 + "ct20-T10-R1-Set3"]
ct20_S4_R1 = [ana_5 + "ct20-T0-R1-Set4", ana_5 + "ct20-T2-R1-Set4",
              ana_5 + "ct20-T6-R1-Set4", ana_5 + "ct20-T10-R1-Set4"]
ct20_S5_R1 = [ana_5 + "ct20-T0-R1-Set5", ana_5 + "ct20-T2-R1-Set5",
              ana_5 + "ct20-T6-R1-Set5", ana_5 + "ct20-T10-R1-Set5"]
ct20_S3_R2 = [ana_5 + "ct20-T0-R2-Set3", ana_5 + "ct20-T2-R2-Set3",
              ana_5 + "ct20-T6-R2-Set3", ana_5 + "ct20-T10-R2-Set3"]
ct20_S4_R2 = [ana_5 + "ct20-T0-R2-Set4", ana_5 + "ct20-T2-R2-Set4",
              ana_5 + "ct20-T6-R2-Set4", ana_5 + "ct20-T10-R2-Set4"]
ct20_S5_R2 = [ana_5 + "ct20-T0-R2-Set5", ana_5 + "ct20-T2-R2-Set5",
              ana_5 + "ct20-T6-R2-Set5", ana_5 + "ct20-T10-R2-Set5"]
ltr.lineage_ngs_np2sum(ct20_S3_R1, "ct20_S3_R1")
ltr.lineage_ngs_np2sum(ct20_S4_R1, "ct20_S4_R1")
ltr.lineage_ngs_np2sum(ct20_S5_R1, "ct20_S5_R1")
ltr.lineage_ngs_np2sum(ct20_S3_R2, "ct20_S3_R2")
ltr.lineage_ngs_np2sum(ct20_S4_R2, "ct20_S4_R2")
ltr.lineage_ngs_np2sum(ct20_S5_R2, "ct20_S5_R2")
ltr.lineage_ngs_aggregate(ct20_S3_R1, "ct20_S3_R1", ana_5 + "summary_ct20_S3_R1.csv")
ltr.lineage_ngs_aggregate(ct20_S4_R1, "ct20_S4_R1", ana_5 + "summary_ct20_S4_R1.csv")
ltr.lineage_ngs_aggregate(ct20_S5_R1, "ct20_S5_R1", ana_5 + "summary_ct20_S5_R1.csv")
ltr.lineage_ngs_aggregate(ct20_S3_R2, "ct20_S3_R2", ana_5 + "summary_ct20_S3_R2.csv")
ltr.lineage_ngs_aggregate(ct20_S4_R2, "ct20_S4_R2", ana_5 + "summary_ct20_S4_R2.csv")
ltr.lineage_ngs_aggregate(ct20_S5_R2, "ct20_S5_R2", ana_5 + "summary_ct20_S5_R2.csv")


""" ############################################################################################ """
gen = msa.gen_putative(psearch_hg38_Alu + ".sam", subset=[hg38_Alu_10_1])
msa.get_target_sequences(gen, ana_6 + "hg38_Alu_10_1", hg38[0], genome_savepath, win=500)
# align len read1, align len read2, read for mutation analysis, reverse complement?
al1, al2, ri, rc = 75, 50, 1, False

ct10_100c = ["100c-111", "100c-112", "100c-121", "100c-122",
             "100c-211", "100c-212", "100c-221", "100c-222"]
ct10_250c = ["250c-111", "250c-112", "250c-121", "250c-122",
             "250c-211", "250c-212", "250c-221", "250c-222"]

[ltr.lineage_ngs_fq2sam(lin_211102 + x, hg38[1], ana_6 + x, al1=al1, al2=al2) for x in ct10_100c]
[ltr.lineage_ngs_fq2sam(lin_211102 + x, hg38[1], ana_6 + x, al1=al1, al2=al2) for x in ct10_250c]

[ltr.lineage_ngs_sam2dict(ana_6 + x, ana_6 + "hg38_Alu_10_1.csv", hg38_Alu_10_1) for x in ct10_100c]
[ltr.lineage_ngs_sam2dict(ana_6 + x, ana_6 + "hg38_Alu_10_1.csv", hg38_Alu_10_1) for x in ct10_250c]

[ltr.lineage_ngs_dict2csv(ana_6 + x, hg38_Alu_10_1) for x in ct10_100c]
[ltr.lineage_ngs_dict2csv(ana_6 + x, hg38_Alu_10_1) for x in ct10_250c]

ct10_100c = [ana_6 + "100c-111", ana_6 + "100c-112", ana_6 + "100c-121", ana_6 + "100c-122",
             ana_6 + "100c-211", ana_6 + "100c-212", ana_6 + "100c-221", ana_6 + "100c-222"]
ct10_250c = [ana_6 + "250c-111", ana_6 + "250c-112", ana_6 + "250c-121", ana_6 + "250c-122",
             ana_6 + "250c-211", ana_6 + "250c-212", ana_6 + "250c-221", ana_6 + "250c-222"]
ltr.lineage_ngs_np2sum(ct10_100c, "ct10_100c")
ltr.lineage_ngs_np2sum(ct10_250c, "ct10_250c")
ltr.lineage_ngs_aggregate(ct10_100c, "ct10_100c", ana_6 + "summary_ct10_100c.csv")
ltr.lineage_ngs_aggregate(ct10_250c, "ct10_250c", ana_6 + "summary_ct10_250c.csv")

ltr.lineage_ngs_distance(ct10_100c, "ct10_100c", ana_6 + "summary_ct10_100c.csv")
ltr.lineage_ngs_distance(ct10_250c, "ct10_250c", ana_6 + "summary_ct10_250c.csv")
