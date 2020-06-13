
import src.chipseq as c
import src.mtss as mtss


""" File paths """
desktop = "/Users/rogerzou/Desktop/"
hg38 = "/Users/rogerzou/bioinformatics/hg38_bowtie2/hg38.fa"
labhome = "/Volumes/Lab-Home/rzou4/NGS_data/4_damage/"
ana = labhome + "Alu_analysis/"
mreGGin = labhome + "200206_chipseq/16_AluGG-MRE11_final.bam"
casGGin = labhome + "200206_chipseq/15_AluGG-Cas9_final.bam"
mreTAin = labhome + "200316_chipseq/AluTA-mre11-rep1_rmdup.bam"
casTAin = labhome + "200316_chipseq/AluTA-cas9-rep1_rmdup.bam"
mreCTin = labhome + "200316_chipseq/AluCT-mre11-rep1_rmdup.bam"
casCTin = labhome + "200316_chipseq/AluCT-cas9-rep1_rmdup.bam"
mreGGin_nD = labhome + "200316_chipseq/AluGG-mre11-noD-rep1_rmdup.bam"
mreGGin_PK = labhome + "200316_chipseq/AluGG-mre11-PKi-rep1_rmdup.bam"

""" Sequences """
Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTG" \
      "GCCAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGC" \
      "TGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAG" \
      "ACTCCGTCTC"
AluGG = "CCTGTAGTCCCAGCTACTGG"
AluCT = "CCTGTAGTCCCAGCTACTCT"
AluTA = "CCTGTAGTCCCAGCTACTTA"


""" Starting from Alu, find all sequences with at most 3 mismatches within 9 bases from PAM.
    Determine all putative on-target genomic sites, calculate distances between adjacent targets """
# mtss.get_targets_fasta(desktop + 'basemut.fa', Alu, 9)
# mtss.get_targets_bowtie2(desktop + 'basemut.fa', desktop + 'basemut.sam', hg38)
# mtss.get_targets_stats(desktop + 'basemut.sam', desktop + 'basemut')
# mtss.get_targets_dist(desktop + 'basemut_align.csv', desktop + 'basemut')

""" Determine paired-end read subsets for all putative on-target sites """
alnpath = desktop + 'basemut_align.csv'
# mtss.read_subsets(mtss.target_gen(alnpath, 1250, AluGG), mreGGin, ana + "GG-TARGET_mre11_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 1250, AluCT), mreCTin, ana + "CT-TARGET_mre11_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 1250, AluTA), mreTAin, ana + "TA-TARGET_mre11_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 750, AluGG), casGGin, ana + "GG-TARGET_cas9_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 750, AluCT), casCTin, ana + "CT-TARGET_cas9_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 750, AluTA), casTAin, ana + "TA-TARGET_cas9_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 1250, AluGG), mreGGin_nD, ana + "GG-TARGET-noD_mre11_rs")
# mtss.read_subsets(mtss.target_gen(alnpath, 1250, AluGG), mreGGin_PK, ana + "GG-TARGET-PKi_mre11_rs")
# GG_mre = mtss.load_nparray(ana + "GG-TARGET_mre11_rs.csv")
# CT_mre = mtss.load_nparray(ana + "CT-TARGET_mre11_rs.csv")
# TA_mre = mtss.load_nparray(ana + "TA-TARGET_mre11_rs.csv")
# GG_cas = mtss.load_nparray(ana + "GG-TARGET_cas9_rs.csv")
# CT_cas = mtss.load_nparray(ana + "CT-TARGET_cas9_rs.csv")
# TA_cas = mtss.load_nparray(ana + "TA-TARGET_cas9_rs.csv")
# mtss.mergerows([GG_mre, CT_mre, TA_mre], ana + "TARGET-mre-merged_rs.csv")
# mtss.mergerows([GG_cas, CT_cas, TA_cas], ana + "TARGET-cas-merged_rs.csv")

""" Figure S4 """
mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluTA), ana + "TA-TARGET_mre11_rs_abut.bam", ana + "TA-TARGET_mre11_rs_abut")
mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluTA), ana + "TA-TARGET_cas9_rs_abut.bam", ana + "TA-TARGET_cas9_rs_abut")

""" Determine peak profiles centered at the cut site for all putative on-target sites """
# alnpath = desktop + 'basemut_align.csv'
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluGG), mreGGin, ana + "GG-TARGET_mre11")
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluCT), mreCTin, ana + "CT-TARGET_mre11")
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluTA), mreTAin, ana + "TA-TARGET_mre11")
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluGG), casGGin, ana + "GG-TARGET_cas9")
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluCT), casCTin, ana + "CT-TARGET_cas9")
# mtss.peak_profile(mtss.target_gen(alnpath, 1500, AluTA), casTAin, ana + "TA-TARGET_cas9")
# GG_mre = mtss.load_nparray(ana + "GG-TARGET_mre11_bpeaks.csv")
# CT_mre = mtss.load_nparray(ana + "CT-TARGET_mre11_bpeaks.csv")
# TA_mre = mtss.load_nparray(ana + "TA-TARGET_mre11_bpeaks.csv")
# GG_cas = mtss.load_nparray(ana + "GG-TARGET_cas9_bpeaks.csv")
# CT_cas = mtss.load_nparray(ana + "CT-TARGET_cas9_bpeaks.csv")
# TA_cas = mtss.load_nparray(ana + "TA-TARGET_cas9_bpeaks.csv")
# mtss.mergerows([GG_mre, CT_mre, TA_mre], ana + "TARGET-mre-merged_bp.csv")
# mtss.mergerows([GG_cas, CT_cas, TA_cas], ana + "TARGET-cas-merged_bp.csv")

""" Determine multiple sequence alignments from paired-end reads around putative on-target sites """
# alnpath = desktop + 'basemut_align.csv'
# mtss.find_msa(mtss.target_gen(alnpath, 750, AluGG), casGGin, ana + "GG-TARGET_cas9_msa", hg38)
# mtss.find_msa(mtss.target_gen(alnpath, 750, AluCT), casCTin, ana + "CT-TARGET_cas9_msa", hg38)
# mtss.find_msa(mtss.target_gen(alnpath, 750, AluTA), casTAin, ana + "TA-TARGET_cas9_msa", hg38)
# mtss.find_msa(mtss.target_gen(alnpath, 1250, AluGG), mreGGin, ana + "GG-TARGET_mre11_msa", hg38)
# mtss.find_msa(mtss.target_gen(alnpath, 1250, AluCT), mreCTin, ana + "CT-TARGET_mre11_msa", hg38)
# mtss.find_msa(mtss.target_gen(alnpath, 1250, AluTA), mreTAin, ana + "TA-TARGET_mre11_msa", hg38)
