import pysam
import numpy as np
from Bio import Entrez, SeqIO
from . import chipseq as c
from . import msa as msa
from . import mtss as m
import pickle
import subprocess as sp
import statistics
import re
import os
import csv
import matplotlib as mpl
from matplotlib import pyplot as plt

GENOME = ""
GENOME_LIST = ['hg38', 'hg19', 'mm10']
genome_size, genome_id, genome_seq = None, None, None

CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX', 'chrY']
CHROMHMM = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts',
            '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']
fullflags_paired = ['83', '163', '99', '147', '339', '419', '355', '403', '77', '141']
initflags_paired = ['83', '99', '339', '355']
primflags_paired = ['83', '99']
fullflags_single = ['0', '16', '256', '272', '4']
initflags_single = ['0', '16', '256', '272']
primflags_single = ['0', '16']

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def get_primers_nested(gen, outfile, genome_str, savepath, ct_values, win=350, npr = 3,
                       rng1_min=200, rng1_max = 500, rng2_min = 300, rng2_max = 400,
                       amp_size_min=100, amp_size_max=350):
    """
    :param gen: bowtie2 alignments for each putative protospacer sequence of format:
        - boolean: whether it is a new putative protospacer sequence
        - string: sequence of putative protospacer sequence
        - string: '+' or '-' that correspond to the sense of current protospacer sequence alignment
        - string: chromosome string (e.g. 'chr1') of current protospacer sequence alignment
        - int: coordinate of current protospacer sequence alignment
        - int: total number of alignments for current protospacer sequence
    :param outfile: string path to output file (extension omitted). Four fasta files are generated:
        - outfile + "_out_1/2.fa" for the outer PCR primers (PCR-1)
        - outfile + "_in_1/2.fa" for the inner PCR primers (PCR-2)
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :param savepath: save path to genome sequence that will be downloaded from NCBI if not already
    :param ct_values: number of targets selected
    :param npr: number of pairs of primers that will be sought in each of the nested PCR,
        e.g. npr=3 will search for 3 for/rev primer pairs in PCR-1 and 3 for/rev primer pairs in PCR-2.
    :param win: radius of base pairs to read, centered at each target site
    :params (rng1_min, rng1_max): range to find primers for 1st nested PCR. Primers are restricted
        to be in [0:rng1_min] (left primer) and [rng1_max:2*win] (right primer)
    :params (rng2_min, rng2_max): range to find primers for 2nd nested PCR. Primers are restricted
        to be in [rng1_min:rng2_min] (left primer) and [rng2_max:rng1_max] (right primer)
    :params (amp_size_min, amp_size_max): range of size for the final amplicon after 2nd nested PCR
    """
    amp_size = str(amp_size_min) + "-" + str(amp_size_max)
    proto_i, align_i = -1, -1
    genome_initialized(savepath, genome_str)
    cter = 0
    with open(outfile + "_in_1.fa", 'w') as f1in, open(outfile + "_in_2.fa", 'w') as f2in, open(outfile + "_out_1.fa", 'w') as f1out, open(outfile + "_out_2.fa", 'w') as f2out:
        for new_i, seq_i, sen_i, chr_i, coo_i, tct_i in gen:
            if ((tct_i in ct_values) and (chr_i in CHR)):
                proto_i = proto_i + 1 if new_i else proto_i
                align_i = 0 if new_i else align_i + 1
                seq = get_target_region(chr_i, coo_i, sen_i, win)
                r1out, r2out = get_primer3_primers(seq, num_primers=npr, rng_min=rng1_min, rng_max=rng1_max, prod_size=None)
                r1in, r2in = get_primer3_primers(seq[rng1_min:rng1_max], num_primers=npr, rng_min=rng2_min-rng1_min,
                                                 rng_max=rng2_max-rng1_min, prod_size=amp_size)
                for i in range(npr):
                    if ((r1out[i] != 0) & (r2out[i] != 0)):
                        f1out.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                                     align_i, i, r1out[i]))
                        f2out.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                                     align_i, i, r2out[i]))
                    if ((r1in[i] != 0) & (r2in[i] != 0)):
                        f1in.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                                     align_i, i, r1in[i]))
                        f2in.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                                     align_i, i, r2in[i]))
                if cter % 10000 == 0:
                    print("get_primers_nested(): processed %i samples" % cter)
                cter += 1


def get_primer3_primers(seq, num_primers, rng_min, rng_max, prod_size=None):
    """ Finds primers using primer3
    :param seq: template sequence to find primers
    :param num_primers: number of primers to seek for
    :params (rng_min, rng_max): range of sequence that is considered target region for PCR.
    :param prod_size: desired amplicon size.
    """
    inptprmr3, outprmr3 = "primer3_input.dat", "primer3_output.dat"
    with open(inptprmr3, 'w') as prmr3_seqs:
        prmr3_seqs.write("%s%s\n" % ("SEQUENCE_ID=Proto-", 0))
        prmr3_seqs.write("%s%s\n" % ("SEQUENCE_TEMPLATE=", seq))
        _get_primer3_primers_helper(prmr3_seqs, num_primers, rng_min, rng_max, prod_size)
    cmnd_prmr3 = "primer3_core " + inptprmr3 + "> " + outprmr3
    os.system(cmnd_prmr3)
    primers_list, primers_right, primers_left = [], [0] * num_primers, [0] * num_primers
    cter_l, cter_r = int(0), int(0)
    # Define strings to find the lines that contain the primers information:
    find_l, find_r = "PRIMER_LEFT_" + str(cter_l) + "_SEQUENCE", "PRIMER_RIGHT_" + str(cter_r) + "_SEQUENCE"
    with open(outprmr3, 'r') as prmrs_file:
        for line in prmrs_file:
            row = line.split('=')  # split line by the equal sign
            if row[0] == find_l:
                primers_left[cter_l] = row[1][:-1]
                cter_l += 1
                # update primer left line to find next left primer
                find_l = "PRIMER_LEFT_" + str(cter_l) + "_SEQUENCE"
            if row[0] == find_r:
                primers_right[cter_r] = row[1][:-1]
                cter_r += 1
                # update primer right line to find next right primer
                find_r = "PRIMER_RIGHT_" + str(cter_r) + "_SEQUENCE"
    # clean up!
    cmnd_clean = "rm " + inptprmr3 + " " + outprmr3
    os.system(cmnd_clean)
    return primers_left, primers_right

def _get_primer3_primers_helper(prim3file, num_prims, rng_min, rng_max, prod_size=None):
    """
    Write primer3 input file settings to "prim3file". Settings that can be changed are:
        :param num_primers: number of primers to seek for
        :params (rng_min, rng_max): range of sequence that is considered target region for PCR.
        :param prod_size: desired amplicon size.
    """
    # Write sequence target (primers will bind outside of this)
    prim3_trgt = "SEQUENCE_TARGET=" + str(rng_min) + "," + str(rng_max-rng_min) + "\n"
    prim3file.write(prim3_trgt)
    prim3file.write("PRIMER_TASK=generic\n") # generic: default option for primer3
    prim3file.write("PRIMER_PICK_LEFT_PRIMER=1\n")
    prim3_return = "PRIMER_NUM_RETURN=" + str(num_prims) + "\n"
    prim3file.write(prim3_return) # choose how many primer sets we want
    prim3file.write("PRIMER_PICK_INTERNAL_OLIGO=0\n") # we don't need internal oligos
    prim3file.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
    prim3file.write("PRIMER_OPT_SIZE=20\n")
    prim3file.write("PRIMER_MIN_SIZE=18\n")
    prim3file.write("PRIMER_MAX_SIZE=22\n")
    prim3file.write("PRIMER_MAX_GC=65\n")
    prim3file.write("PRIMER_MIN_GC=35\n")
    prim3file.write("PRIMER_MIN_TM=57\n")
    prim3file.write("PRIMER_MAX_TM=63\n")
    if prod_size is not None:
        prim3_prod_size = "PRIMER_PRODUCT_SIZE_RANGE=" + prod_size + "\n"
        prim3file.write(prim3_prod_size)
    else:
        prim3file.write("PRIMER_PRODUCT_SIZE_RANGE=150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000\n")
    prim3file.write("PRIMER_EXPLAIN_FLAG=1\n")
    prim3file.write("=\n")


def bowtie2_msa_primers(curfile, genome_path, k_count=10):
    """ Run bowtie2 to align paired-end reads in '-k' mode, followed by samtools to sort and index
        the alignment.

    :param curfile: path to base file name/path.
            The PE read files are of format: curfile + '_1.fa' | curfile + '_1.fa'.
            The output initial SAM file is of format: curfile + '_msa.sam'
            The output (unsorted) BAM file is of format: curfile + '_msa.bam'
            The output sorted BAM file is of format: curfile + '_msa_sorted.bam'
            The output BAM index file is of format: curfile + '_msa_sorted.bam.bai'
    :param genome_path: path to genome for bowtie2 to use
    :param k_count: number of distinct, valid alignments for each PE read
    """
    sp.run(['bowtie2', '-f', '-p', '8', '--local', '--score-min', 'L,0,1.5', '-k', str(k_count), '-X', '1000', '--no-mixed', '--no-discordant',
            '-L', '18', '-N', '1', '-x', genome_path[:-3], '-1', curfile + '_1.fa', '-2', curfile + '_2.fa', '-S', curfile + '_msa.sam'])
    sp.run(['samtools', 'view', '-h', '-S', '-b', '-o', curfile + '_msa.bam', curfile + '_msa.sam'])
    sp.run(['samtools', 'sort', '-o', curfile + '_msa_sorted.bam', curfile + '_msa.bam'])
    sp.run(['samtools', 'index', curfile + '_msa_sorted.bam'])


def parse_sam_primers(outfile, label):
    """ Parse the SAM file output from bowtie2 alignment of paired-end ChIP-seq reads to record the
        information about each PE read, number of alignments, uniqueness, and alignment score.
    :param outfile: path to output csv file (extension omitted)

    Outputs csv file with following columns:
     0. original protospacer sequence
     1-2. chromosome and coordinate of cut site
     3. total number of cut sites for specific protospacer sequence
     4. index of protospacer sequence
     5. index of cut site for a particular protospacer sequence
     6. index of PE ChIP-seq read for a particular cut site
     7. '1' if primary alignment of PE read is to the expected cut site, and there is only one
            alignment or multiple alignments but the primary alignment has the highest score
    """
    global fullflags_paired, initflags_paired, primflags_paired
    """ Parse bowtie2 samfile output, determine all primary/secondary alignments """
    csv_out = open(outfile + '.csv', 'w')
    outrow, outstr, scores = None, "", []
    sam_it = open(outfile + '.sam', 'r')
    for read in sam_it:         # read every alignment of SAM file
        if read[0] == '@':      # skip SAM header lines
            continue
        stp = read.strip()
        row = stp.split('\t')           # read each bowtie2 alignment in SAM format
        n_str = row[0].split('_')       # split the information of each aligned PE read
        seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i = n_str[:7]
        # first in pair, properly mapped - one (first) alignment of a new PE-read.
        if row[1] in initflags_paired:
            # get quality scores and location of alignment
            AS = int(re.findall(r'AS:i:[-0-9]+', stp)[0].split(':')[2])
            YS = int(re.findall(r'YS:i:[-0-9]+', stp)[0].split(':')[2])
            sco_c = (AS + YS) / 2.0                         # quality score of alignment
            chr_c, coo_c = row[2], int(row[3])              # chromosome and position of alignment
            str_c = "%s_%i_%0.3f" % (chr_c, coo_c, sco_c)   # summary string of alignment
            # first in pair, primary alignment - the alignment of a new PE-read designated 'primary'
            if row[1] in primflags_paired:
                if outrow is not None:  # save previous alignment set
                    len_scores = len(scores)
                    outrow[-1] = '1' if outrow[-1] and (
                            len_scores == 1 or scores[0] > max(scores[1:])) else '0'
                    outrow.append(str(len_scores))
                    outrow.append(outstr)
                    csv_out.write(','.join(outrow) + '\n')
                # determine if primary alignment matches previously assigned location
                int_c = int(coo_i)
                boo = True if chr_c == chr_i and int_c - 2e3 <= coo_c <= int_c + 2e3 else False
                outrow = [seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i, label, boo]
                outstr, scores = "", []
            outstr += str_c + "|"
            scores.append(sco_c)
        elif row[1] not in fullflags_paired:
            raise ValueError("parse_msa_sam_paired(): unexpected SAM flag: %s" % row[1])
    if outrow is not None:                  # save final alignment set
        len_scores = len(scores)
        outrow[-1] = '1' if outrow[-1] and (len_scores == 1 or scores[0] > max(scores[1:])) else '0'
        outrow.append(str(len_scores))
        outrow.append(outstr)
        csv_out.write(','.join(outrow) + '\n')
    sam_it.close()
    csv_out.close()


def get_stats_primers(outfile, min_prmrs=2, max_num_tgts=100):
    """ Given the output of 'parse_sam_primers()' determine, for each gRNA:
        (1) how many targets there are for that gRNA
        (2) how many targets passed the primer3 filter, with at least min_prmrs number of primers
        (3) how many targets passed the primer3 and bowtie2 filters, with at least min_prmrs number of primers
            when applying the bowtie2 filter, primers with more than one alignment are discarded,
            even if one of the alignments has much better score than the others
    :param outfile: path of output of 'parse_sam_primers_*()', extension omitted
    :param min_prmrs: minimum number of primer sets found in order to consider a target amplifiable
    :param max_num_prms: this number has to be larger than the maximum number of targets per gRNA

    Outputs csv file with following columns:
     0. original gRNA sequence
     1. total number of cut sites for gRNA sequence
     2. number of cut sites with > min_prmrs primers that pass primer3 filter
     3. number of cut sites with > min_prmrs primers that pass primer3 + bowtie2 (uniqueness mapping) filter
    """
    msa_in = m.load_nparray(outfile + "_in_msa.csv")
    msa_out = m.load_nparray(outfile + "_out_msa.csv")
    msa_all = np.concatenate((msa_in, msa_out), axis=0)
    msa_sorted = msa_all[np.lexsort((msa_all[:, 6].astype(int), msa_all[:, 5].astype(int),
                                 msa_all[:, 4].astype(int), -msa_all[:, 3].astype(int)))]
    """
    For verbose purposes:
    csv_merge = open(outfile + '_merge_msa.csv', 'w')
    for i in msa_sorted:
        csv_merge.write(','.join(i) + '\n')
    """
    csv_out = open(outfile + '_stats.csv', 'w')
    proto_prev, tgt_prev, boo_prev, num_prev, prim_prev, seq_prev = None, None, None, None, None, None
    write_row = None
    sumrow, numct = [], []
    num_prms_prim3, num_prms_final = [0]*max_num_tgts, [0]*max_num_tgts
    for msa_i in msa_sorted:
        seq_i, chr_i, coo_i, tnt_i = msa_i[:4]
        proto_i, boo_i, num_i = int(msa_i[4]), msa_i[8], int(msa_i[9])
        tgt_i, prim_i = int(msa_i[5]), int(msa_i[6])
        if proto_i != proto_prev:
            if len(numct) > 0:
                count_prim3 = len([i for i in num_prms_prim3 if i >= min_prmrs])
                count_final = len([i for i in num_prms_final if i >= min_prmrs])
                sumrow[2] = count_prim3
                sumrow[3] = count_final
                if write_row is not None:
                    write_row = np.append(write_row, [[sumrow[0], sumrow[1], str(sumrow[2]), str(sumrow[3])]], axis=0)
                else:
                    write_row = np.array([[sumrow[0], sumrow[1], str(sumrow[2]), str(sumrow[3])]])
            sumrow = [seq_i, tnt_i, 0, 0]
            numct = []
            proto_prev = proto_i
            # prim_prev, boo_prev, num_prev = None, None, None
            num_prms_prim3, num_prms_final = [0]*max_num_tgts, [0]*max_num_tgts
        if ((prim_i == prim_prev) and (seq_i==seq_prev)):
            num_prms_prim3[tgt_i] += 1
            if ((int(boo_i)==int(boo_prev)==1) and (int(num_i)==int(num_prev)==1)):
                num_prms_final[tgt_i] += 1
        numct.append(num_i)
        boo_prev, prim_prev, num_prev, seq_prev = boo_i, prim_i, num_i, seq_i
    if len(numct) > 0:
        count_prim3 = len([i for i in num_prms_prim3 if i >= min_prmrs])
        count_final = len([i for i in num_prms_final if i >= min_prmrs])
        sumrow[2] = count_prim3
        sumrow[3] = count_final
        if write_row is not None:
            write_row = np.append(write_row, [[sumrow[0], sumrow[1], str(sumrow[2]), str(sumrow[3])]], axis=0)
        else:
            write_row = np.array([[sumrow[0], sumrow[1], str(sumrow[2]), str(sumrow[3])]])
        write_row = write_row[np.lexsort((write_row[:,2].astype(int), write_row[:,3].astype(int)))]
    for i in write_row:
        csv_out.write(','.join(i) + '\n')
    csv_out.close()
######################################################################################################################
######################################################################################################################
####################### DON'T TOUCH FUNCTIONS ########################################################################
######################################################################################################################
######################################################################################################################
def genome_initialized(savepath, genome_str):
    """ Initialize downloading of genome GI indices, full genome, and sizes as global variables
    :param savepath: path to save the genome
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    """
    global GENOME, GENOME_LIST, genome_seq, genome_id, genome_size
    if genome_str in GENOME_LIST and genome_str != GENOME:
        GENOME = genome_str
        genome_seq = msa.get_genome_seq(genome_str, savepath)
        genome_id = msa.get_genome_ids(genome_str)
        genome_size = c.get_genome_dict(genome_str)

def get_target_region(chrom, coord, sen_i, win):
    """ Helper function for get_target_sequences(). Returns the sequences given the chromosome,
        coordinate, and desired window width. Also orients sequences along the same sense (+).
    :param chrom: chromosome stirng, e.g., 'chr1'
    :param coord: coordinate integer, e.g., 1050223
    :param sen_i: the orientation/sense of the sequence, either '+' or '-'
    :param win: window width. The sequence region width will be win * 2 + 1.
    :return: a Seq object containing the genome sequence for region of interest.
    """
    global genome_seq
    chr_i = genome_seq[chrom]
    chr_i_len = len(chr_i)
    if sen_i == '+':
        pe_lt = coord - win + 15
        pe_rt = coord + win + 15 + 1
    else:
        pe_lt = coord - win + 5
        pe_rt = coord + win + 5 + 1
    if pe_lt >= 0 and pe_rt < chr_i_len:
        read = chr_i[pe_lt:pe_rt]
        s = read.seq if sen_i == '+' else read.reverse_complement().seq
        return s
    return ""

def rev_complement(seq):
    bases = list(seq)
    bases = [complement[base] for base in bases]
    new_seq = ''.join(bases)
    return new_seq[::-1]