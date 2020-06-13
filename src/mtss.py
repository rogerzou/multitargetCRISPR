# -*- coding: utf-8 -*-
""" # ChIP-seq analysis after Cas9 cleavage by multi-targeting guides
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from collections import defaultdict
import pickle
import pysam
import os
import subprocess as sp
import statistics
import re
import numpy as np
import csv
from scipy import stats
from . import chipseq as c
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.inspection import permutation_importance
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, roc_auc_score
from lib.pyliftover import liftover

from matplotlib import pyplot as pp

OLD_CHAR = "ACGT"
NEW_CHAR = "TGCA"
WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3]


def find_msa(generator, bamfile, outfile, hg38):
    """ For each target site, find alternative alignments for all paired-end reads in window.

    :param generator: generator that outputs target sites in the following tuple format:
                    ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                      cut_i     =   cut site                 (int)
                      sen_i     =   sense/antisense          (+/- str)
                      pam_i     =   PAM                      (str)
                      gui_i     =   genomic target sequence  (str)
                      mis_i     =   # mismatches             (int)
                      guide     =   intended target sequence (str)
    :param bamfile: BAM file with paired-end ChIP-seq reads to analyze using generator target sites
    :param outfile: output file (extension omitted)
    :param hg38: path to hg38 genome (with .fa extension)

    Default -fr alignments
    {83, 163} = first/second in pair + first is reverse + primary alignment
    {99, 147} = first/second in pair + second is reverse + primary alignment
    {339, 419} = first/second in pair + first is reverse + NOT primary alignment
    {355, 403} = first/second in pair + second is reverse + NOT primary alignment

    """
    fullflags = ['83', '163', '99', '147', '339', '419', '355', '403']
    initflags = ['83', '99', '339', '355']
    primflags = ['83', '99']
    bamin = pysam.AlignmentFile(bamfile, 'rb')
    gen_counter = -1
    csv_out = open(outfile + '.csv', 'w')
    for rs, cut, sen, pam, gui, mis, tar in generator:    # iterate over each target site
        """ Genome-wide alignment of all paired-end reads around each target site """
        # get window centered at Cas9 cut site
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        sta_i, end_i = int(sta_i), int(end_i)
        # counter to track progress through target site generator
        gen_counter += 1
        if gen_counter % 10 == 0:
            print(gen_counter)
        # write all reads in target site window to paired fastq files in original '-fr' format
        with open(outfile + '1.fq', 'w') as fq1, open(outfile + '2.fq', 'w') as fq2:
            fq_counter = -1
            for read1, read2 in c.read_pair_generator(bamin, rs):
                fq_counter += 1
                s1, q1, s2, q2 = _find_msa_helper(read1, read2)
                fq1.write("@%s_%05i\n%s\n+\n%s\n" % (rs, gen_counter, s1, q1))
                fq2.write("@%s_%05i\n%s\n+\n%s\n" % (rs, gen_counter, s2, q2))
        # align all paired-end reads to genome (hg38) using bowtie2
        sp.run(['bowtie2', '-p', '4', '--local', '-k', '300', '-X', '1000', '-x', hg38[:-3],
                '-1', outfile + '1.fq', '-2', outfile + '2.fq', '-S', outfile + '.sam'])
        """ Parse bowtie2 samfile output, determine all primary/secondary alignments """
        currow, scores = None, []
        sam_it = open(outfile + '.sam', 'r')
        sam_counter = -1
        for read in sam_it:                             # read every alignment of SAM file
            if read[0] == '@':                          # skip SAM header lines
                continue
            stp = read.strip()
            row = stp.split('\t')                       # read each bowtie2 alignment in SAM format
            if row[1] in initflags:                     # first in pair, properly mapped
                # get quality scores and location of alignment
                AS = int(re.findall(r'AS:i:[-0-9]+', stp)[0].split(':')[2])
                YS = int(re.findall(r'YS:i:[-0-9]+', stp)[0].split(':')[2])
                sco_c = (AS + YS) / 2.0
                chr_c, aln_c = row[2], int(row[3])
                str_c = "%s_%i_%0.3f" % (chr_c, aln_c, sco_c)
                if row[1] in primflags:                 # first in pair, primary alignment
                    # save previous alignment set
                    if currow is not None:
                        num_c = len(scores)
                        max_c = '1' if currow[4] and (num_c == 1 or scores[0] > max(scores[1:])) else '0'
                        currow[4] = str(max_c)
                        currow.insert(5, str(num_c))
                        csv_out.write(','.join(currow) + '\n')
                    # determine if primary alignment matches previously assigned location
                    sam_counter += 1
                    boo = True if chr_c == chr_i and sta_i - 1e3 <= aln_c <= end_i + 1e3 else False
                    currow = [rs, str(cut), str(gen_counter), str(sam_counter), boo]
                    scores = []
                currow.append(str_c)
                scores.append(sco_c)
            elif row[1] not in fullflags:
                raise ValueError("find_msa(): unexpected SAM flag: %s" % row[1])
        # save final alignment set
        if currow is not None:
            num_c = len(scores)
            max_c = '1' if currow[4] and (num_c == 1 or scores[0] > max(scores[1:])) else '0'
            currow[4] = str(max_c)
            currow.insert(5, str(num_c))
            csv_out.write(','.join(currow) + '\n')
        sam_it.close()
    csv_out.close()
    bamin.close()
    os.remove(outfile + '1.fq')
    os.remove(outfile + '2.fq')
    os.remove(outfile + '.sam')


def _find_msa_helper(read1, read2):
    """ Extract read pair sequence and quality from original paired-end ChIP-seq data

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return read1 sequence, read1 quality, read2 sequence, read2 quality from original ChIP-seq data

    """
    if read1.is_reverse and not read2.is_reverse:
        return get_reverse_complement(read1.seq), read1.qual[::-1], read2.seq, read2.qual
    elif read2.is_reverse and not read1.is_reverse:
        return read1.seq, read1.qual, get_reverse_complement(read2.seq), read2.qual[::-1]
    else:
        raise ValueError("_find_msa_helper(): Unexpected paired-end read conditions.")


def get_targets_fasta(outfile, seqstr, numbases):
    """ Local exhaustive search for all possible protospacer sequences, up to 3 mismatches
    at the PAM-proximal side, from an input sequence.

    :param outfile: name of output fasta file (with .fa extension)
    :param seqstr: input sequence as a string
    :param numbases: number of bases counting from PAM-proximal side available to mutate

    """
    seqstr = seqstr.upper()
    unique_ele = set()
    with open(outfile, 'w') as f:       # only write unique sequences to file
        # check all 20bp substrings from input sequence
        for i in range(len(seqstr)-19):
            seq_i = seqstr[i:i+20]
            print("get_targets_fasta() - Sequence #%i: %s" % (i, seq_i))
            for seqmod in _get_targets_fasta_helper(seq_i, numbases):
                if seqmod not in unique_ele:
                    f.write(">%s\n%s\n" % (seqmod, seqmod))
                    unique_ele.add(seqmod)
        # check all 20bp substrings from reverse complement of input sequence
        seqstr = get_reverse_complement(seqstr)
        for i in range(len(seqstr)-19):
            seq_i = seqstr[i:i+20]
            print("get_targets_fasta() - Sequence #%i: %s" % (i, seq_i))
            for seqmod in _get_targets_fasta_helper(seq_i, numbases):
                if seqmod not in unique_ele:
                    f.write(">%s\n%s\n" % (seqmod, seqmod))
                    unique_ele.add(seqmod)


def _get_targets_fasta_helper(init_str, numbases):
    """ Helper function to yield all sequences with up to 3 mismatches from template.
    Uniqueness is not checked, but the GC content is restricted to between 0.4 and 0.7

    :param init_str: template sequence string
    :param numbases: number of bases counting from PAM-proximal side available to mutate

    """
    init_list = list(init_str)
    for i in range(1, numbases - 1):
        for j in range(i + 1, numbases):
            for k in range(j + 1, numbases + 1):
                for nb_i in {'A', 'C', 'G', 'T'}:
                    for nb_j in {'A', 'C', 'G', 'T'}:
                        for nb_k in {'A', 'C', 'G', 'T'}:
                            mod_ijk = list(init_list)
                            mod_ijk[-i] = nb_i
                            mod_ijk[-j] = nb_j
                            mod_ijk[-k] = nb_k
                            seq_str = "".join(mod_ijk)
                            if 0.4 < get_gc(seq_str) < 0.7:      # restrict to 40-70% GC content
                                yield seq_str + "NGG"


def get_gc(seq_str):
    """ Return GC content from a sequence string. """
    return float(seq_str.count('G') + seq_str.count('C')) / len(seq_str)


def get_targets_bowtie2(infile, outfile, hg38):
    """ Run bowtie2 to align 'infile' fasta file to hg38, result in 'outfile' sam file """
    sp.run(['bowtie2', '-k', '1000', '-f', '-x', hg38, '-U', infile, '-S', outfile])


def get_targets_stats(samfile, fileout):
    """ Given a bowtie2 samfile output, outputs:
    - number of alignments and whether it's in a gene for each sequence (*_counts.csv)
    - number of alignments in a REFSEQ gene for each sequence (*_align.csv)

    :param samfile: samfile output from bowtie2
    :param fileout: output file name/path (no extension)

    """
    sam_it = open(samfile, 'r')
    cnt_list, cnt_set, sam_list, p_coun, p_gene, p_read, cter = [], set(), [], 0, 0, "", 0
    for read in sam_it:     # read every line of SAM file
        # counter to track progress
        if cter % 10000 == 0:
            print("get_targets_stats(): Processed %i alignments." % cter)
        cter += 1
        # skip SAM header lines
        if read[0] == '@':
            continue
        # read each bowtie2 alignment in SAM format
        row = read.strip().split('\t')
        if row[1] == '0' or row[1] == '16':         # first alignment for new sequence
            # determine new alignment for new sequence
            chr_i, cor_i = row[2], int(row[3])
            ig = c.is_gene_refseq(chr_i, cor_i)
            # record previous sequence alignment counts, reset to new sequence
            if p_coun != 0:
                if p_read in cnt_set:
                    raise ValueError("get_targets_stats(): duplicate and separate alignments")
                cnt_list.append([p_read, p_coun, p_gene])
                cnt_set.add(p_read)
            p_coun, p_read = 1, row[0]
            p_gene = 1 if ig else 0
            # record new alignment for new sequence
            nam_i, sen_i = (ig[0], ig[1]) if ig else ("", "")
            sense = '+' if row[1] == '0' else '-'
            sam_list.append([row[0], chr_i, cor_i, nam_i, sen_i, sense])
        elif row[1] == '4':                         # no alignments for new sequence
            # record previous sequence alignment counts, reset to no sequence
            if p_coun != 0:
                if p_read in cnt_set:
                    raise ValueError("get_targets_stats(): duplicate and separate alignments")
                cnt_list.append([p_read, p_coun, p_gene])
                cnt_set.add(p_read)
            p_coun, p_gene, p_read = 0, 0, ""
        elif row[1] == '256' or row[1] == '272':    # more alignments for current sequence
            # determine new alignment for current sequence
            chr_i, cor_i = row[2], int(row[3])
            ig = c.is_gene_refseq(chr_i, cor_i)
            # update alignment counts
            p_coun += 1
            p_gene = p_gene + 1 if ig else p_gene
            # record new alignment for current sequence
            nam_i, sen_i = (ig[0], ig[1]) if ig else ("", "")
            sense = '+' if row[1] == '256' else '-'
            sam_list.append([row[0], chr_i, cor_i, nam_i, sen_i, sense])
        else:                                       # unexpected alignment flag
            raise ValueError("get_targets_stats(): unexpected flag value of %s." % row[1])
    # record last sequence alignment counts
    if p_coun != 0:
        if p_read in cnt_set:
            raise ValueError("get_targets_stats(): duplicate and separate alignments")
        cnt_list.append([p_read, p_coun, p_gene])
        cnt_set.add(p_read)
    # finalize and save to file
    sam_it.close()
    cntnp = np.asarray(cnt_list)
    cntnp = cntnp[cntnp[:, 1].astype(int).argsort()[::-1], :]   # sort descending by counts
    np.savetxt(fileout + '_count.csv', cntnp, fmt='%s', delimiter=',')
    np.savetxt(fileout + '_align.csv', np.asarray(sam_list), fmt='%s', delimiter=',')


def get_targets_dist(alignfile, fileout):
    """ Given get_targets_stats() *_align.csv output, outputs:
    - average distance between adajacent putative on-target sites (*_dist.csv)

    :param alignfile: *_align.csv output from get_targets_stats()
    :param fileout: output file name/path (no extension)

    """
    aln_np = load_nparray(alignfile)
    curseq, curind, dist_list = None, [], []
    for i in range(aln_np.shape[0]):       # iterate over alignlist
        row = aln_np[i, :]
        if curseq != row[0]:               # a new sequence
            if curseq is not None:         # and not the very first one in list
                avdist, numalign = _get_targets_dist_helper(aln_np[curind, :])
                dist_list.append([curseq, numalign, avdist])   # calculate avg dist for previous seq
            curind, curseq = [i], row[0]   # for the new sequence, record its indices from alignlist
        else:
            curind.append(i)               # if same sequence as previous row, add its index
    avdist, numalign = _get_targets_dist_helper(aln_np[curind, :])
    dist_list.append([curseq, numalign, avdist])               # calculate avg dist for the last seq
    distnp = np.asarray(dist_list)
    distnp = distnp[distnp[:, 1].astype(int).argsort()[::-1], :]   # sort descending by counts
    np.savetxt(fileout + '_dist.csv', distnp, fmt='%s', delimiter=',')


def _get_targets_dist_helper(aln):
    """ Helper function to calculate average distance between adjacent putative on-target sites

    :param aln: a subset of alignfile (*_align.csv) rows that correspond to all putative on-target
    sites for a particular target sequence

    """
    aln_sorted = aln[np.lexsort((aln[:, 2].astype(int), aln[:, 1]))]
    numalign, avdist = aln_sorted.shape[0], []
    if numalign > 1:
        for i in range(aln_sorted.shape[0]-1):
            aln1, aln2 = aln_sorted[i, :], aln_sorted[i+1, :]
            if aln1[1] == aln2[1]:
                avdist.append(int(aln2[2]) - int(aln1[2]))
    if len(avdist) > 0:
        return statistics.mean(avdist), numalign
    else:
        return None, numalign


def get_reverse_complement(seq):
    """ Return reverse complement of sequence string. """
    return seq.translate(str.maketrans(OLD_CHAR, NEW_CHAR))[::-1]


def load_nparray(array):
    """ Load dataset as numpy array. """
    return np.loadtxt(array, dtype=object, delimiter=',')


def load_chakrabarti():
    """ Load dataset from Chakrabarti et al., Mol Cell, 2019 """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return np.loadtxt(os.path.dirname(dirname) + "/lib/chakrabarti/chakrabarti_mm.csv",
                      delimiter=',', dtype=object, skiprows=1)


def load_liftover():
    """ Load class instance for liftover from hg19 to hg38
    For use with Chakrabarti et al., which aligned to hg19, but hg38 is used for this analysis """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return liftover.LiftOver(open(os.path.dirname(dirname) + "/lib/hg19ToHg38.over.chain"))


def target_gen(alignfile, span_r, guide):
    """ Generator to yield all putative on-target sites for a given protospacer

    :param alignfile: CSV file generated from get_targets_stats() output
    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param guide: on-target protospacer sequence (no PAM)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM, discovered protospacer,
              # mismatches, non-mismatched protospacer )
            - The format is designed to be consistent with other generators that may yield
              mismatched sequences. Note that no mismatched sequences will be outputted by design.

    """
    aln = load_nparray(alignfile)
    hg38size = c.hg38_dict()
    pam_i = 'NGG'
    for i in range(aln.shape[0]):
        row = aln[i, :]
        if row[0] == guide + pam_i:
            chr_i = row[1]
            sen_i = row[5]
            if sen_i == '+':
                cut_i = int(row[2]) + 16
            else:
                cut_i = int(row[2]) + 6
            span_sta = max(1, cut_i - span_r)
            span_end = min(hg38size[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            mis_i = 0
            yield span_rs, cut_i, sen_i, pam_i, guide, mis_i, guide


def macs_gen(peak, span_r, genome, guide, mismatch=20, cent_r=200):
    """ Generator to yield all peaks from macs2 output.

    :param peak: macs2 output file
    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param genome: path to genome (hg38)
    :param guide: on-target protospacer sequence (no PAM)
    :param mismatch: number of mismatches from protospacer to accept
    :param cent_r: radius of window from peak center to search for cut site
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM, discovered protospacer,
            # mismatches, non-mismatched protospacer )

    """
    m_out = np.loadtxt(peak, dtype=object)    # load macs2 narrowPeak output file
    numpeaks = m_out.shape[0]                 # number of peaks from blender
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        chr_i = m_out[i, 0]
        if chr_i in hg38size:
            center = int(m_out[i, 9]) + int(m_out[i, 1])
            cent_sta = max(1, center - cent_r)
            cent_end = min(hg38size[chr_i], center + cent_r)
            cent_rs = "%s:%i-%i" % (chr_i, cent_sta, cent_end)
            cent_faidx = sp.check_output(['samtools', 'faidx', genome, cent_rs]).split()
            seq = (b"".join(cent_faidx[1:]).upper()).decode("utf-8")
            cand = sub_findmis(seq, guide, mismatch)
            if cand is not None and len(cand) > 0:
                cut_i = cent_sta + cand[0][2]
                sen_i = '+' if cand[0][5] == 1 else '-'
                pam_i = cand[0][4]
                gui_i = cand[0][3]
                mis_i = cand[0][1]
                span_sta = max(1, cut_i - span_r)
                span_end = min(hg38size[chr_i], cut_i + span_r)
                span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide


def chakrabarti_generator(span_r, genome):
    """ Generator to yield all target sites from Chakrabarti et al., Mol Cell, 2019.

    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param genome: path to genome (hg38)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM , discovered protospacer,
            # mismatches, non-mismatched protospacer )

    """
    d = load_chakrabarti()
    hg38size = c.hg38_dict()
    lo = load_liftover()
    numtargets = d.shape[0]
    for i in range(numtargets):
        guide = d[i, 1][:-3]
        pam_init = d[i, 1][-3:]
        chr_i = d[i, 3]
        sta_tmp = int(d[i, 4])
        end_tmp = int(d[i, 5])
        sta_tmp = lo.convert_coordinate(chr_i, sta_tmp)[0][1]
        end_tmp = lo.convert_coordinate(chr_i, end_tmp)[0][1]
        cent_sta = max(1, sta_tmp - 250)
        cent_end = min(hg38size[chr_i], end_tmp + 250)
        cent_rs = "%s:%i-%i" % (chr_i, cent_sta, cent_end)
        cent_faidx = sp.check_output(['samtools', 'faidx', genome, cent_rs]).split()
        seq = (b"".join(cent_faidx[1:]).upper()).decode("utf-8")
        cand = sub_findmis(seq, guide, 0)
        if cand is not None and len(cand) > 0:
            cut_i = cent_sta + cand[0][2]
            sen_i = '+' if cand[0][5] == 1 else '-'
            gui_i = cand[0][3]
            mis_i = cand[0][1]
            pam_i = cand[0][4]
            if pam_i != pam_init:
                print("chakrabarti_generator(): PAM DIFFERENT!!!!!")
            span_sta = max(1, cut_i - span_r)
            span_end = min(hg38size[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide
        else:
            print("chakrabarti_generator(): Did not find a particular sequence.")


def sub_findmis(s, matchstr, maxmismatch):
    """ Find the best protospacer (with mismatch tolerance) in a sequence range

    :param s: sequence string in which to find match
    :param matchstr: match sequence string
    :param maxmismatch: mismatch tolerance (int)
    """
    sublen = len(matchstr)
    len_full = sublen + 3
    candidates = []
    # check mismatches for original string
    for i in range(len(s)-sublen-6):
        i += 3
        # sense
        substr_1 = s[i:i + sublen]
        pam_1 = s[i + sublen:i + sublen + 3]
        cutsite_1 = i + sublen - 4
        str_1 = substr_1 + pam_1
        matchstr_1 = matchstr + 'NGG'
        wlist = [x + 1 if str_1[j] != matchstr_1[j] else x for j, x in enumerate([0] * len_full)]
        mismatches = sum(wlist) - 1
        if mismatches <= maxmismatch:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_1, substr_1, pam_1, 1))
        # anti-sense
        substr_0 = get_reverse_complement(s[i:i + sublen])
        pam_0 = get_reverse_complement(s[i - 3:i])
        cutsite_0 = i + 3
        str_0 = substr_0 + pam_0
        matchstr_0 = matchstr + 'NGG'
        wlist = [x + 1 if str_0[j] != matchstr_0[j] else x for j, x in enumerate([0] * len_full)]
        mismatches = sum(wlist) - 1
        if mismatches <= maxmismatch:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_0, substr_0, pam_0, 0))

    candidates.sort(key=lambda y: y[0])
    return candidates


def peak_profile(generator, bamfilein, fileout):
    """ For each target location outputted from generator, output the flanking region in
    WIG and CSV format

    :param bamfilein:
    :param generator:
    :param fileout:

    """
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    chr_old, csv_peaks = None, []
    wlist_all = []
    for rs, cut, sen, pam, gui, mis, guide in generator:
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        sta_i = int(sta_i)
        end_i = int(end_i)
        wlist = [0] * (end_i - sta_i + 1)
        for read1, read2 in c.read_pair_generator(bamin, rs):
            read = c.read_pair_align(read1, read2)
            wlist = [x + 1 if read[0] - sta_i <= j <= read[-1] - sta_i else x for j, x in
                     enumerate(wlist)]
        wlist = [x / bamin.mapped * 1E6 for x in wlist]
        wlist_all.append([chr_i, sta_i] + wlist)
        csv_peaks.append(wlist) if sen == '+' else csv_peaks.append(wlist[::-1])
    bamin.close()
    with open(fileout + "_bpeaks.wig", 'w') as wigout:
        wlist_all = np.asarray(wlist_all)
        wlist_all = wlist_all[np.lexsort((wlist_all[:, 1].astype(int), wlist_all[:, 0])), :]
        chr_prev = None
        for i in range(wlist_all.shape[0]):
            row = wlist_all[i, :]
            chr_i = row[0]
            if chr_prev != chr_i:
                chr_prev = chr_i
                wigout.write("variableStep\tchrom=%s\n" % chr_i)
            sta_i = int(row[1])
            wlist = row[2:].astype(float)
            for j, x in enumerate(wlist):
                wigout.write("%i\t%0.5f\n" % (sta_i + j, x))
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',')


def read_subsets(generator, filein, fileout):
    bamin = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    bamspan = pysam.AlignmentFile(fileout + "_span.bam", 'wb', template=bamin)
    bamabut = pysam.AlignmentFile(fileout + "_abut.bam", 'wb', template=bamin)
    bamelse = pysam.AlignmentFile(fileout + "_else.bam", 'wb', template=bamin)
    LS = []
    for rs, cut, sen, pam, gui, mis, tar in generator:
        ctM, ctN, ctL, ctR, ctT = 0, 0, 0, 0, 0
        for read1, read2 in c.read_pair_generator(bamin, rs):
            read = c.read_pair_align(read1, read2)
            if not read:
                continue
            ctT += 1
            if read[0] < cut < read[3]:         # fragments that span
                ctM += 1
                bamspan.write(read1)
                bamspan.write(read2)
            elif cut + 5 >= read[0] >= cut:     # fragments that begin 5bp of cleavage site
                ctR += 1
                ctN += 1
                bamabut.write(read1)
                bamabut.write(read2)
            elif cut - 5 <= read[-1] <= cut:    # fragments that end 5bp of cleavage site
                ctL += 1
                ctN += 1
                bamabut.write(read1)
                bamabut.write(read2)
            else:                               # other fragments
                bamelse.write(read1)
                bamelse.write(read2)
        fpm = bamin.mapped / 2E6      # fragments per millon
        ctM /= fpm                  # fragments that span
        ctN /= fpm                  # fragments that don't span
        ctL /= fpm                  # fragments that don't span, on left
        ctR /= fpm                  # fragments that don't span, on right
        ctT /= fpm                  # count total number of reads
        chr_i, sta_i, end_i = c.region_string_split(rs)
        ig = c.is_gene_refseq(chr_i, cut)
        cUp, cDo, cPr, cDi = None, None, None, None
        if sen == '+':
            cDi, cPr = ctL, ctR
        if sen == '-':
            cDi, cPr = ctR, ctL
        if ig and ig[1] == '+':
            cUp, cDo = ctL, ctR
        if ig and ig[1] == '-':
            cUp, cDo = ctR, ctL
        if ig:
            LS.append((rs, cut, sen, gui, pam, tar, mis, ig[0], ig[1], ctT, ctM, ctN, cUp, cDo, cPr, cDi))
        else:
            LS.append((rs, cut, sen, gui, pam, tar, mis, "", "", ctT, ctM, ctN, "", "", cPr, cDi))
    bamin.close()
    bamspan.close()
    bamabut.close()
    bamelse.close()
    file_array = [fileout + "_span.bam", fileout + "_abut.bam", fileout + "_else.bam"]
    for file_i in file_array:
        pysam.sort("-o", file_i, file_i)
        os.system("samtools index " + file_i)
    header = "region string, cut site, Cas9 sense, observed target sequence, PAM, " \
             "expected target sequence, mismatches, RefSeq, RefSeq sense, " \
             "Ctotal, Cspan, Cend, Cupstream, Cdownstream, Cproximal, Cdistal"
    LS = np.asarray(LS)
    np.savetxt(fileout + ".csv", LS, fmt='%s', delimiter=',', header=header)


def read_counts(generator, filein, fileout=None):
    bam = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    list_stat = []
    for rs, cut, sen, pam, gui, mis, tar in generator:
        ct_rpm = bam.count(region=rs) / bam.mapped * 1E6
        list_stat.append((rs, cut, sen, gui, mis, ct_rpm))
    bam.close()
    list_stat = np.asarray(list_stat)
    if fileout:
        np.savetxt(fileout, list_stat, fmt='%s', delimiter=',')
    return list_stat


def mergesubsetcounts(subset, countlists, fileout):
    for x in countlists:
        subset = np.column_stack((subset, x[:, -1]))
    np.savetxt(fileout, subset, fmt='%s', delimiter=',')
    return subset


def mergerows(files, fileout):
    merged = files[0]
    for i in range(len(files)-1):
        merged = np.row_stack((merged, files[i+1]))
    np.savetxt(fileout, merged, fmt='%s', delimiter=',')
    return merged


def getXy_nomismatch(data):
    fl = np.loadtxt(data, dtype=object, delimiter=',')
    fl = fl[fl[:, 5] == '0', :]   # get non-mismatched columns only
    np.savetxt(data[0:-4] + "_0.csv", fl, fmt='%s', delimiter=',')
    y = fl[:, 8].astype(float)
    integer_encoded = LabelEncoder().fit_transform(fl[:, 4])
    onehot_encoded = OneHotEncoder(sparse=False).fit_transform(integer_encoded.reshape(-1, 1))
    X = np.column_stack((onehot_encoded, fl[:, 15:].astype(float)))
    labels = ["onehot1", "onehot2", "onehot3", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac",
              "H3K36me3", "DNase I", "MNase-seq", "ATAC-seq", "RNA-seq"]
    return X, y, labels


def getXy_all(data):
    fl = np.loadtxt(data, dtype=object, delimiter=',')
    y = fl[:, 8].astype(float)
    integer_encoded = LabelEncoder().fit_transform(fl[:, 4])
    onehot_encoded = OneHotEncoder(sparse=False).fit_transform(integer_encoded.reshape(-1, 1))
    X = np.column_stack((onehot_encoded, fl[:, 5].astype(float), fl[:, 15:]))
    labels = ["onehot1", "onehot2", "onehot3", "# mismatches", "H3K4me1", "H3K4me3", "H3K9me3",
              "H3K27ac", "H3K36me3", "DNase I", "MNase-seq", "ATAC-seq", "RNA-seq"]
    return X, y, labels


def getXy_chak(data, epi=2, index=0):
    fl = np.loadtxt(data, dtype=object, delimiter=',')
    if index == 0:
        y = fl[:, 9].astype(float)
    elif index == 1:
        y = fl[:, 10].astype(int)
    elif index == 2:
        y = fl[:, 11]
    else:
        raise ValueError("getXy_chak: index value should be {0, 1, 2}")
    seqs = fl[:, 1]
    nfeat = 21
    one_hot_encoded = np.zeros((fl.shape[0], nfeat*4))
    mapping = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
    for i in range(fl.shape[0]):
        seq = seqs[i]
        for j in range(nfeat):
            one_hot_encoded[i, 4*j:4*(j+1)] = mapping[seq[j]]
    epigenetic_encoded = fl[:, 12:]
    onehotlabels = []
    for j in range(nfeat):
        onehotlabels[4*j:4*(j+1)] = [str(20-j)+'_A', str(20-j)+'_C', str(20-j)+'_G', str(20-j)+'_T']

    if epi == 0:
        X = one_hot_encoded
        labels = onehotlabels
    elif epi == 1:
        X = epigenetic_encoded
        labels = ["H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K36me3", "DNase I", "MNase-seq",
                  "ATAC-seq", "RNA-seq"]
    elif epi == 2:
        X = np.column_stack((one_hot_encoded, epigenetic_encoded))
        labels = onehotlabels + \
                 ["H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K36me3", "DNase I", "MNase-seq",
                  "ATAC-seq", "RNA-seq"]
    else:
        raise ValueError("getXy_chak: epi value should be {0, 1, 2}")
    return X, y, labels


def pca(X, num_components=6):
    pca = PCA(n_components=num_components)
    X = pca.fit_transform(X)
    return X


def NeuralNetworkTrainDefault(X, y, modelfile, solver='lbfgs', classifier=False,
                              alpha=0.1, hidden_layer_sizes=(8,)):
    # Define and fit base regressor
    params = {'solver': solver, 'max_iter': 5000, 'random_state': 42,
              'hidden_layer_sizes': hidden_layer_sizes, 'alpha': alpha}
    if classifier:
        mlp = MLPClassifier(**params)
    else:
        mlp = MLPRegressor(**params)
    mlp.fit(X, y)
    evaluate(X, y, mlp, classifier)
    # save and return base estimator
    pickle.dump(mlp, open(modelfile, 'wb'))
    return mlp


def NeuralNetworkTrainGridCV(X, y, modelfile, params, classifier=False):
    # Define classifier or regressor class
    if classifier:
        mlp = MLPClassifier(max_iter=5000, random_state=42)
    else:
        mlp = MLPRegressor(max_iter=5000, random_state=42)
    # Find best hyperparameters and then refit on all training data
    mlp_grid = GridSearchCV(estimator=mlp, param_grid=params, cv=10, verbose=2, n_jobs=8)
    mlp_grid.fit(X, y)
    # get the best estimator and calculate results (correlation coefficient)
    best_grid = mlp_grid.best_estimator_
    evaluate(X, y, best_grid, classifier)
    # save and return best estimator
    pickle.dump(best_grid, open(modelfile, 'wb'))
    return best_grid


def RandomForestTrainDefault(X, y, modelfile):
    # Define and fit base regressor:
    rf = RandomForestRegressor(n_estimators=500, random_state=42)
    rf.fit(X, y)
    # calculate results (correlation coefficient)
    CORREL(rf, X, y)
    # save and return base estimator
    pickle.dump(rf, open(modelfile, 'wb'))
    return rf


def RandomForestTrainGridCV(X, y, modelfile):
    # Define base regressor:
    rf = RandomForestRegressor()
    # Define search space:
    params = {
        'n_estimators': [1400],
        'criterion': ['mse', 'mae'],
        'max_depth': [int(x) for x in np.linspace(10, 110, num=11)],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
        'max_features': ['auto'],
        'bootstrap': [True, False]
    }
    # Random search of parameters, using 5 fold cross validation,
    # search across 100 different combinations, and use all available cores
    rf_random = GridSearchCV(estimator=rf, param_grid=params, cv=10, verbose=2, n_jobs=8)
    rf_random.fit(X, y)
    # get the best estimator and calculate results (correlation coefficient)
    best_random = rf_random.best_estimator_
    CORREL(best_random, X, y)
    # save and return best estimator
    pickle.dump(best_random, open(modelfile, 'wb'))
    return best_random


def ModelTest(X, y, modelfile, classifier=False):
    model = pickle.load(open(modelfile, 'rb'))
    print(modelfile)
    print(model)
    evaluate(X, y, model, classifier)
    np.savetxt(modelfile[0:-4] + "_out.csv", np.column_stack((y, model.predict(X))), fmt='%s', delimiter=',')


def evaluate(X_test, y_test, model, classifier):
    y_pred = model.predict(X_test)
    if classifier:
        y_prob = model.predict_proba(X_test)
        print(confusion_matrix(y_test, y_pred))
        AUC_ROC(y_test, y_prob)
    else:
        CORREL(y_test, y_pred)


def FeatureImportance(X, y, modelfile, labels, count=10):
    model = pickle.load(open(modelfile, 'rb'))
    labels = np.asarray(labels, dtype=object)

    result = permutation_importance(model, X, y, n_repeats=10, random_state=42, n_jobs=8)
    sorted_idx = result.importances_mean.argsort()

    fig, ax = pp.subplots()
    importances_i = result.importances[sorted_idx].T[:, -count:]
    labels_i = labels[sorted_idx][-count:]
    ax.boxplot(importances_i, vert=False, labels=labels_i)
    ax.set_title("Permutation Importances (test set)")
    fig.tight_layout()
    pp.show()


def data_split(X, y, test_size=0.3):
    return train_test_split(X, y, test_size=test_size, shuffle=True, random_state=42)


def MAPE(y_test, y_pred):
    errors = abs(y_pred - y_test)
    mape = 100 * np.mean(errors / y_test)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy


def CORREL(y_test, y_pred):
    corr = np.corrcoef(y_test, y_pred)
    print('Model Performance')
    print('Correlation coefficient = {:0.2f}%.'.format(corr[1, 0]))
    return


def AUC_ROC(y_test, y_prob):
    macro_roc_auc_ovo = roc_auc_score(y_test, y_prob, multi_class="ovo", average="macro")
    weighted_roc_auc_ovo = roc_auc_score(y_test, y_prob, multi_class="ovo", average="weighted")
    macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
    weighted_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="weighted")
    print("One-vs-One ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
          "(weighted by prevalence)"
          .format(macro_roc_auc_ovo, weighted_roc_auc_ovo))
    print("One-vs-Rest ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
          "(weighted by prevalence)"
          .format(macro_roc_auc_ovr, weighted_roc_auc_ovr))


""" deprecated functions """

def blender_generator(mre11, radius):
    """ Generator

    :param mre11: bender output file (DISCOVER-seq)
    :param radius: radius of analysis window surrounding cut site

    """
    b_out = np.loadtxt(mre11, dtype=object, skiprows=1)    # load blender output file
    numpeaks = b_out.shape[0]                              # number of peaks from blender
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        chr_i = re.split('[:-]', b_out[i, 0])[0]
        cut_i = int(b_out[i, 1])
        sta = max(1, cut_i - radius)
        end = min(hg38size[chr_i], cut_i + radius)
        rs_i = "%s:%i-%i" % (chr_i, sta, end)
        sen_i = 1 if b_out[i, 4] == 'sense' else 0
        pam_i = b_out[i, 5]
        gui_i = b_out[i, 6]
        mis_i = b_out[i, 7]
        yield rs_i, cut_i, sen_i, pam_i, gui_i, mis_i


def blender2csv(mre11, bamfilein, fileout, radius):
    """ For each peak obtained with blender (from DISCOVER-seq), obtain values and statistics
    for a specified window, saved to CSV files

    :param mre11: bender output file (DISCOVER-seq)
    :param bamfilein: input BAM file. All possible output files will be taken from this file
    :param fileout:
    :param radius: radius of analysis window surrounding cut site

    """
    blender_output = np.loadtxt(mre11, dtype=object, skiprows=1)    # load blender output file
    bamin = pysam.AlignmentFile(bamfilein, 'rb')             # BAM file for analysis
    numpeaks = blender_output.shape[0]                              # number of peaks from blender
    csv_peaks = np.zeros((numpeaks, radius * 2 + 1))
    csv_info = np.zeros((numpeaks, 5), dtype=object)
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        [chr_i, sta_i, end_i] = re.split('[:-]', blender_output[i, 0])
        cut_i = int(blender_output[i, 1])
        sta = max(1, cut_i - radius)
        end = min(hg38size[chr_i], cut_i + radius)
        region_string_i = "%s:%i-%i" % (chr_i, sta, end)
        wlist = [0] * (end - sta + 1)
        for readi in bamin.fetch(region=region_string_i):
            read = [readi.positions[0] + 1, readi.positions[-1] + 1]
            wlist = [x + 1 if read[0] - sta <= j <= read[-1] - sta else x for j, x in
                     enumerate(wlist)]
        csv_peaks[i, :] = wlist
        csv_info[i, :] = [chr_i, sta_i, end_i, max(wlist), sum(wlist)]
        c.status_statement(i, numpeaks, 50, chr_i)
    bamin.close()
    np.savetxt(fileout + "_bpeaks.csv", csv_peaks, fmt='%s', delimiter=',')
    np.savetxt(fileout + "_binfo.csv", csv_info, fmt='%s', delimiter=',')


def blender_peaks(mre11, bamfilein, fileout, radius, bamBool=False, wigBool=False, csvBool=True):
    """ For each peak obtained with blender (from DISCOVER-seq), display only the flanking region
    from a BAM file, outputting BAM, WIG, and CSV files

    :param mre11: bender output file (DISCOVER-seq)
    :param bamfilein: input BAM file (without extension). All possible output files will be taken
                    from this template
    :param fileout: TODO
    :param radius: radius of analysis window surrounding cut site
    :param bamBool: boolean to indicate BAM file results output
    :param wigBool: boolean to indicate WIG file results output
    :param csvBool: boolean to indicate CSV file results output

    """
    blender_output = np.loadtxt(mre11, dtype=object, skiprows=1)    # load blender output file
    bamin = pysam.AlignmentFile(bamfilein + ".bam", 'rb')             # BAM file for analysis
    numpeaks = blender_output.shape[0]                              # number of peaks from blender
    if bamBool:
        bamoutdir = fileout + "_bpeaks.bam"
        bamout = pysam.AlignmentFile(bamoutdir, 'wb', template=bamin)
    if wigBool:
        wigout = open(fileout + "_bpeaks.wig", 'w')
        chr_old = None
    if csvBool:
        csv_peaks = np.zeros((numpeaks, radius * 2 + 1))
        csv_info = np.zeros((numpeaks, 5), dtype=object)
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        [chr_i, sta_i, end_i] = re.split('[:-]', blender_output[i, 0])
        cut_i = int(blender_output[i, 1])
        sta = max(1, cut_i - radius)
        end = min(hg38size[chr_i], cut_i + radius)
        region_string_i = "%s:%i-%i" % (chr_i, sta, end)
        if bamBool:                         # write to BAM file
            for read in bamin.fetch(region=region_string_i):
                bamout.write(read)
        if wigBool or csvBool:              # write to WIG file
            wlist = [0] * (end - sta + 1)
            for read1, read2 in c.read_pair_generator(bamin, region_string_i):
                read = c.read_pair_align(read1, read2)
                wlist = [x + 1 if read[0] - sta <= j <= read[-1] - sta else x for j, x in
                         enumerate(wlist)]
            if wigBool:
                if chr_i != chr_old:
                    wigout.write("variableStep\tchrom=%s\n" % chr_i)
                    chr_old = chr_i
                for j, x in enumerate(wlist):
                    wigout.write("%i\t%i\n" % (sta + j, x))
            if csvBool:
                csv_peaks[i, :] = wlist
                csv_info[i, :] = [chr_i, sta_i, end_i, max(wlist), sum(wlist)]
        c.status_statement(i, numpeaks, 50, chr_i)
    bamin.close()
    if csvBool:
        csv_index = np.flipud(np.argsort(csv_info[:, 3]))
        csv_peaks = csv_peaks[csv_index, :]
        csv_info = csv_info[csv_index, :]
        np.savetxt(fileout + "_bpeaks.csv", csv_peaks, fmt='%s', delimiter=',')
        np.savetxt(fileout + "_binfo.csv", csv_info, fmt='%s', delimiter=',')
        np.savetxt(fileout + "_bindex.csv", csv_index, fmt='%s', delimiter=',')
    if bamBool:
        bamout.close()
        pysam.sort("-o", bamoutdir, bamoutdir)
        os.system("samtools index " + bamoutdir)
    if wigBool:
        wigout.close()
