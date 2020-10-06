# -*- coding: utf-8 -*-
""" Functions for ChIP-seq analysis after Cas9 cleavage by multi-targeting guides
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

import pysam
import os
import subprocess as sp
import statistics
import re
import numpy as np
from . import chipseq as c


OLD_CHAR = "ACGT"
NEW_CHAR = "TGCA"
WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3]
CHROMHMM = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts',
            '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']


def get_reverse_complement(seq):
    """ Return reverse complement of sequence string. """
    return seq.translate(str.maketrans(OLD_CHAR, NEW_CHAR))[::-1]


def load_nparray(array):
    """ Load dataset as numpy array. """
    return np.loadtxt(array, dtype=object, delimiter=',')


def load_npheader(array):
    """ Load the header (first line) of numpy array file. """
    with open(array) as f:
        head = f.readline()
    return head[2:-1]


def sort_mmtype_column(datafile, collist):
    """ Given processed data file, generate separate csv file for each column listed in collist.
        For each individual csv file, each column corresponds to a different value of 'mmtype'

        :param datafile: processed CSV data file, i.e. from mergesubsetcounts() or read_subsets()
        :param collist: list of column labels from the processed data file, such that each column
                        label is split to its own csv file; the data are drawn from this column

        Creates a folder named after the datafile, with each file in the folder corresponding to a
        different column from collist
    """
    header, data = load_npheader(datafile), load_nparray(datafile)
    head = header.split(', ')
    ind_mmtype = head.index('mm_type')
    uni_mmtype = np.unique(data[:, ind_mmtype])
    if not os.path.exists(datafile[:-4]):
        os.makedirs(datafile[:-4])
    for col in collist:
        outlist = []
        col_ind = head.index(col)
        for t_i in uni_mmtype:
            outlist.append(data[data[:, ind_mmtype] == t_i, col_ind])
        outnp = np.empty((max([len(x) for x in outlist]), len(outlist)), dtype=object)
        for i, list_i in enumerate(outlist):
            outnp[0:len(list_i), i] = list_i
        np.savetxt(os.path.join(datafile[:-4], "%s.csv" % col), outnp,
                   fmt='%s', delimiter=',', header=",".join(uni_mmtype))


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
    :param hg38: path to hg38 genome (with .fa extension) for use by bowtie2

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
    with open(outfile + ".fa", 'w') as f:       # only write unique sequences to file
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


def get_targets_bowtie2(inf, hg38):
    """ Run bowtie2 to align 'inf' fasta file to hg38, result in 'outfile' sam file """
    sp.run(['bowtie2', '-k', '1000', '-f', '-x', hg38[:-3], '-U', inf + ".fa", '-S', inf + ".sam"])


def get_targets_stats(samfile):
    """ Given a bowtie2 samfile output, outputs:
    - number of alignments and whether it's in a gene for each sequence (*_counts.csv)
    - list of each alignment with the following columns (*_align.csv):
        protospacer + PAM | chr | coord | gene name | gene orientation | protospacer orientation

    :param samfile: samfile output from bowtie2 (no extension)

    """
    sam_it = open(samfile + ".sam", 'r')
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
    np.savetxt(samfile + '_count.csv', cntnp, fmt='%s', delimiter=',')
    np.savetxt(samfile + '_align.csv', np.asarray(sam_list), fmt='%s', delimiter=',')


def get_targets_dist(alignfile, fileout):
    """ Given get_targets_stats() *_align.csv output, outputs:
    - average distance between adjacent putative on-target sites (*_dist.csv)

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


def macs_gen(peak, span_r, genome, guide, mismatch=20, cent_r=200, fenr=0):
    """ Generator to yield all peaks from macs2 output.

    :param peak: macs2 output file
    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param genome: path to genome (hg38 - with .fa extension)
    :param guide: on-target protospacer sequence (no PAM)
    :param mismatch: number of mismatches from protospacer to accept
    :param cent_r: radius of window from peak center to search for cut site
    :param fenr: minimum fold enrichment of peak to be considered
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM, discovered protospacer,
            # mismatches, non-mismatched protospacer )

    """
    m_out = np.loadtxt(peak, dtype=object)    # load macs2 narrowPeak output file
    numpeaks = m_out.shape[0]                 # number of peaks from blender
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        chr_i, fenr_i = m_out[i, 0], float(m_out[i, 6])
        if chr_i in hg38size and fenr_i >= fenr:
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
                if pam_i in {'NGG', 'AGG', 'CGG', 'GGG', "TGG"}:
                    yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide


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
    """ For each target location outputted from generator, output the enrichment from a BAM file as
        (1) CSV file with each row one target location, column is enrichment values in a window
        centered at each target location
        (2) WIG file with the local (window-bounded) enrichment at each target location

    :param bamfilein: path to input BAM file
    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param fileout: path to output file name (excludes extension)

    Results in two files (WIG and CSV), described above.
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


def read_kinetics(subset_list, fileout):
    """ Given a list of processed data files that correspond to different time points, output a
        new file with the relevant data merged into one file

    :param subset_list: list of directory paths to output files of read_subsets(), where the
    column of interest is indexed 9 (starting from 0)
    :param fileout: String file path of output (no extension)
    """
    list_gen, csv_subs, num_kin, num_gen = [], [], len(subset_list), 0
    header, endindex = None, 9
    for ind, subl in enumerate(subset_list):
        if header is None:
            header = ",".join(load_npheader(subl).split(',')[:endindex])
        header += ", timepoint_%i" % (ind + 1)
        r = load_nparray(subl)
        num_gen = len(r)
        csv_subs.append(r)
    for i in range(num_gen):
        irow = csv_subs[0][i]
        gen_list = irow[:endindex].tolist()
        for j in range(num_kin):
            gen_list.append(float(csv_subs[j][i][9]))
        list_gen.append(gen_list)
    np.savetxt(fileout + ".csv", np.asarray(list_gen), fmt='%s', delimiter=',', header=header)


def read_chromhmm(generator, fileout):
    """ Determine chromatin state for each cut site in generator from obtainign the consensus
        ChromHMM annotation from different cell types.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param fileout: path to output file name (please include .csv as extension)
    """
    list_stat = []
    for rs, cut, sen, pam, gui, mis, tar in generator:
        chr_i, sta_i, end_i = c.region_string_split(rs)
        anno_i = c.get_chromhmm_annotation(chr_i, cut)
        active_i = c.get_chromhmm_active(anno_i)
        list_stat.append((rs, cut, sen, gui, mis, anno_i, active_i))
    list_stat = np.asarray(list_stat)
    if fileout:
        np.savetxt(fileout, list_stat, fmt='%s', delimiter=',')
    return list_stat


def read_subsets(generator, filein, fileout):
    """ For each target location from generator, divides the reads from filein BAM file into ones
        spanning cut site, abutting cut site, or neither as separate BAM files. Also, calculates
        total enrichment in a window centered at each target location, as well as enrichment for
        different subsets of reads, including 'left' and 'right' of the cut site in the sense
        orientation.
        Determines if each target location resides in genes - if so, determines the direction, i.e.
        'left' or 'right' of the cut site in gene orientation.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param filein: path to input BAM file for analysis
    :param fileout: path to output file name (excluding extensions)

    OUTPUT: For all target sites, within window specified by generator, outputs
    (1) fileout_span.bam: BAM file with reads that span each target site
    (2) fileout_abut.bam: BAM file with reads that abut (i.e. within 5 bp) of each target site
    (3) fileout_else.bam: BAM file with reads that neither span nor abut
    (4) fileout.csv: CSV file with information on each target site, including target sequence,
    mismatch status, total enrichment, subset enrichment (span, abut, else, abut-left, abut-right),
    orientation relative to sense strand, presence on gene, orientation relative to transcription.
    """
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
    """ For each target location from generator, determine total enrichment within target-centered
        window.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param filein: path to input BAM file for analysis
    :param fileout: path to output file name (excluding extensions)

    For all target sites, within window specified by generator, outputs CSV file with information
    on each target site, including target sequence, mismatch status, and total enrichment in window.
    """
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


def read_mismatch(generator, fileout=None):
    """ For each target location from generator, determine mismatch status, which is outputted from
        the function get_two_mismatches_loc from chipseq.py

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param fileout: path to output file name (excluding extensions)

    For all target sites, outputs CSV file with information on mismatch state.
    """
    list_stat = []
    for rs, cut, sen, pam, gui, mis, tar in generator:
        mm = c.get_two_mismatches_loc(gui, tar)
        mm.append("%02i_%s" % (mis, c.get_two_mismatches_dist(mm)))
        gen_i = [rs, cut, sen, gui, mis] + mm
        list_stat.append(gen_i)
    list_stat = np.asarray(list_stat)
    if fileout:
        np.savetxt(fileout, list_stat, fmt='%s', delimiter=',')
    return list_stat


def mergesubsetcounts(subset, countlists, num_cols, fileout, head=None):
    """ Given one output of read_subsets() and >=1 of read_counts(), merge the outputs in one file.
        For each read_counts() output, num_cols is a list that stores the # of columns starting from
        the end to take for the merged file.

        :param subset: a read_subset() output to merge. Its entirety will be added first to the
            merged CSV output file
        :param countlists: list of read_counts() outputs to merge
        :param num_cols: integer array the same length as countlists, where num_cols[i] indicates
            the number of columns (counting from the end) to take from countlists[i] to add to the
            merged output file
        :param fileout: the path of the merged CSV file (.csv extension should be included)
        :param head: (optional) the header of the merged CSV file

    """
    for x, ncol in zip(countlists, num_cols):
        subset = np.column_stack((subset, x[:, x.shape[1]-ncol:]))
    if head:
        np.savetxt(fileout, subset, fmt='%s', delimiter=',', header=head)
    else:
        np.savetxt(fileout, subset, fmt='%s', delimiter=',')
    return subset


def mergerows(files, fileout, head=None):
    """ Given a list of data files, merge the rows for each of them.

        :param files: list of files, generally from the output of read_subsets(), read_counts(), or
            mergesubsetcounts()
        :param fileout: the path of the merged CSV file (.csv extension should be included)
        :param head: (optional) the header of the merged CSV file

    """
    merged = files[0]
    for i in range(len(files)-1):
        merged = np.row_stack((merged, files[i+1]))
    if head:
        np.savetxt(fileout, merged, fmt='%s', delimiter=',', header=head)
    else:
        np.savetxt(fileout, merged, fmt='%s', delimiter=',')
    return merged

