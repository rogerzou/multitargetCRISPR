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
import re
import numpy as np
from . import chipseq as c

WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3]


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
        return c.get_reverse_complement(read1.seq), read1.qual[::-1], read2.seq, read2.qual
    elif read2.is_reverse and not read1.is_reverse:
        return read1.seq, read1.qual, c.get_reverse_complement(read2.seq), read2.qual[::-1]
    else:
        raise ValueError("_find_msa_helper(): Unexpected paired-end read conditions.")


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
        substr_0 = c.get_reverse_complement(s[i:i + sublen])
        pam_0 = c.get_reverse_complement(s[i - 3:i])
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


def peak_profile_wide(generator, bamfilein, fileout, span_rad=2000000, res=100, wind_rad=10000):
    """ For each target location from generator, calculates enrichment at specified 'resolution'
        with sliding window of specified 'radius'. Outputs the enrichment from a BAM file as:
        (1) CSV file with each row one target location, column is enrichment values in a window
        centered at each target location, at a specific genomic resolution.
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
    :param span_rad:
    :param res:
    :param wind_rad:

    Results in two files (WIG and CSV), described above.
    """
    hg38size = c.hg38_dict()
    bamin = pysam.AlignmentFile(bamfilein, 'rb')
    chr_old, csv_peaks = None, []
    wlist_all = []
    numrows = int(span_rad * 2 / res) + 1
    for rs, cut, sen, pam, gui, mis, guide in generator:
        chr_i = re.split('[:-]', rs)[0]
        sta_i = cut - span_rad
        end_i = cut + span_rad
        if sta_i - wind_rad >= 0 and end_i + wind_rad < hg38size[chr_i]:
            wlist = [0] * numrows
            for row_i in range(numrows):
                center = sta_i + row_i * res
                rs_i = "%s:%i-%i" % (chr_i, center - wind_rad, center + wind_rad)
                wlist[row_i] = bamin.count(region=rs_i) / bamin.mapped * 1E6
            wlist_all.append([chr_i, sta_i] + wlist)
            csv_peaks.append([chr_i, sta_i, gui + pam, mis] + wlist)
    bamin.close()
    _peak_profile_helper(wlist_all, res, fileout)
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',')


def peak_profile_bp_resolution(generator, bamfilein, fileout):
    """ For each target location from generator, calculates enrichment at each base pair as the
        number of fragments that 'span' the base, i.e. the base is either (1) sequenced by either
        end of paired-end sequencing, or (2) not sequenced but spanned by the imputed DNA fragment.
        Output the enrichment from a BAM file as:
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
        wlist = wlist if sen == '+' else wlist[::-1]
        csv_peaks.append([chr_i, sta_i, gui + pam, mis] + wlist)
    bamin.close()
    _peak_profile_helper(wlist_all, 1, fileout)
    np.savetxt(fileout + "_bpeaks.csv", np.asarray(csv_peaks), fmt='%s', delimiter=',')


def _peak_profile_helper(wlist_all, resolution, fileout):
    """ Helper function for peak_profile_bp_resolution() or peak_profile_wide(). Writes all peak
        profiles to a wiggle file.

    :param wlist_all:
    :param resolution:
    :param fileout:
    """
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
                wigout.write("%i\t%0.5f\n" % (sta_i + j * resolution, x))


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

