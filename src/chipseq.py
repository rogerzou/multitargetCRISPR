# -*- coding: utf-8 -*-
""" Useful functions for ChIP-seq analysis after Cas9 cleavage
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from collections import defaultdict
import pysam
import os
import re
import numpy as np
import csv
from scipy import stats

OLD_CHAR = "ACGT"
NEW_CHAR = "TGCA"
ACTIVE = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh',
          '10_TssBiv', '11_BivFlnk', '12_EnhBiv']
REFSEQ = None
REF_INDEX = {}
REFSEQ_hg = None
e114_bed, e116_bed, e117_bed, e123_bed = None, None, None, None
e114_ind, e116_ind, e117_ind, e123_ind = {}, {}, {}, {}
CHROMHMM_hg = None


def get_reverse_complement(seq):
    """ Return reverse complement of sequence string. """
    return seq.translate(str.maketrans(OLD_CHAR, NEW_CHAR))[::-1]


def status_statement(current, final, count, chromosome=None):
    """ Print progress statements for long processes

    :param current: current iteration number
    :param final: total number of iterations
    :param count: number of print statements total
    :param chromosome: print chromosome information of current progress

    """
    if current % int(final/count) == 0:
        if chromosome is None:
            print("Processed %i out of %i" % (current, final))
        else:
            print("Processed %i out of %i in %s" % (current, final, chromosome))


def region_string_split(rs):
    return re.split('[:-]', rs)


def get_two_mismatches_loc(g_ob, g_ex):
    """ Given a observed protospacer vs expected protospacer, if 1-2 mismatches, output the position
        of mismatches starting from the PAM-distal end.

        :param g_ob: observed protospacer
        :param g_ex: expected protospacer
        :return: [position A, position B] corresponds to two mismatches at positions A and B
                 [position A, -1] corresponds to a site with 1 mismatch at position A
                 [-1, -1] corresponds to either no mismatches or >2 mismatches
    """
    wlist = [x + 1 if g_ob[j] != g_ex[j] else x for j, x in enumerate([0] * len(g_ex))]
    glen, nmm = len(g_ex), sum(wlist)
    if 1 <= nmm <= 2:
        mminds = [i.start() for i in re.finditer('1', "".join(list(map(str, wlist))))]
        if nmm == 1:
            return [glen - mminds[0], -1]
        else:
            return [glen - mminds[0], glen - mminds[1]]
    else:
        return [-1, -1]


def get_two_mismatches_dist(mm_positions):
    """ Given mismatch positions outputted from get_two_mismatches_loc(), label type of mismatch

        :param mm_positions: ooutput from get_two_mismatches_loc()
        :return: 'null' if no mismatches or > two mismatches
                 'dist' if 1-2 mismatches, all PAM-distal (positions 12-20)
                 'prox' if 1-2 mismatches, all PAM-proximal (positions 1-11)
                 'mix' if 2 mismatches, one PAM-proximal, one PAM-distal
    """
    if mm_positions[0] == -1 and mm_positions[1] == -1:        # neither 1 nor 2 mismatches
        return 'null'
    elif mm_positions[1] == -1:                            # 1 mismatch
        if mm_positions[0] >= 12:
            return 'dist'
        elif 0 <= mm_positions[0] < 12:
            return 'prox'
        else:
            return ValueError("Error")
    else:                                                   # 2 mismatches
        if mm_positions[0] >= 12 and mm_positions[1] >= 12:
            return 'dist'
        elif 0 <= mm_positions[0] < 12 and 0 <= mm_positions[1] < 12:
            return 'prox'
        else:
            return 'mix'


def is_chromhmm(hg):
    global CHROMHMM_hg
    if any(v is None for v in [e114_bed, e116_bed, e117_bed, e123_bed]) or CHROMHMM_hg != hg:
        return False
    return True


def chromhmm_initialize(hg):
    """ Initialize chromosomal indexing of chromhmm bed file for rapid querying.
    :param hg: 'hg19' or 'hg38' """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    global e114_bed, e116_bed, e117_bed, e123_bed, \
        e114_ind, e116_ind, e117_ind, e123_ind, CHROMHMM_hg
    CHROMHMM_hg = hg
    e114_bed, e114_ind = bed_indexing(os.path.dirname(dirname) +
                                      "/lib/chromhmm/E114_15_coreMarks_%s_dense.bed" % hg)
    e116_bed, e116_ind = bed_indexing(os.path.dirname(dirname) +
                                      "/lib/chromhmm/E116_15_coreMarks_%s_dense.bed" % hg)
    e117_bed, e117_ind = bed_indexing(os.path.dirname(dirname) +
                                      "/lib/chromhmm/E117_15_coreMarks_%s_dense.bed" % hg)
    e123_bed, e123_ind = bed_indexing(os.path.dirname(dirname) +
                                      "/lib/chromhmm/E123_15_coreMarks_%s_dense.bed" % hg)


def get_chromhmm_annotation(hg, chromosome, coordinate):
    """ Determine the consensus (mode) chromatin state from ChromHMM of 4 cell lines.

    :param hg: 'hg19' or 'hg38'
    :param chromosome: chromosome of location
    :param coordinate: coordinate of location
    :return: Name and direction (sense) of gene if it exists as a tuple, None if doesn't exist
    """
    if not is_chromhmm(hg):
        chromhmm_initialize(hg)
    e114_row = bed_getrow(e114_bed, e114_ind, chromosome, coordinate)
    e116_row = bed_getrow(e116_bed, e116_ind, chromosome, coordinate)
    e117_row = bed_getrow(e117_bed, e117_ind, chromosome, coordinate)
    e123_row = bed_getrow(e123_bed, e123_ind, chromosome, coordinate)
    lst = _get_chromhmm_annotation_helper(e114_row, e116_row, e117_row, e123_row)
    if len(lst) > 0:
        return max(set(lst), key=lst.count)
    else:
        return None


def _get_chromhmm_annotation_helper(e114, e116, e117, e123):
    lst = [None, None, None, None]
    if e114 is not None:
        lst[0] = e114[3]
    if e116 is not None:
        lst[1] = e116[3]
    if e117 is not None:
        lst[2] = e117[3]
    if e123 is not None:
        lst[3] = e123[3]
    return [x for x in lst if x is not None]


def get_chromhmm_active(annotation):
    return 'active' if annotation in ACTIVE else 'inactive'


def bed_indexing(bedfile):
    """ Index bed file by chromosome - facilitates search only within a specific chromosome.

    :return bed: numpy array holding generic BED file
    :return bed_index: indexing dict with keys as chromosomes, values as start and end indexes from
                       bed file.
    """
    bed = np.loadtxt(bedfile, dtype=object)
    bed_index = {}
    chr_prev = None
    sta_prev = None
    for i in range(bed.shape[0]):
        chr_i = bed[i, 0]
        if chr_prev != chr_i:
            if i != 0:
                bed_index[chr_prev] = (sta_prev, i-1)
            chr_prev = chr_i
            sta_prev = i
        if i == bed.shape[0]-1:
            bed_index[chr_i] = (sta_prev, i)
    return bed, bed_index


def bed_getrow(bed, bed_index, chromosome, coordinate):
    """ Output specific row of bed file that corresponds to the correct chromosome and coordinate.

    :param bed: bed file
    :param bed_index: index file
    :param chromosome: chromosome of location
    :param coordinate: coordinate of location
    :return: Row of bed file
    """
    if chromosome in bed_index:
        chr_sta, chr_end = bed_index[chromosome]
        for i in range(chr_sta, chr_end+1):
            row_i = bed[i, :]
            row_chr = row_i[0]
            row_sta = int(row_i[1])
            row_end = int(row_i[2])
            if chromosome != row_chr:
                raise ValueError("bed_getrow(): Indexing of bedfile has errors.")
            if row_sta <= coordinate <= row_end:
                return row_i
    else:
        print("bed_getrow(): queried chromosome %s is not in bedfile." % chromosome)
    return None


def refseq_initialize(hg):
    """ Initialize global database of RefSeq gene annotations with indices for fast querying

    :param hg: 'hg38' or 'hg19'
    REFSEQ: global numpy array holding BED file of RefSeq gene annotations (hg38)
    REF_INDEX: indexing dict with keys as chromosomes, values as start and end indexes from REFSEQ
               numpy matrix for each chromosome key
    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    d = np.loadtxt(os.path.dirname(dirname) + "/lib/refseq_%s.bed" % hg, dtype=object)
    global REFSEQ, REF_INDEX, REFSEQ_hg
    REFSEQ_hg = hg
    REFSEQ = d[:, [0, 1, 2, 3, 5]]
    chr_prev = None
    sta_prev = None
    for i in range(REFSEQ.shape[0]):
        chr_i = REFSEQ[i, 0]
        if chr_prev != chr_i:
            if i != 0:
                REF_INDEX[chr_prev] = (sta_prev, i-1)
            chr_prev = chr_i
            sta_prev = i
        if i == REFSEQ.shape[0]-1:
            REF_INDEX[chr_i] = (sta_prev, i)


def is_refseq(hg):
    global REFSEQ, REFSEQ_hg
    if REFSEQ is not None and REFSEQ_hg == hg:
        return True
    return False


def is_gene_refseq(hg, chromosome, coordinate):
    """ Check if a specific chromosomal location resides in an annotated RefSeq gene.

    :param hg: either 'hg38' or 'hg19'
    :param chromosome: chromosome of location
    :param coordinate: coordinate of location
    :return: Name and direction (sense) of gene if it exists as a tuple, None if doesn't exist
    """
    global REFSEQ, REF_INDEX
    if not is_refseq(hg):   # if RefSeq is not initialized for a specific genome, initialize again
        refseq_initialize(hg)
    if chromosome in REF_INDEX:
        chr_sta, chr_end = REF_INDEX[chromosome]
        for i in range(chr_sta, chr_end+1):
            gene_i = REFSEQ[i, :]
            gene_chr = gene_i[0]
            gene_sta = int(gene_i[1])
            gene_end = int(gene_i[2])
            if chromosome != gene_chr:
                raise ValueError("is_gene_refseq(): Indexing of REFSEQ has errors.")
            if gene_sta <= coordinate <= gene_end:
                return gene_i[3], gene_i[4]
    else:
        print("is_gene_refseq(): queried chromosome %s is not in REFSEQ." % chromosome)
    return None


def hg_dict(genome_str):
    """ Return dict that holds the number of base pairs for each chromosome in hg38 (human)

    :param genome_str: either 'hg19' or 'hg38'
    :return: dict with keys as chromosomes, values as maximum coordinate of each chromosome key
    """
    d = {}
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s.sizes" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            d[row[0]] = int(row[1])
    return d


def hg_generator(genome_str):
    """ Generate the number of base pairs for each chromosome in hg38 (human)

    :param genome_str: either 'hg19' or 'hg38'
    :return generator that outputs the number of base pairs for each chromosome in hg38
            in the format: [chr7, 159345973]
    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/%s.sizes" % genome_str, 'r') as f:
        for row in csv.reader(f, delimiter='\t'):
            yield row


def read_pair_generator(bam, region_string=None):
    """ Generate read pairs in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.

    :param bam: pysam AlignmentFile loaded with BAM file that contains paired-end reads
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :return generator of read pairs with the following format: [read1, read2]

    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_pair_align(read1, read2):
    """ Extract read pair locations as a fragment oriented in increasing chromosome coordinates

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return 4-item array in the following format: [fragA-start, fragA-end, fragB-start, fragB-end]
            with monotonically increasing chromosome coordinates

    """
    r1pos = [x+1 for x in read1.positions]
    r2pos = [x+1 for x in read2.positions]
    if read1.mate_is_reverse and r1pos[0] < r2pos[0]:  # read1 is earlier
        read = [r1pos[0], r1pos[-1], r2pos[0], r2pos[-1]]
    elif read2.mate_is_reverse and r2pos[0] < r1pos[0]:  # read2 is earlier
        read = [r2pos[0], r2pos[-1], r1pos[0], r1pos[-1]]
    else:
        read = []
        # print("Skipping read pair from error in alignment.")
    # print("%s--%s>  <%s--%s" % tuple(read))
    return read


def get_read_subsets(filein, fileout, region_string, target, bool_array=None):
    """ Categorize fragments (paired reads) by its relationship to the target (cut) site

    For a target (cut) site, subset the fragments within region_range to
    - fragments that span the cut site
    - fragments that do not -> fragments that start from the cut site and span left vs right

    :param filein: BAM file containing paired-end reads
    :param fileout: base output file name with extension (.bam) omitted - function generates
                    multiple files where each file only contains reads from a specific category
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :param target: coordinate within region_string that denotes the 4th nucleotide from PAM
    :param bool_array: boolean array indicating whether to output the subsetted reads in the
                       following order: [M, N, L, R, S1, S2, I]

    """

    if bool_array is None:
        bool_array = [True, True, False, False]
    bam = pysam.AlignmentFile(filein, 'rb')
    countM = 0      # fragments that span
    countN = 0      # fragments that don't span
    countL = 0      # fragments that don't span, on left
    countR = 0      # fragments that don't span, on right
    counter = 0     # count total number of reads
    fileM = fileout + "_M.bam"
    fileN = fileout + "_N.bam"
    fileL = fileout + "_L.bam"
    fileR = fileout + "_R.bam"
    readsM = pysam.AlignmentFile(fileM, "wb", template=bam)
    readsN = pysam.AlignmentFile(fileN, "wb", template=bam)
    readsL = pysam.AlignmentFile(fileL, "wb", template=bam)
    readsR = pysam.AlignmentFile(fileR, "wb", template=bam)
    for read1, read2 in read_pair_generator(bam, region_string):
        read = read_pair_align(read1, read2)
        if not read:
            continue
        counter += 1
        if read[0] < target < read[3]:          # fragments that span
            countM += 1
            readsM.write(read1)
            readsM.write(read2)
        elif target + 5 >= read[0] >= target:   # fragments that begin 5bp of cleavage site
            countR += 1
            readsR.write(read1)
            readsR.write(read2)
            countN += 1
            readsN.write(read1)
            readsN.write(read2)
        elif target - 5 <= read[-1] <= target:   # fragments that end 5bp of cleavage site
            countL += 1
            readsL.write(read1)
            readsL.write(read2)
            countN += 1
            readsN.write(read1)
            readsN.write(read2)
    readsM.close()
    readsN.close()
    readsL.close()
    readsR.close()
    file_array = [fileM, fileN, fileL, fileR]
    for i, boo in enumerate(bool_array):
        file_i = file_array[i]
        if boo:
            pysam.sort("-o", file_i, file_i)
            os.system("samtools index " + file_i)
        else:
            os.system("rm " + file_i)
    print("%i span | %i end/start | %i window total | %i mapped total"
          % (countM, countN, counter, bam.mapped))
    bam.close()


def to_wiggle_pairs(filein, fileout, region_string, endcrop=False):
    """ Constructs fragment pile-ups in wiggle format using paired-end information

    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param region_string: region of interest, formatted like this example: chr7:5527160-5532160
    :param endcrop: if True, don't count the first and last read of a fragment - the output for each
                    position is effectively a safe underestimate for the number of spanning reads

    """
    [chr_i, sta, end] = region_string_split(region_string)
    sta = int(sta)
    end = int(end)
    wlist = [0] * (end-sta+1)
    wig = open(fileout + ".wig", "w")
    wig.write("variableStep\tchrom=%s\n" % chr_i)
    bam = pysam.AlignmentFile(filein, 'rb')
    for read1, read2 in read_pair_generator(bam, region_string):
        read = read_pair_align(read1, read2)
        if not read:
            continue
        if endcrop:
            wlist = [x+1 if read[0]-sta < i < read[-1]-sta else x for i, x in enumerate(wlist)]
        else:
            wlist = [x+1 if read[0]-sta <= i <= read[-1]-sta else x for i, x in enumerate(wlist)]
    for i, x in enumerate(wlist):
        wig.write("%i\t%i\n" % (sta + i, x))
    wig.close()
    bam.close()


def to_wiggle_windows(genome, filein, fileout, window, chromosome=None, generator=None):
    """ Outputs wiggle file that counts reads per million in each window

    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param filein: BAM file that contains paired-end reads
    :param fileout: base output file name with extension (.wig) omitted
    :param window: size of window (i.e. if 500, genome is divided into 500bp windows, function
                outputs the number of reads in each window)
    :param chromosome: array of chromosome strings to limit analysis to particular chromosomes
                       (i.e. ['chr7', 'chr8']). If not set, then all chromosomes are processed.
    :param generator: species-specific generator that outputs the number of base pairs for each
                      chromosome in the format: [chr7, 159345973]. If not set, then hg38 is used.

    """
    if not generator:
        generator = hg_generator(genome[0])
    bam = pysam.AlignmentFile(filein, 'rb')
    treads = float(bam.mapped)
    wig = open(fileout + ".wig", "w")
    for row in generator:  # iterate over each chromosome
        if chromosome is None or (chromosome is not None and row[0] in chromosome):
            chr_i = row[0]
            wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
            numwins = int(int(row[1]) / window)  # number of windows
            for i in range(numwins):              # iterate over each bin (total - 1)
                start1 = i * window
                finish1 = (i + 1) * window - 1
                rpm = bam.count(chr_i, start1, finish1) / treads * 1E6
                wig.write("%0.6f\n" % rpm)
                status_statement(i, numwins, 20, chr_i)
    wig.close()
    bam.close()


def ttest_two(samp_file, ctrl_file, fileout, p=0.01):
    """ From to_bins() outputs, perform significance testing between the bins of each window,
        between experimental and control conditions

    - T-test is passed if bin values of experimental condition windows is significantly higher than
    those of the control condition (one-sided with Bonferroni correction)
    - Outputs CSV file of binning results with t-test result
    - Outputs WIG file indicating windows with significant t-test results

    :param samp_file: output of to_bins() as CSV file for experimental condition
    :param ctrl_file: output of to_bins() as CSV file for control as comparison
    :param fileout: base output file name with extension (.csv) omitted; each row corresponds
                    to each window, each column corresponds to the number of reads for each bin
    :param p: p value cut-off

    """
    np_samp = np.loadtxt(samp_file, delimiter=',', dtype=object)
    np_ctrl = np.loadtxt(ctrl_file, delimiter=',', dtype=object)
    count = np_samp.shape[0]
    if np_samp.shape[0] != np_ctrl.shape[0]:
        raise ValueError("Number of windows between sample and control files are different!")
    ot = np.zeros((count, 9), dtype=object)         # array to hold info for each bin
    wig = open(fileout + "_ttest.wig", "w")
    chr_i = None
    prev = 0
    chr_ctrl = np_ctrl[:, 0]                        # chr
    chr_samp = np_samp[:, 0]
    ran_ctrl = np_ctrl[:, 1:3].astype(int)          # start and end coordinates
    ran_samp = np_samp[:, 1:3].astype(int)
    bin_ctrl = np_ctrl[:, 3:].astype(int)           # bins
    bin_samp = np_samp[:, 3:].astype(int)
    for i in range(count):
        if chr_samp[i] != chr_ctrl[i] or ran_samp[i, 0] != ran_ctrl[i, 0] \
                or ran_samp[i, 1] != ran_ctrl[i, 1]:
            raise ValueError("Between sample and control files, the chr, start, or end coordinates "
                             "are different!")
        else:
            if chr_samp[i] != chr_i:
                chr_i = chr_samp[i]
                wig.write("variableStep\tchrom=%s\n" % chr_i)
            ot[i, 0:3] = np_samp[i, 0:3]                    # chr, start, and end coordinates
            ot[i, 3] = np.sum(bin_ctrl[i, :])               # sum of ctrl
            ot[i, 4] = np.sum(bin_samp[i, :])               # sum of sample
            ttest = stats.ttest_rel(bin_samp[i, :], bin_ctrl[i, :])
            ot[i, 5] = ttest[0]                             # t test
            ot[i, 6] = ttest[1]                             # t test
            start_i = ran_samp[i, 0]
            if ot[i, 5] > 0 and ot[i, 6] / 2 < p / count:   # one-sided t-test with Bonferroni
                ot[i, 7] = 1
                ot[i, 8] = start_i - prev
                prev = start_i
                wig.write("%i\t%i\n" % (start_i, 1))
            else:
                ot[i, 7] = np.nan
                ot[i, 8] = np.nan
                wig.write("%i\t%i\n" % (start_i, 0))
            status_statement(i, count, 20, chr_i)
    wig.close()
    np.savetxt(fileout + "_ttest.csv", ot, fmt='%s', delimiter=',')


def ttest_span(samp_file, fileout, chrs, cuts, names, skipl):
    """ Calculate the width of broad peaks, given center positions and a constant skip length, from
        ttest of binning analysis; outputs in broadPeak format.

    :param samp_file: input file that is the output of "ttest_two()"
    :param fileout: output file (extension omitted)
    :param chrs: array of chromosomes that correspond to each cut site
    :param cuts: array of coordinates that denote cut site
    :param names: array of names for each cut site
    :param skipl: skip length

    """
    np_samp = np.loadtxt(samp_file, delimiter=',', dtype=object)
    f = open(fileout + ".broadPeak", 'w')
    for chr_i, cut_i, name_i in zip(chrs, cuts, names):     # iterate over each expected peak
        cutrow = None                   # row of np_samp with expected cut site
        minrow = None                   # row of np_samp that is the furthest upstream of cut site
        maxrow = None                   # row of np_samp that is the furthest downstream of cut site
        for i in range(np_samp.shape[0]):           # assign values to minrow, maxrow, cutrow
            if np_samp[i, 7] != "nan":
                if minrow is None and np_samp[i, 0] == chr_i:
                    minrow = i
                    maxrow = i
                if maxrow is not None and np_samp[i, 0] == chr_i:
                    maxrow = i
                if np_samp[i, 0] == chr_i and int(np_samp[i, 1]) <= cut_i <= int(np_samp[i, 2]):
                    cutrow = i
        if cutrow is not None:                      # determine peak span if peak is found
            # check in the upstream (rev) direction
            revrow_tmp = cutrow
            revrow = cutrow
            revct = 0
            while revrow_tmp > minrow and revct <= skipl:
                revrow_tmp -= 1
                if np_samp[revrow_tmp, 7] == 'nan':
                    revct += int(np_samp[revrow_tmp, 2]) - int(np_samp[revrow_tmp, 1]) + 1
                else:
                    revct = 0
                    revrow = revrow_tmp
            fwdrow_tmp = cutrow
            fwdrow = cutrow
            fwdct = 0
            # check in the downstream (fwd) direction
            while fwdrow_tmp < maxrow and fwdct <= skipl:
                fwdrow_tmp += 1
                if np_samp[fwdrow_tmp, 7] == 'nan':
                    fwdct += int(np_samp[fwdrow_tmp, 2]) - int(np_samp[fwdrow_tmp, 1]) + 1
                else:
                    fwdct = 0
                    fwdrow = fwdrow_tmp
            # get coordinates in both upstream and downstream direction, write to file
            revval = int(np_samp[revrow, 1])
            fwdval = int(np_samp[fwdrow, 2]) + 1
            width = fwdval - revval
            f.write("%s\t%i\t%i\t%s\t1\t.\t%i\t-1\t-1\n" % (chr_i, revval, fwdval, name_i, width))
        else:
            print("Cut site has no enrichment")
    f.close()


def avgwig(file1, file2, fileout):
    """ Given two wig files of same structure, output the average of both as a wig file

    :param file1: first wig file
    :param file2: second wig file
    :param fileout: wig file (no extension) that outputs the average of the two input files

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".wig", "w") as outwig:
        for line1, line2 in zip(f1, f2):
            if line1.split("\t")[0] == "fixedStep" or line2.split("\t")[0] == "fixedStep":
                if line1 == line2:
                    outwig.write(line1)
                else:
                    raise ValueError("The files are too different to merge!")
            else:
                val1 = float(line1.split("\n")[0])
                val2 = float(line2.split("\n")[0])
                outwig.write("%0.5f\n" % np.mean([val2, val1]))


def avgspan(file1, file2, fileout):
    """ Given two broadPeak files listing peak spans, output the average as a broadPeak file

    :param file1: broadPeak file 1 (e.x. replicate 1)
    :param file2: broadPeak file 2 (e.x. replicate 2)
    :param fileout: output broadPeak file (extension omitted)

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".broadPeak", "w") as outwg:
        for line1, line2 in zip(f1, f2):
            spl1 = line1.split('\t')
            spl2 = line2.split('\t')
            outline = spl1
            outline[1] = str(int(np.mean([int(spl1[1]), int(spl2[1])])))
            outline[2] = str(int(np.mean([int(spl1[2]), int(spl2[2])])))
            outline[6] = str(int(outline[2]) - int(outline[1]) + 1)
            outline[8] = outline[8].split('\n')[0]
            outwg.write('\t'.join(outline) + '\n')


def percentchange(file1, file2, fileout, cutoff=0.2):
    """ Given two wig files, calculate the percent change from file1 to file2

    :param file1: first wig file
    :param file2: second wig file
    :param fileout: output wig file (extension omitted)
    :param cutoff: only values greater than cutoff will be analyzed, to prevent dividing by small
                   values, which would be numerically unstable

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".wig", "w") as outwig:
        for l1, l2 in zip(f1, f2):
            if l1.split("\t")[0] == "fixedStep" or l2.split("\t")[0] == "fixedStep":
                if l1 == l2:
                    outwig.write(l1)
                else:
                    raise ValueError("The files are too different to compare!")
            else:
                val1 = float(l1.split("\n")[0])
                val2 = float(l2.split("\n")[0])
                if val1 <= cutoff:
                    outwig.write("%0.5f\n" % 0)
                else:
                    outwig.write("%0.5f\n" % ((val2 - val1) / val1))


def absolutechange(file1, file2, fileout):
    """ Given two wig files, calculate the absolute change from file1 to file2

    :param file1: first wig file
    :param file2: second wig file
    :param fileout: output wig file (extension omitted

    """
    with open(file1, "r") as f1, open(file2, "r") as f2, open(fileout + ".wig", "w") as outwig:
        for l1, l2 in zip(f1, f2):
            if l1.split("\t")[0] == "fixedStep" or l2.split("\t")[0] == "fixedStep":
                if l1 == l2:
                    outwig.write(l1)
                else:
                    raise ValueError("The files are too different to compare!")
            else:
                val1 = float(l1.split("\n")[0])
                val2 = float(l2.split("\n")[0])
                outwig.write("%0.5f\n" % (val2 - val1))
