# miscellaneous functions for ChIP-seq analysis after Cas9 cleavage

from collections import defaultdict
import pysam
import os
import re
import numpy as np
import csv
from scipy import stats


def bam_pairs_to_wiggle(infile, outfile, region_range):
    """
    Converts paired-end reads as sam/bam files from infile to bigwig outfile within
    specific region_range
    """
    [chr, sta, end] = re.split('[:-]', region_range)
    sta = int(sta)
    end = int(end)
    wlist = [0] * (end-sta+1)
    f = open(outfile, "w")
    f.write("variableStep\tchrom=%s\n" % chr)
    bam = pysam.AlignmentFile(infile, 'rb')
    for read1, read2 in read_pair_generator(bam, region_range):
        read = read_pair_align(read1, read2)
        wlist = [x+1 if read[0]-sta <= i <= read[-1]-sta else x for i, x in enumerate(wlist)]
        # print("%i\t%i" % (read[0]-sta, read[-1]-sta))
    for i, x in enumerate(wlist):
        f.write("%i\t%i\n" % (sta + i, x))
    f.close()
    bam.close()


def tobigwig(base, fname, label, window):
    bam1 = pysam.AlignmentFile(base + fname, 'rb')
    with open(base + "hg38.sizes", 'r') as f:   # parse sizes for each chromosome of hg38
        counter = 0
        outpath = base + "analysis/" + "%s_%i" % (label, window)
        wig = open(outpath + "_subbin.wig", "w")
        for row in csv.reader(f, dialect='excel', delimiter='\t'):  # iterate over each chromosome
            if counter < 24:
                print("Counter: %i" % counter)
                counter += 1
                chr_i = row[0]
                wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
                count = int(int(row[1]) / window)  # number of bins
                for i in range(count):              # iterate over each bin (total - 1)
                    start1 = i * window
                    finish1 = (i + 1) * window - 1
                    wig.write("%i\n" % bam1.count(chr_i, start1, finish1))
        wig.close()


def subbinsave(base, fname, label, window, res, chr=None):
    bam1 = pysam.AlignmentFile(base + fname, 'rb')
    with open(base + "hg38.sizes", 'r') as f:   # parse sizes for each chromosome of hg38
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        counter = 0
        for row in reader:  # iterate over each chromosome
            if chr is None or (chr is not None and row[0] in chr):  # checks for chr #
                outpath = base + "analysis/" + "%s_%i_%i" % (label, window, res)
                range_j = window / res  # number of sub-bins per bin for t-test
                outfile = np.zeros((1, 3 + range_j), dtype=object)  # arrays to hold sub-bin info
                wig = open(outpath + "_subbin.wig", "w")
                if counter < 24:
                    counter += 1
                    count = int(int(row[1]) / window)       # number of bins
                    chr_i = row[0]
                    wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
                    outfile_i = np.zeros((count, 3 + range_j), dtype=object)  # arrays to hold sub-bin info
                    for i in range(count):              # iterate over each bin (total - 1)
                        start1 = i * window
                        finish1 = (i + 1) * window - 1
                        for j in range(range_j):        # get all sub-bin counts for each bin
                            outfile_i[i, 0] = chr_i
                            outfile_i[i, 1] = start1
                            outfile_i[i, 2] = finish1
                            outfile_i[i, j+3] = bam1.count(chr_i, start1+j*res, start1+(j+1)*res-1)
                        if i % 100 == 0:
                            print("Making sub-bin file | %s: %i of %i" % (chr_i, i, count))
                        np.append(outfile, outfile_i, axis=0)
                        wig.write("%i\n" % bam1.count(chr_i, start1, finish1))
                np.savetxt(outpath + "_subbin.csv", outfile, delimiter=',')
                wig.close()
    return outpath


def ttest_one(subbinfile, window, p=0.01):
    """
    Bins ONE file by window, performs significance testing between consecutive bins using
    sub-bins by res
    Outputs csv file of binning results with t-test results
    Outputs bigwig file of binning results
    """
    outpath = re.split('[.]', subbinfile)[0]
    subbins = np.loadtxt(subbinfile, delimiter=',')
    count = subbins.shape[0]
    ot = np.zeros((count, 7), dtype=float)  # array to hold info for each bin
    for i in range(count):
        if 0 < i < count-1:
            start1 = i * window
            ot[i, 0] = start1                                       # start coordinate
            ot[i, 1] = np.sum(subbins[i, :])
            ttest1 = stats.ttest_ind(subbins[i, :], subbins[i - 1, :])
            ttest2 = stats.ttest_ind(subbins[i, :], subbins[i + 1, :])
            ot[i, 2] = ttest1[0]   # t test
            ot[i, 3] = ttest1[1]   # t test
            ot[i, 4] = ttest2[0]   # t test
            ot[i, 5] = ttest2[1]   # t test
            ot[i, 6] = 1 if (ot[i, 2] > 0 and ot[i, 3]/2 < p / count) or \
                            (ot[i, 4] > 0 and ot[i, 5]/2 < p / count) else np.nan   # significance
            if i % 100 == 0:
                print("Parsing sub-bins | %i of %i" % (i, count))
    np.savetxt(outpath + "_ttest-one.csv", ot, delimiter=',')


def ttest_two(filesamp, filectrl, window, p=0.01):
    """
    Bins ONE file by window, performs significance testing between consecutive bins using
    sub-bins by res
    Outputs csv file of binning results with t-test results
    Outputs bigwig file of binning results
    """
    outpath = re.split('[.]', filesamp)[0]
    subbinssamp = np.loadtxt(filesamp, delimiter=',')
    subbinsctrl = np.loadtxt(filectrl, delimiter=',')
    count = subbinssamp.shape[0]
    ot = np.zeros((count, 7), dtype=float)  # array to hold info for each bin
    prev = 0
    for i in range(count):
            start1 = i * window
            ot[i, 0] = start1                                       # start coordinate
            ot[i, 1] = np.sum(subbinsctrl[i, :])
            ot[i, 2] = np.sum(subbinssamp[i, :])
            ttest = stats.ttest_rel(subbinssamp[i, :], subbinsctrl[i, :])
            ot[i, 3] = ttest[0]   # t test
            ot[i, 4] = ttest[1]   # t test
            if ot[i, 3] > 0 and ot[i, 4] / 2 < p / count:
                ot[i, 5] = 1
                ot[i, 6] = start1 - prev
                prev = start1
            else:
                ot[i, 5] = np.nan
                ot[i, 6] = np.nan
            if i % 100 == 0:
                print("Parsing sub-bins | %i of %i" % (i, count))
    np.savetxt(outpath + "_ttest-two.csv", ot, delimiter=',')


# assuming hg38
def binning2ttest(base, fname1, fname2, label, window, res, chr=None, p=0.01):
    """
    Bins TWO files by window, performs significance testing from sub-binning by res
    Outputs csv file of binning results with t-test results
    Outputs bigwig file of binning results
    """
    bam1 = pysam.AlignmentFile(base + fname1, 'rb')
    bam2 = pysam.AlignmentFile(base + fname2, 'rb')
    with open(base + "hg38.sizes", 'r') as f:   # parse sizes for each chromosome of hg38
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        counter = 0
        for row in reader:  # iterate over each chromosome
            if chr is None or (chr is not None and row[0] in chr):  # checks for chr #
                if counter < 24:
                    counter += 1
                    count = int(int(row[1]) / window)       # number of bins
                    ot = np.zeros((count, 6), dtype=float)  # array to hold info for each bin
                    chr_i = row[0]
                    outpath = base + "analysis/" + label + "_" + chr_i
                    wig1 = open(outpath + "_1.wig", "w")
                    wig1.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
                    wig2 = open(outpath + "_2.wig", "w")
                    wig2.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
                    for i in range(count):                  # iterate over each bin
                        range_j = window / res              # number of sub-bins per bin for t-test
                        outfile1 = np.zeros((range_j,), dtype=int)  # arrays to hold sub-bin info
                        outfile2 = np.zeros((range_j,), dtype=int)
                        start = i * window
                        finish = (i + 1) * window - 1
                        for j in range(range_j):                    # get all sub-bin counts for each bin
                            outfile1[j] = bam1.count(chr_i, start + j * res,
                                                     start + (j + 1) * res - 1)
                            outfile2[j] = bam2.count(chr_i, start + j * res,
                                                     start + (j + 1) * res - 1)
                        ot[i, 0] = start                                    # start coordinate
                        ot[i, 1] = bam1.count(chr_i, start, finish)         # number of reads for file1
                        ot[i, 2] = bam2.count(chr_i, start, finish)         # number of reads for file2
                        ot[i, 3] = stats.ttest_ind(outfile1, outfile2)[0]   # t test
                        ot[i, 4] = stats.ttest_ind(outfile1, outfile2)[1]   # t test
                        ot[i, 5] = 1 if ot[i, 3] > 0 and ot[i, 4] / 2 < p / count else np.nan    # significance
                        if i % 100 == 0:
                            print("%s: %i of %i" % (chr_i, i, count))
                        wig1.write("%i\n" % (ot[i, 1]))
                        wig2.write("%i\n" % (ot[i, 2]))
                    np.savetxt(outpath + ".csv", ot, delimiter=',')
                else:
                    break


def binning(base, fname, label, window, chr=None):
    """
    """
    bam = pysam.AlignmentFile(base + fname, 'rb')
    with open(base + "hg38.sizes", 'r') as f:   # parse sizes for each chromosome of hg38
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        counter = 0
        for row in reader:  # iterate over each chromosome
            if chr is None or (chr is not None and row[0] in chr):  # checks for chr #
                if counter < 24:
                    counter += 1
                    count = int(int(row[1]) / window)       # number of bins
                    ot = np.zeros((count, 6), dtype=float)  # array to hold info for each bin
                    chr_i = row[0]
                    outpath = base + "analysis/" + label + "_" + chr_i
                    wig = open(outpath + ".wig", "w")
                    wig.write("fixedStep\tchrom=%s\tstart=0 step=%i\n" % (chr_i, window))
                    for i in range(count):                  # iterate over each bin
                        start = i * window
                        finish = (i + 1) * window - 1
                        ot[i, 0] = start                                    # start coordinate
                        ot[i, 1] = bam.count(chr_i, start, finish)         # number of reads for file1
                        if i % 100 == 0:
                            print("%s: %i of %i" % (chr_i, i, count))
                        wig.write("%i\n" % (ot[i, 1]))
                    np.savetxt(outpath + ".csv", ot, delimiter=',')
                else:
                    break
