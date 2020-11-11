# -*- coding: utf-8 -*-
""" Hi-C analysis of ChIP-seq after multi-targeting Cas9
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

import h5py
import pysam
import os
import re
import numpy as np
from scipy import sparse
from . import chipseq as c


def get_span_width(generator, genome, f_test, f_ctrl, outpath, w_rad=10000, skip=5000, false_ct=10):
    """ Determine width of 53BP1 or gH2AX enrichment by comparing test sample to negative control
        sample. Extending from the cut site at fixed intervals, enrichment width on either end of
        the cut site is defined to be where there are under 'false_ct' evaluations of negative
        control sample enrichment that is higher than test sample enrichment.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg38', path/to/hg38.fa]
    :param f_test: test sample BAM file
    :param f_ctrl: negative control BAM file
    :param outpath: path to output BED file (.bed extension omitted)
    :param w_rad: radius of window of enrichment evaluation at each site
    :param skip: number of bases to skip per evaluation of enrichment over control
    :param false_ct: maximum number of times control sample has higher enrichment than test sample
                     for a region to be included in enrichment span width centered at the cut site
    """
    hgsize = c.hg_dict(genome[0])
    outbed = open(outpath + ".bed", 'w')
    outnpy = []
    bam_test, bam_ctrl = pysam.AlignmentFile(f_test, 'rb'), pysam.AlignmentFile(f_ctrl, 'rb')
    for rs, cut, sen, pam, gui, mis, guide in generator:
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        index_neg, count_neg, width_neg = 0, 0, 0
        while True:
            index_neg -= skip
            ind_lt_neg, ind_rt_neg = cut + index_neg - w_rad, cut + index_neg + w_rad
            if ind_lt_neg >= 0:
                rs_neg = chr_i + ":" + str(ind_lt_neg) + "-" + str(ind_rt_neg)
                rpm_neg_test = bam_test.count(region=rs_neg) / bam_test.mapped * 1E6
                rpm_neg_ctrl = bam_ctrl.count(region=rs_neg) / bam_ctrl.mapped * 1E6
                if rpm_neg_test <= rpm_neg_ctrl:
                    count_neg += 1
                if count_neg >= false_ct:
                    break
            else:
                break
        index_pos, count_pos, width_pos = 0, 0, 0
        while True:
            index_pos += skip
            ind_lt_pos, ind_rt_pos = cut + index_pos - w_rad, cut + index_pos + w_rad
            if ind_rt_pos <= hgsize[chr_i]:
                rs_pos = chr_i + ":" + str(ind_lt_pos) + "-" + str(ind_rt_pos)
                rpm_pos_test = bam_test.count(region=rs_pos) / bam_test.mapped * 1E6
                rpm_pos_ctrl = bam_ctrl.count(region=rs_pos) / bam_ctrl.mapped * 1E6
                if rpm_pos_test <= rpm_pos_ctrl:
                    count_pos += 1
                if count_pos >= false_ct:
                    break
            else:
                break
        span_rs = chr_i + ":" + str(cut + index_neg) + "-" + str(cut + index_pos)
        enrich_test = bam_test.count(region=span_rs) / bam_test.mapped * 1E6
        enrich_ctrl = bam_ctrl.count(region=span_rs) / bam_ctrl.mapped * 1E6
        bed_1, bed_2, bed_3 = chr_i, str(cut + index_neg), str(cut + index_pos)
        bed_4, bed_5, bed_6 = chr_i + ":" + str(cut), "%0.6f" % (enrich_test - enrich_ctrl), "+"
        bed_7, bed_8 = str(sta_i), str(end_i)
        outbed.write("\t".join([bed_1, bed_2, bed_3, bed_4, bed_5, bed_6, bed_7, bed_8]) + "\n")
        outnpy.append([rs, chr_i, str(cut), bed_2, bed_3, str(int(bed_3) - int(bed_2))])
    header = "region_string, chromosome, cut coordinate, " \
             "start coordinate, end coordinate, span width"
    np.savetxt(outpath + ".csv", np.asarray(outnpy), fmt='%s', delimiter=',', header=header)
    bam_test.close()
    bam_ctrl.close()
    outbed.close()


def gen_filter_dist(generator, distance):
    """ Given cut site generator, filter for cut sites that are separated by at least 'distance'.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param distance: minimum acceptable distance between adjacent cut sites on the same chromosome

    :yield another generator with cut sites too close to adjacent ones filtered out
    """
    chr_prev = None
    gen_prev, gen_curr = None, None
    for gen in generator:
        rs, cut = gen[0], gen[1]
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        if chr_prev != chr_i:
            if gen_prev and gen_prev[2]:
                yield gen_prev[0]
            chr_prev = chr_i
            gen_prev = (gen, cut, True)
        else:
            if cut - gen_prev[1] > distance:
                if gen_prev[2]:
                    yield gen_prev[0]
                gen_prev = (gen, cut, True)
            else:
                gen_prev = (gen, cut, False)


def rao_fourCseq_gen(generator, path_out, path_hic, kb_resolution, radius):
    """ Determine 4C-seq profiles using Hi-C data from Rao et al., 2014 at all viewpoints from a
    cut site generator. All data is written to a merged wiggle file.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param path_out: path to wiggle file (extension omitted) to write 4C-seq profile data
    :param path_hic: path to Hi-C "root" path from Rao et al., 2014 (e.g. "K562")
    :param kb_resolution: [integer] Hi-C resolution in kilobases {5, 10, 25, 50, 100, 250, 500}
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file
    """
    chr_prev = None
    chr_vals = None
    wigout = open(path_out + ".wig", 'w')
    for rs, cut, sen, pam, gui, mis, tar in generator:    # iterate over each target site
        chr_i = re.split('[:-]', rs)[0]
        if chr_prev != chr_i:
            print("rao_fourCseq_gen(): Processing %s." % chr_i)
            if chr_prev:            # save the first to second-to-last chromosome
                wigout.write("variableStep\tchrom=%s\n" % chr_prev)
                chr_vals = sorted(list(chr_vals), key=lambda x: x[0])
                for val in chr_vals:
                    wigout.write("%i\t%0.5f\n" % val)
            chr_vals = set(_rao_fourCseq_helper(path_hic, kb_resolution, chr_i, cut, radius))
            chr_prev = chr_i
        else:
            chr_vals |= set(_rao_fourCseq_helper(path_hic, kb_resolution, chr_i, cut, radius))
    # save last chromosome
    if chr_prev:
        wigout.write("variableStep\tchrom=%s\n" % chr_prev)
        chr_vals = sorted(list(chr_vals), key=lambda x: x[0])
        for val in chr_vals:
            wigout.write("%i\t%0.5f\n" % val)
    wigout.close()


def rao_fourCseq_single(path_out, path_hic, kb_resolution, chromosome, coordinate, radius=None):
    """ Determine 4C-seq profile at a single viewpoint using Hi-C data from Rao et al., 2014

    :param path_out: path to wiggle file (extension omitted) to write 4C-seq profile data
    :param path_hic: path to Hi-C "root" path from Rao et al., 2014 (e.g. "K562")
    :param kb_resolution: [integer] Hi-C resolution in kilobases {5, 10, 25, 50, 100, 250, 500}
    :param chromosome: [string] chromosome of viewpoint (e.g. "chr7")
    :param coordinate: [integer] coordinate of view point (e.g. 5529660)
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file
    """
    wigout = open(path_out + "rao_%i_%s_%i.wig" % (kb_resolution, chromosome, coordinate), 'w')
    wigout.write("variableStep\tchrom=%s\n" % chromosome)
    outvals = _rao_fourCseq_helper(path_hic, kb_resolution, chromosome, coordinate, radius)
    for val in outvals:
        wigout.write("%i\t%0.5f\n" % val)
    wigout.close()


def _rao_fourCseq_helper(path_hic, kb_resolution, chromosome, coordinate, radius=None):
    """ Determine 4C-seq profile at a specific viewpoint using Hi-C data from Rao et al., 2014
    Helper file with input being an open wiggle file in which to enter 4C-seq profile data

    :param path_hic: path to Hi-C "root" path from Rao et al., 2014 (e.g. "K562")
    :param kb_resolution: [integer] Hi-C resolution in kilobases {5, 10, 25, 50, 100, 250, 500}
    :param chromosome: [string] chromosome of viewpoint (e.g. "chr7")
    :param coordinate: [integer] coordinate of view point (e.g. 5529660)
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file

    :return outvals: [array] of (coordinate, value) tuples that correspond to the coordinate and
                     values for display in wiggle format.
    """
    path_file = os.path.join(path_hic, "%ikb_resolution_intrachromosomal" % kb_resolution, chromosome,
                             "MAPQGE30", "%s_%ikb.RAWobserved" % (chromosome, kb_resolution))
    round_coord = myround(coordinate, kb_resolution*1000)
    dM = np.loadtxt(path_file)  # load distance matrix
    outvals = []
    for i in range(dM.shape[0]):
        dm_i = dM[i, :]
        if dm_i[0] == round_coord or dm_i[1] == round_coord:
            dm_view, dm_coor, dm_valu = None, None, None
            if dm_i[0] == round_coord:
                dm_view, dm_coor, dm_valu = dm_i[0], dm_i[1], dm_i[2]
            elif dm_i[1] == round_coord:
                dm_view, dm_coor, dm_valu = dm_i[1], dm_i[0], dm_i[2]
            if not radius or not (dm_coor < coordinate - radius or dm_coor > coordinate + radius):
                outvals.append((dm_coor, dm_valu))
    return outvals


def h5_fourCseq_gen(generator, path_out, path_hic, radius):
    """ Determine 4C-seq profiles using Hi-C data from Rao et al., 2014 at all viewpoints from a
    cut site generator. All data is written to a merged wiggle file.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param path_out: path to wiggle file (extension omitted) to write 4C-seq profile data
    :param path_h5: path to Hi-C distance matrix from TODO
    :param kb_resolution: [integer] Hi-C resolution in kilobases {5, 10, 25, 50, 100, 250, 500}
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file
    """
    chr_prev = None
    chr_vals = None
    wigout = open(path_out + ".wig", 'w')
    for rs, cut, sen, pam, gui, mis, tar in generator:    # iterate over each target site
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        if chr_prev != chr_i:
            print("%s %i" % (chr_i, cut))
            if chr_prev:            # save the first to second-to-last chromosome
                wigout.write("variableStep\tchrom=%s\n" % chr_prev)
                chr_vals = sorted(list(chr_vals), key=lambda x: x[0])
                for val in chr_vals:
                    wigout.write("%i\t%0.5f\n" % val)
            chr_vals = set(_h5_fourCseq_helper(path_hic, chr_i, cut, radius))
            chr_prev = chr_i
        else:
            chr_vals |= set(_h5_fourCseq_helper(path_hic, chr_i, cut, radius))
    # save last chromosome
    if chr_prev:
        wigout.write("variableStep\tchrom=%s\n" % chr_prev)
        chr_vals = sorted(list(chr_vals), key=lambda x: x[0])
        for val in chr_vals:
            wigout.write("%i\t%0.5f\n" % val)
    wigout.close()


def h5_fourCseq_single(path_out, path_h5, chromosome, coordinate, radius=None):
    """ Determine 4C-seq profile at a single viewpoint using Hi-C data from Rao et al., 2014

    :param path_out: path to wiggle file (extension omitted) to write 4C-seq profile data
    :param path_h5: path to Hi-C distance matrix from TODO
    :param kb_resolution: [integer] Hi-C resolution in kilobases {5, 10, 25, 50, 100, 250, 500}
    :param chromosome: [string] chromosome of viewpoint (e.g. "chr7")
    :param coordinate: [integer] coordinate of view point (e.g. 5529660)
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file
    """
    wigout = open(path_out + "h5_%s_%i.wig" % (chromosome, coordinate), 'w')
    wigout.write("variableStep\tchrom=%s\n" % chromosome)
    outvals = _h5_fourCseq_helper(path_h5, chromosome, coordinate, radius)
    for val in outvals:
        wigout.write("%i\t%0.5f\n" % val)
    wigout.close()


def _h5_fourCseq_helper(path_h5, chromosome, coordinate, radius=None):
    """ Determine 4C-seq profile at a specific viewpoint using Hi-C data from Rao et al., 2014
    Helper file with input being an open wiggle file in which to enter 4C-seq profile data

    :param path_h5: path to Hi-C h5 file
    :param kb_resolution: [integer] Hi-C resolution in kilobases associated with h5 file
    :param chromosome: [string] chromosome of viewpoint (e.g. "chr7")
    :param coordinate: [integer] coordinate of view point (e.g. 5529660)
    :param radius: [integer] 4C-seq profile radius centered at coordinate to write on wiggle file

    :return outvals: [array] of (coordinate, value) tuples that correspond to the coordinate and
                     values for display in wiggle format.
    """
    h5 = h5py.File(path_h5, 'r')  # load distance matrix
    matrix = h5.get('matrix')
    intervals = h5.get('intervals')
    chrList = np.array(intervals.get('chr_list')).astype(str)
    staList = np.array(intervals.get('start_list'))
    endList = np.array(intervals.get('end_list'))
    data = matrix.get('data')
    indices = matrix.get('indices')
    indptr = matrix.get('indptr')
    dM = sparse.csr_matrix((data, indices, indptr))
    outvals = []
    cut_ind, sta_ind, end_ind = None, None, None
    thisChrInds = np.where(chrList.astype(str) == chromosome)[0]
    for i in thisChrInds:
        if staList[i] <= coordinate <= endList[i]:
            cut_ind = i
            break
    for i in thisChrInds:
        dm_coor, dm_valu = staList[i], dM[cut_ind, i]
        if not radius or not (dm_coor < coordinate - radius or dm_coor > coordinate + radius):
            outvals.append((dm_coor, dm_valu))
    return outvals


def myround(x, base):
    return base * round(x/base)

