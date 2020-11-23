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
from scipy import sparse, stats
from . import chipseq as c
from . import mtss as m

CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX']


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
        chr_i = re.split('[:-]', rs)[0]
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
        if chr_i in CHR:
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
    path_file = os.path.join(path_hic, "%ikb_resolution_intrachromosomal" % kb_resolution,
                             chromosome,
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


def absolute_change_from_cutsite(generator, genome, f_test, f_ctrl, outpath,
                                 span_rad=50000, res=5000, wind_rad=2000000):
    """ Determine absolute change in enrichment between test and ctrl at positions moving away from
        the cut site, either in the positive (downstream) or negative (upstream) direction.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg19', path/to/hg19.fa]
    :param f_test: test sample BAM file
    :param f_ctrl: negative control BAM file
    :param outpath: path to output BED file (.bed extension omitted)
    :param span_rad: radius of window of enrichment evaluation at each site
    :param res: number of bases to skip per evaluation of enrichment over control
    :param wind_rad: radius in bp centered at the cut site to consider

    output: csv file with first two columns corresponding to the chr and coordinate of the cut site,
            followed by each column listing enrichment values taken from coordinates closest to the
            cut site to furthest from the cut site, i.e. increasing genomic coordinates for
            downstream, or decreasing genomic coordinates for upstream.
    """
    hgsize = c.hg_dict(genome[0])
    bam_test, bam_ctrl = pysam.AlignmentFile(f_test, 'rb'), pysam.AlignmentFile(f_ctrl, 'rb')
    n_iter = int(wind_rad / res) + 1
    all_neg, all_pos, indices_neg, indices_pos = [], [], None, None
    for rs, cut, sen, pam, gui, mis, guide in generator:
        chr_i = re.split('[:-]', rs)[0]
        if chr_i in CHR:
            # Get absolute change in the downstream direction
            index_neg, tracker_neg, indices_neg = 0, [chr_i, cut], ['None', 'None']
            for i in range(n_iter):
                indices_neg.append(str(index_neg / 1E6))
                ind_lt_neg, ind_rt_neg = cut + index_neg - span_rad, cut + index_neg + span_rad
                if ind_lt_neg >= 0:
                    rs_neg = chr_i + ":" + str(ind_lt_neg) + "-" + str(ind_rt_neg)
                    rpm_neg_test = bam_test.count(region=rs_neg) / bam_test.mapped * 1E6
                    rpm_neg_ctrl = bam_ctrl.count(region=rs_neg) / bam_ctrl.mapped * 1E6
                    tracker_neg.append(rpm_neg_test - rpm_neg_ctrl)
                else:
                    tracker_neg.append(None)
                index_neg -= res
            all_neg.append(tracker_neg)
            # Get absolute change in the upstream direction
            index_pos, tracker_pos, indices_pos = 0, [chr_i, cut], ['None', 'None']
            for i in range(n_iter):
                indices_pos.append(str(index_pos / 1E6))
                ind_lt_pos, ind_rt_pos = cut + index_pos - span_rad, cut + index_pos + span_rad
                if ind_rt_pos < hgsize[chr_i]:
                    rs_pos = chr_i + ":" + str(ind_lt_pos) + "-" + str(ind_rt_pos)
                    rpm_pos_test = bam_test.count(region=rs_pos) / bam_test.mapped * 1E6
                    rpm_pos_ctrl = bam_ctrl.count(region=rs_pos) / bam_ctrl.mapped * 1E6
                    tracker_pos.append(rpm_pos_test - rpm_pos_ctrl)
                else:
                    tracker_pos.append(None)
                index_pos += res
            all_pos.append(tracker_pos)
    np.savetxt(outpath + "_achange_neg.csv", np.asarray(all_neg), fmt='%s', delimiter=',',
               header=", ".join(indices_neg))
    np.savetxt(outpath + "_achange_pos.csv", np.asarray(all_pos), fmt='%s', delimiter=',',
               header=", ".join(indices_pos))
    bam_test.close()
    bam_ctrl.close()


def wigvals_from_cutsite(generator, genome, f_wig, outpath, res=5000, wind_rad=2000000):
    """ Determine insulation scores at positions moving away from the cut site, either in the
        positive (downstream) or negative (upstream) direction.

    :param generator: generator that outputs target sites in the following tuple format:
                ( span_rs   =   region string in "chr1:100-200" format, centered at cut site
                  cut_i     =   cut site                 (int)
                  sen_i     =   sense/antisense          (+/- str)
                  pam_i     =   PAM                      (str)
                  gui_i     =   genomic target sequence  (str)
                  mis_i     =   # mismatches             (int)
                  guide     =   intended target sequence (str)
    :param genome: [genome name, path to genome with .fa extension], i.e. ['hg19', path/to/hg19.fa]
    :param f_wig: test sample BAM file
    :param outpath: path to output BED file (.bed extension omitted)
    :param res: number of bases to skip per evaluation
    :param wind_rad: radius in bp centered at the cut site to consider

    output: csv file with first two columns corresponding to the chr and coordinate of the cut site,
            followed by each column listing insulation scores taken from coordinates closest to the
            cut site to furthest from the cut site, i.e. increasing genomic coordinates for
            downstream, or decreasing genomic coordinates for upstream.
    """
    hgsize = c.hg_dict(genome[0])
    wig = Wig(f_wig)
    n_iter = int(wind_rad / res) + 1
    all_neg, all_pos, indices_neg, indices_pos = [], [], None, None
    for rs, cut, sen, pam, gui, mis, guide in generator:
        chr_i = re.split('[:-]', rs)[0]
        if chr_i in CHR:
            # Get insulation score in the downstream direction
            index_neg, tracker_neg, indices_neg = 0, [chr_i, cut], ['None', 'None']
            for i in range(n_iter):
                indices_neg.append(str(index_neg / 1E6))
                coor_neg = cut + index_neg
                if coor_neg >= 0:
                    tracker_neg.append(wig.get_value(chr_i, coor_neg))
                else:
                    tracker_neg.append(None)
                index_neg -= res
            all_neg.append(tracker_neg)
            # Get insulation score in the upstream direction
            index_pos, tracker_pos, indices_pos = 0, [chr_i, cut], ['None', 'None']
            for i in range(n_iter):
                indices_pos.append(str(index_pos / 1E6))
                coor_pos = cut + index_pos
                if coor_pos < hgsize[chr_i]:
                    tracker_pos.append(wig.get_value(chr_i, coor_pos))
                else:
                    tracker_pos.append(None)
                index_pos += res
            all_pos.append(tracker_pos)
    np.savetxt(outpath + "_wigvals_neg.csv", np.asarray(all_neg), fmt='%s', delimiter=',',
               header=", ".join(indices_neg))
    np.savetxt(outpath + "_wigvals_pos.csv", np.asarray(all_pos), fmt='%s', delimiter=',',
               header=", ".join(indices_pos))


def merged_from_cutsite(f_data_pos, f_data_neg, outpath):
    """ Merge both downstream (pos) and upstream (neg) outputs from absolute_change_from_cutsite()
        or wigvals_from_cutsite()

    :param f_data_pos: downstream data from absolute_change_from_cutsite()/wigvals_from_cutsite()
    :param f_data_neg: upstream data from absolute_change_from_cutsite()/wigvals_from_cutsite()
    :param outpath: path to output file, extension omitted; function will add "_merged.csv"

    output: csv file with first two columns corresponding to the chr and coordinate of the cut site,
            followed by each column listing values taken from coordinates most negative (upstream)
            to positive (downstream).
    """
    ind_pos = m.load_npheader(f_data_pos).split(', ')
    ind_neg = m.load_npheader(f_data_neg).split(', ')
    data_pos, data_neg = m.load_nparray(f_data_pos), m.load_nparray(f_data_neg)
    all_track_merged = np.hstack((data_neg[:, :2], np.fliplr(data_neg[:, 3:]), data_pos[:, 2:]))
    hmerged = ', '.join(np.hstack((ind_neg[:2], np.flip(ind_neg[3:]), ind_pos[2:])))
    np.savetxt(outpath + "_merged.csv", all_track_merged, fmt='%s', delimiter=',', header=hmerged)


def derivative_from_cutsite(f_data, outpath, rad_npts=5):
    """ Computes the derivative of outputs from absolute_change_from_cutsite() or
        wigvals_from_cutsite()

    :param f_data: input, i.e. output from absolute_change_from_cutsite() or wigvals_from_cutsite()
    :param outpath: path to output file, extension omitted; function will add "_delta.csv"
    :param rad_npts: number of points from either side (i.e. radius) to consider when
                     calculating derivative using stats.linregress()

    output: csv file with first two columns corresponding to the chr and coordinate of the cut site,
            followed by each column listing derivative values taken from coordinates closest to the
            cut site to furthest from the cut site, i.e. increasing genomic coordinates for
            downstream, or decreasing genomic coordinates for upstream.
    """
    header = m.load_npheader(f_data)
    data = m.load_nparray(f_data)
    outdata = []
    for data_i in data:
        if 'None' not in data_i[0]:
            lendata_i = len(data_i) - 2
            out_d_i = np.zeros_like(data_i, dtype=object)
            out_d_i[:2] = data_i[:2]
            data_i = data_i[2:]
            for j in range(lendata_i):
                left_ind = max(0, j - rad_npts)
                right_ind = min(lendata_i, j + rad_npts + 1)
                regress_d = data_i[left_ind:right_ind]
                if 'None' in regress_d:
                    out_d_i[j + 2] = 'None'
                else:
                    out_d_i[j + 2] = stats.linregress(np.arange(right_ind - left_ind),
                                                      regress_d.astype(float))[0]
            outdata.append(out_d_i)
        else:
            outdata.append(data_i)
    np.savetxt(outpath + "_delta.csv", np.asarray(outdata), fmt='%s', delimiter=',', header=header)


def getXy_insulation(f_ins_pos, f_ins_neg, f_readcts, f_span, outpath, ctonly=False):
    """
    :param f_ins_pos: downstream data from wigvals_from_cutsite()
    :param f_ins_neg: upstream data from wigvals_from_cutsite()
    :param f_readcts: read_counts() data (i.e. from mre11)
    """
    data_pos = m.load_nparray(f_ins_pos)[:, 2:]
    data_neg = m.load_nparray(f_ins_neg)[:, 2:]
    feat_1 = m.load_nparray(f_readcts)[:, 5]
    y = m.load_nparray(f_span)[:, 5]
    datanone = np.logical_or(data_pos == 'None', data_neg == 'None')
    y = y[~datanone.any(axis=1)].astype(float) / 1E6
    feat_1 = feat_1[~datanone.any(axis=1)].astype(float)
    data_pos = data_pos[~datanone.any(axis=1)].astype(float)
    data_neg = data_neg[~datanone.any(axis=1)].astype(float)
    feat_2 = np.mean(np.concatenate((data_pos[:, :5], data_neg[:, :5]), axis=1), axis=1)
    data = np.concatenate((data_pos, data_neg), axis=1)
    num_l0 = np.sum(np.ones_like(data), axis=1, where=(data <= 0))
    num_g0 = np.sum(np.ones_like(data), axis=1, where=(data > 0))
    feat_3 = num_l0 / (num_l0 + num_g0)
    feat_4 = np.min(data, axis=1)
    feat_5 = np.max(data, axis=1)
    if ctonly:
        X = feat_1.reshape(-1, 1)
    else:
        X = np.transpose(np.vstack((feat_1, feat_2, feat_3, feat_4, feat_5)))
    np.savetxt(outpath + "_features.csv",
               np.transpose(np.vstack((y, feat_1, feat_2, feat_3, feat_4, feat_5))),
               fmt='%s', delimiter=',', header="y, rc, I_at_cut, %_I_neg, I_min, I_max")
    return X, y


def dtw(s, t, window=10):
    """ Implementation of Dynamic Time Warping from
        https://towardsdatascience.com/dynamic-time-warping-3933f25fcdd
    """
    n, m = len(s), len(t)
    w = np.max([window, abs(n - m)])
    dtw_matrix = np.zeros((n + 1, m + 1))

    for i in range(n + 1):
        for j in range(m + 1):
            dtw_matrix[i, j] = np.inf
    dtw_matrix[0, 0] = 0

    for i in range(1, n + 1):
        for j in range(np.max([1, i - w]), np.min([m, i + w]) + 1):
            dtw_matrix[i, j] = 0

    for i in range(1, n + 1):
        for j in range(np.max([1, i - w]), np.min([m, i + w]) + 1):
            cost = abs(s[i - 1] - t[j - 1])
            # take last min from a square box
            last_min = np.min(
                [dtw_matrix[i - 1, j], dtw_matrix[i, j - 1], dtw_matrix[i - 1, j - 1]])
            dtw_matrix[i, j] = cost + last_min
    return dtw_matrix


def dtw_randomize(f_data, f_score, outpath, numrand=100):
    """ Computes Dynamic Time Warping value between
        (1) ChIP-seq enrichment - output of absolute_change_from_cutsite()
        (2) insulation scores - output of wigvals_from_cutsite()
        Permutation of insulation scores to match with ChIP-seq enrichment is performed after the
        first iteration (# of iterations specified by numrand variable)

    :param f_data: path to ChIP-seq enrichent output from absolute_change_from_cutsite()
    :param f_score: path to insulation score output from wigvals_from_cutsite()
    :param outpath: path to output file, extension omitted; function will add "_dtw-rand.csv"
    :param numrand: (int >= 1) number of iterations to perform, where the first one has no
                    permutation of insulation score, while the remaining ones do.

    output: csv file with first two columns recording chromosome and coordinate of cut site, then
            subsequent columns contain values from dynamic time warping. Iterations after the first
            one have the insulation scores permuted. Rows correspond to cut sites.
    """
    data, score = m.load_nparray(f_data), m.load_nparray(f_score)
    header = "chr, coord, " + ', '.join(['iter_%i' % n for n in range(numrand)])
    dtwscore = []
    datanone = data == 'None'
    data, score = data[~datanone.any(axis=1)], score[~datanone.any(axis=1)]  # remove rows with None
    for n in range(numrand):
        if n > 0:
            np.random.shuffle(score)
        if n % 10 == 0:
            print("dtw_randomize(): Processing %i / %i." % (n, numrand))
        for i, (data_i, score_i) in enumerate(zip(data, score)):
            dtw_matrix = dtw(data_i[2:].astype(float), score_i[2:].astype(float))
            if n == 0:
                dtwscore.append([data_i[0], data_i[1], dtw_matrix.flat[-1]])
            else:
                dtwscore[i].append(dtw_matrix.flat[-1])
    np.savetxt(outpath + "_dtw-rand.csv", np.asarray(dtwscore), fmt='%s', delimiter=',',
               header=header)


def categorize_by_insulation_randomize(f_data, f_score, outpath, numrand=100):
    """ Determines average values of input data (output of absolute_change_gen()) that correspond to
        insulation scores (output of wigvals_from_cutsite()) either <=0 or >0.
        Permutation of insulation scores to match with ChIP-seq enrichment is performed after the
        first iteration (# of iterations specified by numrand variable)

    :param f_data: path to ChIP-seq enrichent output from absolute_change_from_cutsite()
    :param f_score: path to insulation score output from wigvals_from_cutsite()
    :param outpath: path to output file, extension omitted; function will add "_cat-rand.csv"
    :param numrand: (int >= 1) number of iterations to perform, where the first one has no
                    permutation of insulation score, while the remaining ones do.

    output: csv file with first two columns recording chromosome and coordinate of cut site, then
            two column each corresponds to average values of input data where insulation scores
            are <= 0 or >0. Iterations after the first one have the insulation scores permuted.
            Rows correspond to cut sites.
    """
    data, score = m.load_nparray(f_data), m.load_nparray(f_score)
    header = "chr, coord, " + ', '.join(['>0_%i, <=_%i' % (n, n) for n in range(numrand)])
    catscore = []
    datanone = data == 'None'
    data, score = data[~datanone.any(axis=1)], score[~datanone.any(axis=1)]  # remove rows with None
    for n in range(numrand):
        if n > 0:
            np.random.shuffle(score)
        if n % 10 == 0:
            print("categorize_by_insulation_randomize(): Processing %i / %i." % (n, numrand))
        for i, (data_i, score_i) in enumerate(zip(data, score)):
            greater = data_i[2:][score_i[2:].astype(float) > 0].astype(float)
            lesser = data_i[2:][score_i[2:].astype(float) <= 0].astype(float)
            if n == 0:
                catscore.append([data_i[0], data_i[1], np.mean(greater), np.mean(lesser)])
            else:
                catscore[i].append(np.mean(greater))
                catscore[i].append(np.mean(lesser))
    np.savetxt(outpath + "_cat-rand.csv", np.asarray(catscore), fmt='%s', delimiter=',',
               header=header)


class Wig:
    """ Class that allows easy querying of wiggle files. Supports both fixedStep and variableStep.
        Assumes canonical formatting of form:
            variableStep chrom=chrN
        or
            fixedStep chrom=chrN start=position step=stepInterval
        where each word is separated by a space.
    """
    def __init__(self, wigfile):
        """ Inputs a variableStep or fixedStep wigfile. Writes to two dictionaries that contain
            the values and indices of the wiggle for easy querying.
        """
        self.file = wigfile
        self.D, self.I = {}, {}
        fixed_step = True
        with open(wigfile) as f:
            start_i, step_i, iter_i = None, None, None
            for line in f:
                line = line.strip().split(' ')
                if line[0] == "track":
                    continue
                elif line[0] == "fixedStep":
                    chr_i = line[1].split('=')[1]
                    start_i = int(line[2].split('=')[1])
                    step_i = int(line[3].split('=')[1])
                    iter_i = 0
                    if chr_i not in self.D:
                        self.D[chr_i], self.I[chr_i] = [], []
                    else:
                        raise ValueError("Wig(): wig file has duplicate chromosome annotations!")
                elif line[0] == "variableStep":
                    fixed_step = False
                    chr_i = line[1].split('=')[1]
                    if chr_i not in self.D:
                        self.D[chr_i], self.I[chr_i] = [], []
                    else:
                        raise ValueError("Wig(): wig file has duplicate chromosome annotations!")
                else:
                    if fixed_step:
                        self.I[chr_i].append(start_i + iter_i * step_i)
                        self.D[chr_i].append(float(line[0]))
                        iter_i += 1
                    else:
                        self.I[chr_i].append(int(line[0]))
                        self.D[chr_i].append(float(line[1]))

    def get_value(self, chromosome, coordinate):
        """ Function call that will output the value at any genomic coordinate. For coordinates
            that are between two coordinates actually specified by the wiggle file, linear
            interpolation is employed.
        """
        if chromosome in self.D:
            Dchr, Ichr = self.D[chromosome], self.I[chromosome]
            for i, ind_i in enumerate(Ichr):
                if ind_i >= coordinate:
                    i1, i2 = i-1, i
                    if i1 < 0:
                        return Dchr[i2]
                    else:
                        val1, val2 = Dchr[i1], Dchr[i2]
                        ind1, ind2 = Ichr[i1], Ichr[i2]
                        return val1 + (val2 - val1) / (ind2 - ind1) * (coordinate - ind1)
            return Dchr[-1]
        else:
            print("Wig() get_value(): queried chromosome %s is not found." % chromosome)
            return None
