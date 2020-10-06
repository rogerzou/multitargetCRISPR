# -*- coding: utf-8 -*-
""" Hi-C analysis of ChIP-seq after multi-targeting Cas9
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

import h5py
import os
import re
import numpy as np
from scipy import sparse


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
        [chr_i, sta_i, end_i] = re.split('[:-]', rs)
        if chr_prev != chr_i:
            print("%s %i" % (chr_i, cut))
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

