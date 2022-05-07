# -*- coding: utf-8 -*-
""" Anallysis of Chakrabarti et al., Mol Cell, 2019
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

import os
import subprocess as sp
import numpy as np
from . import mtss as m
from . import chipseq as c
from lib.pyliftover import liftover


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


def chakrabarti_generator(span_r, genome):
    """ Generator to yield all target sites from Chakrabarti et al., Mol Cell, 2019.

    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param genome: path to genome (hg38 - with .fa extension)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM , discovered protospacer,
            # mismatches, non-mismatched protospacer )

    """
    d = load_chakrabarti()
    hgsize = c.get_genome_dict(genome[0])
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
        cent_end = min(hgsize[chr_i], end_tmp + 250)
        cent_rs = "%s:%i-%i" % (chr_i, cent_sta, cent_end)
        cent_faidx = sp.check_output(['samtools', 'faidx', genome[1], cent_rs]).split()
        seq = (b"".join(cent_faidx[1:]).upper()).decode("utf-8")
        cand = m.sub_findmis(seq, guide, 0)
        if cand is not None and len(cand) > 0:
            cut_i = cent_sta + cand[0][2]
            sen_i = '+' if cand[0][5] == 1 else '-'
            gui_i = cand[0][3]
            mis_i = cand[0][1]
            pam_i = cand[0][4]
            if pam_i != pam_init:
                print("chakrabarti_generator(): PAM DIFFERENT!!!!!")
            span_sta = max(1, cut_i - span_r)
            span_end = min(hgsize[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide
        else:
            print("chakrabarti_generator(): Did not find a particular sequence.")


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
