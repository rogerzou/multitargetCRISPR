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
import re
import numpy as np
import sys
sys.path.append("..") # Adds higher directory to python modules path.
import csv
from scipy import stats
from . import chipseq as c
import string
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, train_test_split
from sklearn.svm import SVR
from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.inspection import permutation_importance
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import normalize, LabelEncoder, OneHotEncoder
from lib.pyliftover import liftover

from matplotlib import pyplot as pp

OLD_CHAR = "ACGT"
NEW_CHAR = "TGCA"
# WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 0.8, 0.8, 0.8, 0.8, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
WEIGHT = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


def get_reverse_complement(seq):
    return seq.translate(str.maketrans(OLD_CHAR, NEW_CHAR))[::-1]


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


def macs_generator(peak, span_r, genome, guide, mismatch=10, cent_r=100):
    """ Generator

    :param peak: macs2 output file
    :param span_r:
    :param radius: radius of window from peak center to search for cut site
    :param genome:
    :param guide:
    :param mismatch:

    """
    m_out = np.loadtxt(peak, dtype=object)    # load macs2 narrowPeak output file
    numpeaks = m_out.shape[0]                 # number of peaks from blender
    hg38size = c.hg38_dict()
    for i in range(numpeaks):
        chr_i = m_out[i, 0]
        if chr_i in hg38size:
            # center = int(np.mean([int(m_out[i, 1]), int(m_out[i, 2])]))
            center = int(m_out[i, 9]) + int(m_out[i, 1])
            span_sta = max(1, center - span_r)
            span_end = min(hg38size[chr_i], center + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            cent_sta = max(1, center - cent_r)
            cent_end = min(hg38size[chr_i], center + cent_r)
            cent_rs = "%s:%i-%i" % (chr_i, cent_sta, cent_end)
            cent_faidx = sp.check_output(['samtools', 'faidx', genome, cent_rs]).split()
            seq = (b"".join(cent_faidx[1:]).upper()).decode("utf-8")
            cand = sub_findmis(seq, guide, mismatch)
            if cand is not None and len(cand) > 0:
                cut_i = cent_sta + cand[0][2]
                sen_i = '+' if cand[0][4] == 1 else '-'
                pam_i = "NGG"
                gui_i = cand[0][3]
                mis_i = cand[0][1]
                yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide


def chakrabarti_generator(span_r, genome):
    dirname, filename = os.path.split(os.path.abspath(__file__))
    d = np.loadtxt(os.path.dirname(dirname) + "/lib/chakrabarti_mm.csv",
                   delimiter=',', dtype=object, skiprows=1)
    numtargets = d.shape[0]
    hg38size = c.hg38_dict()
    lo = liftover.LiftOver(open(os.path.dirname(dirname) + "/lib/hg19ToHg38.over.chain"))
    for i in range(numtargets):
        guide = d[i, 1][:-3]
        pam_i = d[i, 1][-3:]
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
            sen_i = '+' if cand[0][4] == 1 else '-'
            gui_i = cand[0][3]
            mis_i = cand[0][1]
            span_sta = max(1, cut_i - span_r)
            span_end = min(hg38size[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            yield span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide
        else:
            print("oh no")


def sub_findmis(s, matchstr, maxmismatch):
    sublen = len(matchstr)
    candidates = []

    # check mismatches for original string
    for i in range(len(s)-sublen-6):
        i += 3
        # sense
        substr_1 = s[i:i + sublen]
        pam_1 = s[i + sublen:i + sublen + 3]
        cutsite_1 = i + sublen - 4
        wlist = [x + 1 if substr_1[j] != matchstr[j] else x for j, x in enumerate([0] * sublen)]
        mismatches = sum(wlist)
        if mismatches <= maxmismatch and pam_1 in ['AGG', 'CGG', 'TGG', 'GGG']:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_1, substr_1, 1))
        # anti-sense
        substr_0 = get_reverse_complement(s[i:i + sublen])
        pam_0 = get_reverse_complement(s[i - 3:i])
        cutsite_0 = i + 4
        wlist = [x + 1 if substr_0[j] != matchstr[j] else x for j, x in enumerate([0] * sublen)]
        mismatches = sum(wlist)
        if mismatches <= maxmismatch and pam_0 in ['AGG', 'CGG', 'TGG', 'GGG']:
            score = sum([a*b for a, b in zip(wlist, WEIGHT)])
            candidates.append((score, mismatches, cutsite_0, substr_0, 0))

    candidates.sort(key=lambda y: y[0])
    return candidates


def alu_read_subsets(generator, filein, fileout=None):
    bam = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    LS = []
    for rs, cut, sen, pam, gui, mis, guide in generator:
        ctM, ctN, ctL, ctR, ctT = 0, 0, 0, 0, 0
        for read1, read2 in c.read_pair_generator(bam, rs):
            read = c.read_pair_align(read1, read2)
            if not read:
                continue
            ctT += 1
            if read[0] < cut < read[3]:  # fragments that span
                ctM += 1
            elif cut + 5 >= read[0] >= cut:  # fragments that begin 5bp of cleavage site
                ctR += 1
                ctN += 1
            elif cut - 5 <= read[-1] <= cut:  # fragments that end 5bp of cleavage site
                ctL += 1
                ctN += 1
        fpm = bam.mapped / 2E6      # fragments per millon
        ctM /= fpm                  # fragments that span
        ctN /= fpm                  # fragments that don't span
        ctL /= fpm                  # fragments that don't span, on left
        ctR /= fpm                  # fragments that don't span, on right
        ctT /= fpm                  # count total number of reads
        if pam == 'NGG':
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
                LS.append((rs, cut, sen, gui, guide, mis, ig[0], ig[1], ctT, ctM, ctN, cUp, cDo, cPr, cDi))
            else:
                LS.append((rs, cut, sen, gui, guide, mis, "", "", ctT, ctM, ctN, "", "", cPr, cDi))
    bam.close()
    header = "region_string, cut_site, sense, matched guide, original guide, mismatch, gene, " \
             "orientation, Ctotal, Cspan, Cend, Cupstream, Cdownstream, Cproximal, Cdistal"
    LS = np.asarray(LS)
    if fileout:
        np.savetxt(fileout, LS, fmt='%s', delimiter=',', header=header)
    return LS


def alu_read_counts(generator, filein, fileout=None):
    bam = pysam.AlignmentFile(filein, 'rb')             # BAM file for analysis
    list_stat = []
    for rs, cut, sen, pam, gui, mis in generator:
        if pam == 'NGG':
            ct_rpm = bam.count(region=rs) / bam.mapped * 1E6
            list_stat.append((rs, cut, sen, gui, mis, ct_rpm))
    bam.close()
    list_stat = np.asarray(list_stat)
    if fileout:
        np.savetxt(fileout, list_stat, fmt='%s', delimiter=',')
    return list_stat


def load_nparray(array):
    return np.loadtxt(array, dtype=object, delimiter=',')


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


def pca(X, num_components=6):
    pca = PCA(n_components=num_components)
    X = pca.fit_transform(X)
    return X


def NeuralNetworkTrainDefault(X, y, modelfile, solver='lbfgs', alpha=0.1, hidden_layer_sizes=(8,)):
    # Define and fit base regressor
    mlp = MLPRegressor(solver=solver, max_iter=5000, random_state=42,
                       hidden_layer_sizes=hidden_layer_sizes, alpha=alpha)
    mlp.fit(X, y)
    # calculate results (correlation coefficient)
    CORREL(mlp, X, y)
    # save and return base estimator
    pickle.dump(mlp, open(modelfile, 'wb'))
    return mlp


def NeuralNetworkTrainGridCV(X, y, modelfile):
    # Define base regressor:
    mlp = MLPRegressor(max_iter=5000, random_state=42)
    # Define search space:
    params = {

        'solver': ['lbfgs'],
        'alpha': [0.001, 0.01, 0.1, 1, 2],
        'hidden_layer_sizes': [(5,), (6,)]
    }
    # Find best hyperparameters and then refit on all training data
    mlp_grid = GridSearchCV(estimator=mlp, param_grid=params, cv=5, verbose=2, n_jobs=8)
    mlp_grid.fit(X, y)
    # get the best estimator and calculate results (correlation coefficient)
    best_grid = mlp_grid.best_estimator_
    CORREL(best_grid, X, y)
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


def ModelTest(X, y, modelfile):
    model = pickle.load(open(modelfile, 'rb'))
    print(modelfile)
    print(model)
    CORREL(model, X, y)
    np.savetxt(modelfile[0:-4] + "_out.csv", np.column_stack((y, model.predict(X))), fmt='%s', delimiter=',')


def FeatureImportance(X, y, modelfile, labels):
    model = pickle.load(open(modelfile, 'rb'))
    labels = np.asarray(labels, dtype=object)

    result = permutation_importance(model, X, y, n_repeats=10, random_state=42, n_jobs=8)
    sorted_idx = result.importances_mean.argsort()

    fig, ax = pp.subplots()
    ax.boxplot(result.importances[sorted_idx].T,
               vert=False, labels=labels[sorted_idx])
    ax.set_title("Permutation Importances (test set)")
    fig.tight_layout()
    pp.show()


def data_split(X, y, test_size=0.3):
    return train_test_split(X, y, test_size=test_size, shuffle=True, random_state=42)


def MAPE(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy


def CORREL(model, test_features, test_labels):
    corr = np.corrcoef(test_labels, model.predict(test_features))
    print('Model Performance')
    print('Correlation coefficient = {:0.2f}%.'.format(corr[1, 0]))
    return


""" deprecated functions """

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

