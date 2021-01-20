""" Analysis of gH2AX-high and gH2AX-low cut sites (based on ChIP-seq data)
"""
import numpy as np
import pandas as pd
import seaborn as sns
from tabulate import tabulate
from scipy import stats
import matplotlib.pyplot as plt


def load_files(workdir, path_gh2ax, path_atac, path_dnase, path_h3k4me3, path_h3k4me1, path_h3k9me3, path_h3k27ac,
               path_h3k36me3):
    """ Read csv files of gH2AX and histone marks containing RPM ChIP-seq enrichment in 5kb bins +/- 2Mb per cut site

    :param workdir: absolute path to directory containing all ChIP-seq csv files
    :param path_gh2ax: path to gH2AX csv, relative to workdir
    :param path_atac: path to ATAC-seq csv, relative to workdir
    :param path_dnase: path to DNase-seq csv, relative to workdir
    :param path_h3k4me3: path to H3K4me3 csv, relative to workdir
    :param path_h3k4me1: path to H3K4me1 csv, relative to workdir
    :param path_h3k9me3: path to H3K9me3 csv, relative to workdir
    :param path_h3k27ac: path to H3K27ac csv, relative to workdir
    :param path_h3k36me3: path to H3K36me3 csv, relative to workdir

    :return: pandas dataframes of csv files

    """
    # insulation scores, gh2ax
    gh2ax = pd.read_csv(workdir + path_gh2ax)
    # epigenetics
    atac = pd.read_csv(workdir + path_atac)
    dnase = pd.read_csv(workdir + path_dnase)
    h3k4me3 = pd.read_csv(workdir + path_h3k4me3)
    h3k4me1 = pd.read_csv(workdir + path_h3k4me1)
    h3k9me3 = pd.read_csv(workdir + path_h3k9me3)
    h3k27ac = pd.read_csv(workdir + path_h3k27ac)
    h3k36me3 = pd.read_csv(workdir + path_h3k36me3)

    return gh2ax, atac, dnase, h3k4me3, h3k4me1, h3k9me3, h3k27ac, h3k36me3


def process_files(gh2ax_df, atac_df, dnase_df, h3k4me3_df, h3k4me1_df, h3k9me3_df, h3k27ac_df, h3k36me3_df):
    """ Generator to process ChIP-seq dataframes to ensure correct structure and type

    :param gh2ax_df: dataframe of gH2AX RPM enrichment
    :param atac_df: dataframe of ATAC-seq RPM enrichment
    :param dnase_df: dataframe of DNase-seq RPM enrichment
    :param h3k4me3_df: dataframe of H3K4me3 RPM enrichment
    :param h3k4me1_df: dataframe of H3K4me1 RPM enrichment
    :param h3k9me3_df: dataframe of H3K9me3 RPM enrichment
    :param h3k27ac_df: dataframe of H3K27ac RPM enrichment
    :param h3k36me3_df: dataframe of H3K36me3 RPM enrichment

    :return generator giving NaN-free and float type dataframes
    """
    for df in [gh2ax_df, atac_df, dnase_df, h3k4me3_df, h3k4me1_df, h3k9me3_df, h3k27ac_df, h3k36me3_df]:
        df = df.replace(to_replace='None', value=0)  # replace none w/ zeros
        df.iloc[:, 2:] = df.iloc[:, 2:].astype(float)  # change all numbers to float type
        yield df


def gh2ax_enrichment(gh2ax_float_df):
    """ Extract top 10% and bottom 10% of cut sites using gH2AX RPM enrichment

    :param gh2ax_float_df: dataframe of gH2AX RPM enrichment (float type)
    :return indices (row #) of the original gH2AX dataframe for top 10% and bottom 10% of cut sites (by gH2AX signal)
    """
    rowlist = []                                        # list of rows in gH2AX dataframe
    for row in range(gh2ax_float_df.shape[0]):
        rowlist.append(list(gh2ax_float_df.iloc[row, 2:]))
    avgrow = [np.mean(i) for i in rowlist]              # Compute average gH2AX RPM signal per cut site
    pct10_num = round(gh2ax_float_df.shape[0] * 0.10)   # Number of sites with top 10% gH2AX signal
    ind_top10pct = np.argsort(avgrow)[-pct10_num:]      # Get original indices of top 10% gH2AX-enriched cut sites
    ind_bot10pct = np.argsort(avgrow)[:pct10_num]       # Get original indices of bottom 10% gH2AX-enriched cut sites

    return ind_top10pct, ind_bot10pct


def avg_chipseq(gh2ax_df, atac_df, dnase_df, h3k4me3_df, h3k4me1_df, h3k9me3_df, h3k27ac_df, h3k36me3_df):
    """ Compute average RPM enrichment per cut site for each histone mark/ATAC-seq/DNase-seq

    :param gh2ax_df: processed dataframe of gH2AX RPM enrichment
    :param atac_df: processed dataframe of ATAC-seq RPM enrichment
    :param dnase_df: processed dataframe of DNase-seq RPM enrichment
    :param h3k4me3_df: processed dataframe of H3K4me3 RPM enrichment
    :param h3k4me1_df: processed dataframe of H3K4me1 RPM enrichment
    :param h3k9me3_df: processed dataframe of H3K9me3 RPM enrichment
    :param h3k27ac_df: processed dataframe of H3K27ac RPM enrichment
    :param h3k36me3_df: processed dataframe of H3K36me3 RPM enrichment

    :return generator of numpy arrays of average RPM enrichment per cut site for all chromatin features

    """
    avg_atac_arr = []
    avg_h3k4me1_arr = []
    avg_h3k4me3_arr = []
    avg_h3k27ac_arr = []
    avg_h3k36me3_arr = []
    avg_h3k9me3_arr = []
    avg_dnase_arr = []

    dicts = [atac_df, h3k4me1_df, h3k4me3_df, h3k27ac_df, h3k36me3_df, h3k9me3_df, dnase_df]
    arr = [avg_atac_arr, avg_h3k4me1_arr, avg_h3k4me3_arr, avg_h3k27ac_arr, avg_h3k36me3_arr, avg_h3k9me3_arr, avg_dnase_arr]

    # Calculate average RPM absolute change +/- 2Mb each cut site for all chipseq data
    for i, j in zip(arr, dicts):
        for row in range(gh2ax_df.shape[0]):
            i.append(np.mean(j.iloc[row, 2:]))
    for k in arr:
        yield np.array(k)


def epigen(top_ind, bot_ind, avg_atac_arr, avg_h3k4me1_arr, avg_h3k4me3_arr, avg_h3k27ac_arr, avg_h3k36me3_arr,
           avg_h3k9me3_arr, avg_dnase_arr):
    """ Relate gH2AX-high and gH2AX-low regions to epigenetics

    :param top_ind: indices (row #) of the original gH2AX dataframe for top 10% of cut sites (by gH2AX signal)
    :param bot_ind: indices (row #) of the original gH2AX dataframe for bottom 10% of cut sites (by gH2AX signal)
    :param avg_atac_arr: numpy array of average ATAC-seq RPM enrichment per cut site
    :param avg_h3k4me1_arr: numpy array of average H3K4me1 RPM enrichment per cut site
    :param avg_h3k4me3_arr: numpy array of average H3k4me3 RPM enrichment per cut site
    :param avg_h3k27ac_arr: numpy array of average H3k27ac RPM enrichment per cut site
    :param avg_h3k36me3_arr: numpy array of average H3K36me3 RPM enrichment per cut site
    :param avg_h3k9me3_arr: numpy array of average H3k9me3 RPM enrichment per cut site
    :param avg_dnase_arr: numpy array of average DNase-seq RPM enrichment per cut site

    :return: - dataframe containing average RPM enrichment for shapes 1 and 2 for all epigenetic features
             - lists containing avg RPM enrichment of epigenetic features (by shape)

    - Note that 'shape 1' and 'shape 2' refer to two extremes on the gH2AX shape spectrum (i.e. gH2AX-high, gH2AX-low)
    """

    shape1_h3k4me3 = []
    shape2_h3k4me3 = []
    shape1_h3k27ac = []
    shape2_h3k27ac = []
    shape1_h3k36me3 = []
    shape2_h3k36me3 = []
    shape1_h3k4me1 = []
    shape2_h3k4me1 = []
    shape1_h3k9me3 = []
    shape2_h3k9me3 = []
    shape1_atac = []
    shape2_atac = []
    shape1_dnase = []
    shape2_dnase = []
    shape1_2_all_feat = []
    shape1_2_all = []

    # Add average RPM enrichment per cut site per 'shape' to arrays
    shape1 = [shape1_h3k4me3, shape1_h3k4me1, shape1_h3k36me3, shape1_h3k27ac, shape1_h3k9me3, shape1_atac,
              shape1_dnase]
    shape2 = [shape2_h3k4me3, shape2_h3k4me1, shape2_h3k36me3, shape2_h3k27ac, shape2_h3k9me3, shape2_atac,
              shape2_dnase]
    avg = [avg_h3k4me3_arr, avg_h3k4me1_arr, avg_h3k36me3_arr, avg_h3k27ac_arr,
           avg_h3k9me3_arr, avg_atac_arr, avg_dnase_arr]

    for i in top_ind:
        for j, k in zip(shape1, avg):
            j.append(k[i])
    for a in bot_ind:
        for b, c in zip(shape2, avg):
            b.append(c[a])
    for x in range(len(shape1)):
        shape1_2_all_feat.append(shape1[x])
        shape1_2_all_feat.append(shape2[x])
    for y, z in zip(shape1, shape2):
        shape1_2_all.extend(y + z)

    shape_labels = (['shape 1'] * len(shape1_h3k4me3) + ['shape 2'] * len(shape2_h3k4me3)) * 7
    shape1_2 = np.column_stack([shape1_2_all, shape_labels])
    shape1_2_df = pd.DataFrame(shape1_2, columns=['Per cut site average RPM enrichment', 'shape'])
    n = int(shape1_2_df.shape[0] / 7)
    marks_str = ['h3k4me3'] * n + ['h3k4me1'] * n + ['h3k36me3'] * n + ['h3k27ac'] * n + ['h3k9me3'] * n + ['atac'] * n + \
                ['dnase'] * n
    shape1_2_df['histone_marks/chromatin features'] = marks_str
    shape1_2_df.iloc[:, 0] = shape1_2_df.iloc[:, 0].astype(float)

    return shape1_2_df, shape1_2_all_feat


def mre11_groups(topind, botind, workdir, mre11_dir):
    """ Group MRE11 peaks by enrichment

    :param topind: indices (row #) of the original gH2AX dataframe for top 10% of cut sites (by gH2AX signal)
    :param botind: indices (row #) of the original gH2AX dataframe for bottom 10% of cut sites (by gH2AX signal)
    :param workdir: absolute path to directory containing all ChIP-seq csv files
    :param mre11_dir: path to MRE11 ChIP-seq csv, relative to workdir
    """
    mre11 = pd.read_csv(workdir + mre11_dir)

    mre11_shape1 = mre11.loc[list(topind)]
    mre11_shape2 = mre11.loc[list(botind)]
    mre11_shape1_rpm = list(mre11_shape1[' counts_RPM'])
    mre11_shape2_rpm = list(mre11_shape2[' counts_RPM'])
    mre11_shape2_rpm.remove(max(mre11_shape2_rpm))          # Remove first artifact from shape 2 (gH2AX low)
    mre11_shape2_rpm.remove(max(mre11_shape2_rpm))          # Remove second artifact from shape 2 (gH2AX low)
    mre11_shape1_2_rpm = mre11_shape1_rpm + mre11_shape2_rpm

    shape_labels_mre11 = ['shape 1'] * len(mre11_shape1_rpm) + ['shape 2'] * len(mre11_shape2_rpm)
    shape1_2_mre11 = np.column_stack([mre11_shape1_2_rpm, shape_labels_mre11])
    shape1_2_df_mre11 = pd.DataFrame(shape1_2_mre11, columns=['Per cut site MRE11 RPM enrichment', 'shape'])
    shape1_2_df_mre11.iloc[:, 0] = shape1_2_df_mre11.iloc[:, 0].astype(float)

    # Save boxplot
    plt.figure()
    sns.set_theme(style="whitegrid")
    ax_boxplt_mre11 = sns.boxplot(x='shape', y='Per cut site MRE11 RPM enrichment', order=['shape 1', 'shape 2'],
                                  data=shape1_2_df_mre11, palette="Set3")
    ax_boxplt_mre11.figure.savefig(workdir + 'boxplot_mre11.png')

    # t-test
    mre11_ttest = stats.ttest_ind(mre11_shape1_rpm, mre11_shape2_rpm)
    print(tabulate([['mre11', str(mre11_ttest[0]), str(mre11_ttest[1])]],
                   headers=['Chromatin feature', 't-statistic', 'p-value']))


def boxplot(shape1_2_df, workdir):
    """ Box plot of average epigenetic RPM enrichment for shapes/shapes 1 and 2, saved to working directory

    :param shape1_2_df: dataframe containing average RPM enrichment for shapes 1 and 2 for all epigenetic features
    :param workdir: absolute path to directory containing all ChIP-seq csv files
    """
    plt.figure()
    sns.set_theme(style="whitegrid")
    ax_boxplt = sns.boxplot(x='histone_marks/chromatin features', y='Per cut site average RPM enrichment', hue='shape',
                            data=shape1_2_df, palette="Set3")
    ax_boxplt.figure.savefig(workdir + 'boxplot_epigen.png')


def ttest(shape1_2_all_feat):
    """ Two-sample t test of epigenetic RPM enrichment between shapes/shapes 1 and 2; results printed to screen

    :param shape1_2_all_feat: dataframe containing average RPM enrichment for shapes 1 and 2 for all epigenetic features
    """
    ttest_results = []
    for i in range(0, len(shape1_2_all_feat) - 1, 2):
        ttest_results.extend(stats.ttest_ind(shape1_2_all_feat[i], shape1_2_all_feat[i + 1]))

    print(tabulate([['h3k4me3', str(ttest_results[0]), str(ttest_results[1])],
                    ['h3k4me1', str(ttest_results[2]), str(ttest_results[3])],
                    ['h3k36me3', str(ttest_results[4]), str(ttest_results[5])],
                    ['h3k27ac', str(ttest_results[6]), str(ttest_results[7])],
                    ['h3k9me3', str(ttest_results[8]), str(ttest_results[9])],
                    ['atac', str(ttest_results[10]), str(ttest_results[11])],
                    ['dnase', str(ttest_results[12]), str(ttest_results[13])]],
                   headers=['Chromatin feature', 't-statistic', 'p-value']))


def to_bed(gh2ax_df, top_ind, bot_ind, workfolder, span_rad=50E3):
    """ Write shape 1 and 2 genomic coordinates to BED files, saved in working directory

    :param gh2ax_df: processed dataframe of gH2AX RPM enrichment
    :param top_ind: indices (row #) of shape 1 (top 10% by enrichment) in original gH2AX dataframe
    :param bot_ind: indices (row #) of shape 2 (bottom 10% by enrichment) in original gH2AX dataframe
    :param workfolder: absolute path to directory containing all ChIP-seq csv files
    :param span_rad: radius (bp) of region centered at cut site
    """
    with open(workfolder + 'dtw_kmeans_gh2ax/top10pct.bed', 'w') as f1:
        for row in top_ind:
            chr_i_merged = gh2ax_df.iloc[row, 0]
            reg_sta = gh2ax_df.iloc[row, 1] - span_rad
            reg_end = gh2ax_df.iloc[row, 1] + span_rad
            if reg_sta < 0:
                raise ValueError('Genomic coordinates negative! Please adjust radius.')
            f1.write(chr_i_merged + '\t' + str(int(reg_sta)) + '\t' + str(int(reg_end)) + '\n')

    with open(workfolder + 'dtw_kmeans_gh2ax/bot10pct.bed', 'w') as f2:
        for row_i in bot_ind:
            chr_i_merged = gh2ax_df.iloc[row_i, 0]
            reg_sta = gh2ax_df.iloc[row_i, 1] - span_rad
            reg_end = gh2ax_df.iloc[row_i, 1] + span_rad
            if reg_sta < 0:
                raise ValueError('Genomic coordinates negative! Please adjust radius.')
            f2.write(chr_i_merged + '\t' + str(int(reg_sta)) + '\t' + str(int(reg_end)) + '\n')

