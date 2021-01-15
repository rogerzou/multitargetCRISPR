""" Dynamic time warping and k-means clustering of gH2AX ChIP-seq data
"""

import numpy as np
import pandas as pd
import seaborn as sns
from tabulate import tabulate
from scipy import stats
import matplotlib.pyplot as plt
from tslearn.utils import to_time_series
from sklearn.model_selection import train_test_split
from tslearn.clustering import TimeSeriesKMeans


def load_files(workdir, path_gh2ax, path_atac, path_dnase, path_h3k4me3, path_h3k4me1, path_h3k9me3, path_h3k27ac, path_h3k36me3):
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
        df = df.replace(to_replace='None', value=0)     # replace none w/ zeros
        df.iloc[:, 2:] = df.iloc[:, 2:].astype(float)   # change all numbers to float type
        yield df


def preprocessing(gh2ax_df_float):
    """ Convert gH2AX data to tslearn format

    :param gh2ax_df_float: preprocessed gH2AX dataframe from last function
    :return time-series-formatted gH2AX data, filler data of genomic coordinates
    """

    gh2ax_ts = to_time_series(gh2ax_df_float.iloc[:, 2:])
    gh2ax_ts_3d = gh2ax_ts.reshape(gh2ax_ts.shape[0], gh2ax_ts.shape[1], 1)
    gh2ax_ts_final = to_time_series(gh2ax_ts_3d)

    # The below block is filler for scikit-learn's train-test split function as we are only interested in gh2ax shape
    X_coord = np.tile(np.array(np.arange(-2, 2 + 0.005, 0.005)), (469, 1))
    X_coord_ts = to_time_series(X_coord)
    X_coord_ts_3d = X_coord_ts.reshape(X_coord_ts.shape[0], X_coord_ts.shape[1], 1)
    X_coord_ts_final = to_time_series(X_coord_ts_3d)

    return gh2ax_ts_final, X_coord_ts_final


def kmeans_DTW(gh2ax_ts_final_processed, x_coord_ts_final_processed, k, max_iterations, random_state, verbose):
    """ Dynamic time warping (DTW)-based k-means clustering of gH2AX data

        :param gh2ax_ts_final_processed: time-series-formatted gH2AX data
        :param x_coord_ts_final_processed: filler data of genomic coordinates
        :param k: k value (i.e. # of clusters)
        :param max_iterations: maximum # of iterations for k-means algorithm

    """
    indices = np.arange(gh2ax_ts_final_processed.shape[0])
    X_train, X_test, y_train, y_test, train_ind, test_ind = train_test_split(gh2ax_ts_final_processed, x_coord_ts_final_processed,
                                                                             indices, test_size=0.33, random_state=42)

    # DTW k-means for gH2AX
    print("DTW k-means for gH2AX for k = " + str(k))
    model = TimeSeriesKMeans(n_clusters=k, metric="dtw",
                             max_iter=max_iterations, random_state=random_state, verbose=verbose)
    model.fit(X_train)
    kmeans_labels_k10 = model.labels_.reshape(-1, 1)
    y_pred = model.fit_predict(kmeans_labels_k10)

    return X_train, train_ind, kmeans_labels_k10, y_pred, model


def plot(train_set, y_predict, dtw_kmeans_model, work_folder, k_val):
    """ Plot DTW-based k-means clustering results
    - Based on code from Romain Tavenard (https://tslearn.readthedocs.io/en/stable/auto_examples/clustering/plot_kmeans.html#sphx-glr-auto-examples-clustering-plot-kmeans-py)

    :param train_set: training set generated by scikit-learn
    :param y_predict: cluster predictions generated by tslearn
    :param dtw_kmeans_model: DTW k-means model generated by tslearn
    :param work_folder: absolute path to directory containing all ChIP-seq csv files
    :param k_val: k value (i.e. # of clusters)

    """

    for yi in range(k_val):
        plt.subplot(k_val, k_val, yi + 1)
        for xx in train_set[y_predict == yi]:
            plt.plot(xx.ravel(), "k-", alpha=.2)
        plt.plot(dtw_kmeans_model.cluster_centers_[yi].ravel(), "r-")
        plt.xlim(0, 803)                              # x axis: array indices in 5kb bins (+/- 2Mb of cut site)
        plt.ylim(-10, 100)                            # y axis: gH2AX RPM enrichment
        plt.text(0.55, 0.85, 'Cluster %d' % (yi + 1), transform=plt.gca().transAxes)
        if yi == 1:
            plt.title("DTW $k$-means")
    plt.savefig(work_folder + 'DTW_kmeans_gh2ax.png')
    plt.close()
    
    
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

    :return numpy arrays of average RPM enrichment per cut site for all chromatin features
    """
    avg_atac_arr = []
    avg_h3k4me1_arr = []
    avg_h3k4me3_arr = []
    avg_h3k27ac_arr = []
    avg_h3k36me3_arr = []
    avg_h3k9me3_arr = []
    avg_dnase_arr = []

    # Calculate average RPM absolute change +/- 2Mb each cut site for all chipseq data
    gh2ax_achange_avg, atac_avg, dnase_avg, ctcf_avg, smc3_avg, h3k4me1_avg, h3k4me3_avg, h3k9me3_avg, h3k27ac_avg, h3k36me3_avg = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

    for row in range(gh2ax_df.shape[0]):
        atac_avg[atac_df.iloc[row, 0] + ":" + str(atac_df.iloc[row, 1])] = np.mean(atac_df.iloc[row, 2:])
        dnase_avg[dnase_df.iloc[row, 0] + ":" + str(dnase_df.iloc[row, 1])] = np.mean(dnase_df.iloc[row, 2:])
        h3k4me1_avg[h3k4me1_df.iloc[row, 0] + ":" + str(h3k4me1_df.iloc[row, 1])] = np.mean(h3k4me1_df.iloc[row, 2:])
        h3k4me3_avg[h3k4me3_df.iloc[row, 0] + ":" + str(h3k4me3_df.iloc[row, 1])] = np.mean(h3k4me3_df.iloc[row, 2:])
        h3k9me3_avg[h3k9me3_df.iloc[row, 0] + ":" + str(h3k9me3_df.iloc[row, 1])] = np.mean(h3k9me3_df.iloc[row, 2:])
        h3k27ac_avg[h3k27ac_df.iloc[row, 0] + ":" + str(h3k27ac_df.iloc[row, 1])] = np.mean(h3k27ac_df.iloc[row, 2:])
        h3k36me3_avg[h3k36me3_df.iloc[row, 0] + ":" + str(h3k36me3_df.iloc[row, 1])] = np.mean(h3k36me3_df.iloc[row, 2:])

    for j in atac_avg.keys():                       # atac
        avg_atac_arr.append(atac_avg[j])
    for k in h3k4me1_avg.keys():                    # h3k4me1
        avg_h3k4me1_arr.append(h3k4me1_avg[k])
    for a in h3k4me3_avg.keys():                    # h3k4me3
        avg_h3k4me3_arr.append(h3k4me3_avg[a])
    for b in h3k27ac_avg.keys():                    # h3k4me3
        avg_h3k27ac_arr.append(h3k27ac_avg[b])
    for c in h3k36me3_avg.keys():                   # h3k36me3
        avg_h3k36me3_arr.append(h3k36me3_avg[c])
    for d in h3k9me3_avg.keys():                    # h3k9me3
        avg_h3k9me3_arr.append(h3k9me3_avg[d])
    for e in dnase_avg.keys():                      # dnase
        avg_dnase_arr.append(dnase_avg[e])

    # convert all to np array
    avg_atac_arr = np.array(avg_atac_arr)
    avg_h3k4me1_arr = np.array(avg_h3k4me1_arr)
    avg_h3k4me3_arr = np.array(avg_h3k4me3_arr)
    avg_h3k27ac_arr = np.array(avg_h3k27ac_arr)
    avg_h3k36me3_arr = np.array(avg_h3k36me3_arr)
    avg_h3k9me3_arr = np.array(avg_h3k9me3_arr)
    avg_dnase_arr = np.array(avg_dnase_arr)

    return avg_atac_arr, avg_h3k4me1_arr, avg_h3k4me3_arr, avg_h3k27ac_arr, avg_h3k36me3_arr, avg_h3k9me3_arr, avg_dnase_arr


def dtw_epigen(train_indices, model_kmeans, avg_atac_arr, avg_h3k4me1_arr, avg_h3k4me3_arr, avg_h3k27ac_arr, avg_h3k36me3_arr, avg_h3k9me3_arr, avg_dnase_arr):
    """ Relate DTW k-means results to epigenetics

    :param train_indices: indices (row #) of the original gH2AX dataframe, mapped to members of training set
    :param model_kmeans: DTW k-means model generated by tslearn
    :param avg_atac_arr: numpy array of average ATAC-seq RPM enrichment per cut site
    :param avg_h3k4me1_arr: numpy array of average H3K4me1 RPM enrichment per cut site
    :param avg_h3k4me3_arr: numpy array of average H3k4me3 RPM enrichment per cut site
    :param avg_h3k27ac_arr: numpy array of average H3k27ac RPM enrichment per cut site
    :param avg_h3k36me3_arr: numpy array of average H3K36me3 RPM enrichment per cut site
    :param avg_h3k9me3_arr: numpy array of average H3k9me3 RPM enrichment per cut site
    :param avg_dnase_arr: numpy array of average DNase-seq RPM enrichment per cut site

    :return: - dataframe containing average RPM enrichment for clusters 1 and 2 for all epigenetic features
             - indices (row #) of cluster 1 and 2 cut sites in original gH2AX dataframe, mapped to members of training set
             - lists containing avg RPM enrichment of epigenetic features (by cluster)

    - Note that 'cluster 1' and 'cluster 2' refer to two extremes on the gH2AX shape spectrum (i.e. gH2AX-high, gH2AX-low)
    """

    kmeans_labels_reshaped = np.ravel(model_kmeans.labels_)    # flatten array
    cluster1_ind = np.array([i for i, x in enumerate(kmeans_labels_reshaped) if x == 4])    # "cluster"/shape with lowest gh2ax enrichment
    cluster2_ind = np.array([i for i, x in enumerate(kmeans_labels_reshaped) if x == 6])    # "cluster"/shape with highest gh2ax enrichment

    train_ind_clust1 = train_indices[cluster1_ind]
    train_ind_clust2 = train_indices[cluster2_ind]

    cluster1_h3k4me3 = []
    cluster2_h3k4me3 = []
    cluster1_h3k27ac = []
    cluster2_h3k27ac = []
    cluster1_h3k36me3 = []
    cluster2_h3k36me3 = []
    cluster1_h3k4me1 = []
    cluster2_h3k4me1 = []
    cluster1_h3k9me3 = []
    cluster2_h3k9me3 = []
    cluster1_atac = []
    cluster2_atac = []
    cluster1_dnase = []
    cluster2_dnase = []

    for i in train_ind_clust1:
        cluster1_h3k4me3.append(avg_h3k4me3_arr[i])
    for j in train_ind_clust2:
        cluster2_h3k4me3.append(avg_h3k4me3_arr[j])

    for i in train_ind_clust1:
        cluster1_h3k4me1.append(avg_h3k4me1_arr[i])
    for j in train_ind_clust2:
        cluster2_h3k4me1.append(avg_h3k4me1_arr[j])

    for i in train_ind_clust1:
        cluster1_h3k27ac.append(avg_h3k27ac_arr[i])
    for j in train_ind_clust2:
        cluster2_h3k27ac.append(avg_h3k27ac_arr[j])

    for i in train_ind_clust1:
        cluster1_h3k36me3.append(avg_h3k36me3_arr[i])
    for j in train_ind_clust2:
        cluster2_h3k36me3.append(avg_h3k36me3_arr[j])

    for i in train_ind_clust1:
        cluster1_h3k9me3.append(avg_h3k9me3_arr[i])
    for j in train_ind_clust2:
        cluster2_h3k9me3.append(avg_h3k9me3_arr[j])

    for i in train_ind_clust1:
        cluster1_atac.append(avg_atac_arr[i])
    for j in train_ind_clust2:
        cluster2_atac.append(avg_atac_arr[j])

    for i in train_ind_clust1:
        cluster1_dnase.append(avg_dnase_arr[i])
    for j in train_ind_clust2:
        cluster2_dnase.append(avg_dnase_arr[j])

    cluster1_2_h3k4me3 = cluster1_h3k4me3 + cluster2_h3k4me3
    cluster1_2_h3k4me1 = cluster1_h3k4me1 + cluster2_h3k4me1
    cluster1_2_h3k36me3 = cluster1_h3k36me3 + cluster2_h3k36me3
    cluster1_2_h3k27ac = cluster1_h3k27ac + cluster2_h3k27ac
    cluster1_2_h3k9me3 = cluster1_h3k9me3 + cluster2_h3k9me3
    cluster1_2_atac = cluster1_atac + cluster2_atac
    cluster1_2_dnase = cluster1_dnase + cluster2_dnase

    cluster1_2_all_feat = [cluster1_h3k4me3, cluster2_h3k4me3, cluster1_h3k4me1, cluster2_h3k4me1,
                           cluster1_h3k36me3, cluster2_h3k36me3, cluster1_h3k27ac, cluster2_h3k27ac,
                           cluster1_h3k9me3, cluster2_h3k9me3, cluster1_atac, cluster2_atac, cluster1_dnase,
                           cluster2_dnase]

    cluster_labels = (['Cluster 1'] * len(cluster1_h3k4me3) + ['Cluster 2'] * len(cluster2_h3k4me3)) * 7

    cluster1_2_all = cluster1_2_h3k4me3 + cluster1_2_h3k4me1 + cluster1_2_h3k36me3 + cluster1_2_h3k27ac + cluster1_2_h3k9me3 \
                     + cluster1_2_atac + cluster1_2_dnase
    cluster1_2 = np.column_stack([cluster1_2_all, cluster_labels])
    cluster1_2_df = pd.DataFrame(cluster1_2, columns=['Per cut site average RPM enrichment', 'cluster'])
    marks_str = ['h3k4me3'] * int(cluster1_2_df.shape[0] / 7) + ['h3k4me1'] * int(cluster1_2_df.shape[0] / 7) \
                  + ['h3k36me3'] * int(cluster1_2_df.shape[0] / 7) + ['h3k27ac'] * int(cluster1_2_df.shape[0] / 7) + \
                  ['h3k9me3'] * int(cluster1_2_df.shape[0] / 7) + ['atac'] * int(cluster1_2_df.shape[0] / 7) + \
                ['dnase'] * int(cluster1_2_df.shape[0] / 7)

    cluster1_2_df['histone_marks/chromatin features'] = marks_str
    cluster1_2_df.iloc[:, 0] = cluster1_2_df.iloc[:, 0].astype(float)

    return cluster1_2_df, train_ind_clust1, train_ind_clust2, cluster1_2_all_feat


def boxplot(cluster1_2_df, workdir):
    """ Box plot of average epigenetic RPM enrichment for clusters/shapes 1 and 2, saved to working directory

    :param cluster1_2_df: dataframe containing average RPM enrichment for clusters 1 and 2 for all epigenetic features
    :param workdir: absolute path to directory containing all ChIP-seq csv files
    """
    sns.set_theme(style="whitegrid")
    ax_boxplt = sns.boxplot(x='histone_marks/chromatin features', y='Per cut site average RPM enrichment', hue='cluster', data=cluster1_2_df, palette="Set3")
    ax_boxplt.figure.savefig(workdir + 'boxplot_epigen.png')


def ttest(cluster1_2_all_feat):
    """ Two-sample t test of epigenetic RPM enrichment between clusters/shapes 1 and 2; results printed to screen

    :param cluster1_2_all_feat: dataframe containing average RPM enrichment for clusters 1 and 2 for all epigenetic features
    """
    h3k4me3_tteset = stats.ttest_ind(cluster1_2_all_feat[0], cluster1_2_all_feat[1])    # p = 0.01019
    h3k4me1_tteset = stats.ttest_ind(cluster1_2_all_feat[2], cluster1_2_all_feat[3])    # p = 0.00232
    h3k36me3_ttest = stats.ttest_ind(cluster1_2_all_feat[4], cluster1_2_all_feat[5])    # p = 0.001966
    h3k27ac_tteset = stats.ttest_ind(cluster1_2_all_feat[6], cluster1_2_all_feat[7])    # p = 0.001587
    h3k9me3_tteset = stats.ttest_ind(cluster1_2_all_feat[8], cluster1_2_all_feat[9])    # ns
    atac_tteset = stats.ttest_ind(cluster1_2_all_feat[10], cluster1_2_all_feat[11])     # p = 0.00029705
    dnase_tteset = stats.ttest_ind(cluster1_2_all_feat[12], cluster1_2_all_feat[13])    # p = 0.00609

    print(tabulate([['h3k4me3', str(h3k4me3_tteset[0]), str(h3k4me3_tteset[1])],
                    ['h3k4me1', str(h3k4me1_tteset[0]), str(h3k4me1_tteset[1])],
                    ['h3k36me3', str(h3k36me3_ttest[0]), str(h3k36me3_ttest[1])],
                    ['h3k27ac', str(h3k27ac_tteset[0]), str(h3k27ac_tteset[1])],
                    ['h3k9me3', str(h3k9me3_tteset[0]), str(h3k9me3_tteset[1])],
                    ['atac', str(atac_tteset[0]), str(atac_tteset[1])],
                    ['dnase', str(dnase_tteset[0]), str(dnase_tteset[1])]],
                   headers=['Chromatin feature', 't-statistic', 'p-value']))


def to_bed(gh2ax_df, train_indices_clust1, train_indices_clust2, workfolder, span_rad=50E3):
    """ Write cluster 1 and 2 genomic coordinates to BED files, saved in working directory

    :param gh2ax_df: processed dataframe of gH2AX RPM enrichment
    :param train_indices_clust1: indices (row #) of cluster 1 in original gH2AX dataframe
    :param train_indices_clust2: indices (row #) of cluster 2 in original gH2AX dataframe
    :param workfolder: absolute path to directory containing all ChIP-seq csv files
    :param span_rad: radius (bp) of region centered at cut site
    """
    with open(workfolder + 'dtw_kmeans_gh2ax/clust1.bed', 'w') as f1:
        for row in train_indices_clust1:
            chr_i_merged = gh2ax_df.iloc[row, 0]
            reg_sta = gh2ax_df.iloc[row, 1] - span_rad
            reg_end = gh2ax_df.iloc[row, 1] + span_rad
            if reg_sta < 0:
                raise ValueError('Genomic coordinates negative! Please adjust radius.')
            f1.write(chr_i_merged + '\t' + str(int(reg_sta)) + '\t' + str(int(reg_end)) + '\n')

    with open(workfolder + 'dtw_kmeans_gh2ax/clust2.bed', 'w') as f2:
        for row_i in train_indices_clust2:
            chr_i_merged = gh2ax_df.iloc[row_i, 0]
            reg_sta = gh2ax_df.iloc[row_i, 1] - span_rad
            reg_end = gh2ax_df.iloc[row_i, 1] + span_rad
            if reg_sta < 0:
                raise ValueError('Genomic coordinates negative! Please adjust radius.')
            f2.write(chr_i_merged + '\t' + str(int(reg_sta)) + '\t' + str(int(reg_end)) + '\n')
