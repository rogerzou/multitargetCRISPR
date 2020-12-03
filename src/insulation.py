"""
Functions for generation of insulation scores from raw normalized Hi-C matrices
"""

import numpy as np
import numpy.ma as ma
import os
from scipy.sparse import csr_matrix

CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX']


def load_matrices(i, matrix_path):
    """ Import normalized raw, binned (e.g. 5kb, 10kb, 50kb) Hi-C contact matrices (in txt format)
        - Bin sizes must be multiples of 5000
        - Modify path to raw matrix as needed
        - NaN values are replaced with zeros
        - Normalized raw matrices were generated from Juicer's Knight-Ruiz (KR) normalization
        (see supplemental information in Rao et al. 2014)

        :param i: [string] chromosome name (e.g. 'chr7')
        :param matrix_path: path to raw normalized Hi-C matrix
        """
    if os.stat(matrix_path).st_size == 0:
        raise ValueError(i + " raw matrix is empty! Juicer's KR normalization may have failed to converge.")
    raw_matrix_files = np.genfromtxt(matrix_path, missing_values='NaN')
    NaN_contacts = np.isnan(raw_matrix_files)
    raw_matrix_files[NaN_contacts] = 0
    raw_matrix_file = raw_matrix_files.astype('int32')
    return raw_matrix_file


def reformat_raw_matrices(raw_matrix_files):
    """ Convert txt format to sparse matrix format. Use csr format (compressed sparse row; triplet format).

        :param raw_matrix_files: [txt] raw Hi-C matrices
    """
    chr_raw_matrix_col_0 = np.array(raw_matrix_files[:, 0])
    chr_raw_matrix_col_1 = np.array(raw_matrix_files[:, 1])
    chr_raw_matrix_col_2 = np.array(raw_matrix_files[:, 2])
    output_matrix = csr_matrix(
        (chr_raw_matrix_col_2,
         (chr_raw_matrix_col_0, chr_raw_matrix_col_1)))
    return output_matrix


def generate_insulation_scores(raw_matrices, large=250000, small=5000):
    """ Aggregate signal within each large sliding square (50 bins x 50 bins submatrix).
        - The large square is slid along the matrix diagonal with an increment of the bin size of the input Hi-C data.
        - The mean # of contacts per large sliding square is assigned to each small square ((bin size) bp x (bin size) bp)
        (e.g. for 5kb-binned Hi-C data, assign mean signal of 250kb x 250kb sliding squares to each 5kb x 5kb square)
        - The first and last (large) bp of the matrix is not assigned an insulation score as the large sliding square
        would exceed matrix dimensions.
        - This algorithm is adapted from Crane et al. 2015 Nature

        :param raw_matrices: normalized raw matrices in csr format
        :param large: [integer] size (bp) of large sliding window
        :param small: [integer] bin size of Hi-C data (resolution) (bp) {5000, 10000, 50000}
                                bin sizes must be multiples of 5000.
    """

    signal_5kb = []
    chr_size = raw_matrices.shape[0]
    for n in range(large, chr_size - large, small):
        sliding_window = raw_matrices[(n - large):(n + small), (n + small):(n + 2 * small + large)]
        sliding_window_sum = np.sum(sliding_window)  # Aggregate signal within sliding square
        sliding_window_avg = sliding_window_sum / (large / small)  # Calculate mean signal per sliding square
        signal_5kb.append(sliding_window_avg)  # Assign mean signal of large sliding square to 5kb bin
    avg_contacts = sum(signal_5kb) / len(signal_5kb)  # Calculate avg # of contacts per chromosome
    # Calculate log2 ratio between the sum of contacts per 5kb bin and the chromosome average
    log2_normalized_contacts = [ma.log2(i_5kb / avg_contacts) for i_5kb in signal_5kb]
    normalized_contacts = np.array(log2_normalized_contacts)  # Convert list of normalized scores to array
    NaN_contacts_normalized = np.isnan(normalized_contacts)
    normalized_contacts[NaN_contacts_normalized] = 0  # Set NaN entries to zero

    return normalized_contacts


def convert_to_wiggle(normalized_scores, wig_i, wigout_path, wigname, sampname, large=250000, small=5000):
    """ Normalized score arrays are converted into a merged wiggle file for ease of viewing

    :param normalized_scores: [array] normalized insulation scores
    :param wig_i: [string] chromosome name (e.g. 'chr7')
    :param wigout_path: wiggle file output path
    :param wigname: name of wiggle file
    :param sampname: name of sample
    :param large: [integer] size (bp) of large sliding window
    :param small: [integer] bin size (resolution) (bp) of Hi-C matrix
                            bin sizes must be multiples of 5000.
    """

    with open(wigout_path + wigname + '_' + sampname + ".wig", 'a') as wig:
        wig.write('track type=wiggle_0\nfixedStep' + ' chrom=' + str(wig_i) + ' start=' + str(large) +
                  ' step=' + str(small) + '\n')
        np.savetxt(wig, normalized_scores, delimiter=' ')
    wig.close()


def gen_insu_scores(samples, maindir, wiggle_name, slid_sq_size, binsize):
    """ Calculate insulation scores for all chromosomes
    :param samples: list of samples or cell lines to be processed
    :param maindir: folder containing all sample folders and folders containing raw normalized Hi-C matrices within sample folders

           maindir/
           |--GM12878
                |--KR_normalized_raw_matrices
           |--K562
                |--KR_normalized_raw_matrices
            etc.

    :param wiggle_name: name of wiggle file
    :param [integer] slid_sq_size: size of sliding square (bp) {250000, 500000}
    :param binsize: [integer] bin size (resolution) of raw normalized Hi-C matrices (bp) {5000, 10000, 25000, 50000}
                              bin sizes must be multiples of 5000.
    """
    for sample in samples:
        for chr_i in CHR:
            raw_matrices = maindir + "Rao_2014/" + str(sample) + "/KR_normalized_raw_matrices/" + str(chr_i) + "_5kb.txt"
            # Load the normalized raw matrices
            print("Loading " + chr_i + " raw matrix for " + sample + "...")
            raw_matrix_files_int = load_matrices(chr_i, raw_matrices)
            # Reformat matrices to csr matrices
            print("Reformatting " + chr_i + " raw matrix for " + sample + "...")
            reformatted_raw_matrices = reformat_raw_matrices(raw_matrix_files_int)
            # Generate insulation scores (large sliding square and bin size can be set here)
            print("Generating " + chr_i + " insulation scores for " + sample + "...")
            normalized_scores_array = generate_insulation_scores(reformatted_raw_matrices, large=slid_sq_size, small=binsize)
            # Convert to wiggle file
            print("Converting " + chr_i + " insulation scores to wiggle file for " + sample + "...")
            convert_to_wiggle(normalized_scores_array, chr_i, main_dir, wiggle_name, sample, large=slid_sq_size, small=binsize)
            # Print progress
            print("Done processing %s of %s" % (chr_i, sample))
            if chr_i == CHR[-1]:
                print("Done processing " + sample + "!")
        if sample == sample[-1]:
            print("Done processing all samples!")
