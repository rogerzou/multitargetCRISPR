"""
Script for dynamic time warping (DTW) and k-means clustering of gH2AX data
"""
import src.dtw_kmeans_gh2ax as d

work_dir = '/Users/jayluo/HiC_analysis/gH2AX_53BP1_files/'
gh2ax_path = 'ALL_gh2ax_WT-3h_hg19_achange_merged.csv'
atac_path = 'epigenetics_chipseq_hg19/ALL_atacseq_bamvals_merged.csv'
dnase_path = 'epigenetics_chipseq_hg19/ALL_dnaseseq_bamvals_merged.csv'
h3k4me3_path = 'epigenetics_chipseq_hg19/ALL_h3k4me3_bamvals_merged.csv'
h3k4me1_path = 'epigenetics_chipseq_hg19/ALL_h3k4me1_bamvals_merged.csv'
h3k9me3_path = 'epigenetics_chipseq_hg19/ALL_h3k9me3_bamvals_merged.csv'
h3k27ac_path = 'epigenetics_chipseq_hg19/ALL_h3k27ac_bamvals_merged.csv'
h3k36me3_path = 'epigenetics_chipseq_hg19/ALL_h3k36me3_bamvals_merged.csv'

gh2ax_df_raw, atac_df_raw, dnase_df_raw, h3k4me3_df_raw, h3k4me1_df_raw, h3k9me3_df_raw, h3k27ac_df_raw, h3k36me3_df_raw\
    = d.load_files(work_dir, gh2ax_path, atac_path, dnase_path, h3k4me3_path, h3k4me1_path, h3k9me3_path, h3k27ac_path, h3k36me3_path)

gh2ax_float, atac_float, dnase_float, h3k4me3_float, h3k4me1_float, h3k9me3_float, h3k27ac_float, h3k36me3_float = \
    d.process_files(gh2ax_df_raw, atac_df_raw, dnase_df_raw, h3k4me3_df_raw, h3k4me1_df_raw, h3k9me3_df_raw, h3k27ac_df_raw, h3k36me3_df_raw)

gh2ax_ts_final_df, X_coord_ts_final_df = d.preprocessing(gh2ax_float)

training_set, training_set_indices, kmeans_labels, y_predicted, kmeans_model = d.kmeans_DTW(gh2ax_ts_final_df, X_coord_ts_final_df,
                                                                                          k=10, max_iterations=10, random_state=42, verbose=1)
d.plot(training_set, y_predicted, kmeans_model, work_dir, k_val=10)

avg_atac_array, avg_h3k4me1_array, avg_h3k4me3_array, avg_h3k27ac_array, avg_h3k36me3_array, avg_h3k9me3_array, avg_dnase_array = \
    d.avg_chipseq(gh2ax_float, atac_float, dnase_float, h3k4me3_float, h3k4me1_float, h3k9me3_float, h3k27ac_float, h3k36me3_float)

clusters_df, train_ind_1, train_ind_2, cluster1_2_all_features = d.dtw_epigen(training_set_indices, kmeans_model, avg_atac_array,
                                                                   avg_h3k4me1_array, avg_h3k4me3_array, avg_h3k27ac_array,
                                                                   avg_h3k36me3_array, avg_h3k9me3_array, avg_dnase_array)

d.boxplot(clusters_df, work_dir)
d.ttest(cluster1_2_all_features)
d.to_bed(gh2ax_float, train_ind_1, train_ind_2, work_dir, span_rad=50E3)

