"""
Script for analysis of gH2AX-high and gH2AX-low cut sites
"""
import src.gh2ax_shapes as g

work_dir = '/Users/jayluo/HiC_analysis/gH2AX_53BP1_files/'
gh2ax_path = 'ALL_gh2ax_WT-3h_hg19_achange_merged.csv'
atac_path = 'epigenetics_chipseq_hg19/ALL_atacseq_bamvals_merged.csv'
dnase_path = 'epigenetics_chipseq_hg19/ALL_dnaseseq_bamvals_merged.csv'
h3k4me3_path = 'epigenetics_chipseq_hg19/ALL_h3k4me3_bamvals_merged.csv'
h3k4me1_path = 'epigenetics_chipseq_hg19/ALL_h3k4me1_bamvals_merged.csv'
h3k9me3_path = 'epigenetics_chipseq_hg19/ALL_h3k9me3_bamvals_merged.csv'
h3k27ac_path = 'epigenetics_chipseq_hg19/ALL_h3k27ac_bamvals_merged.csv'
h3k36me3_path = 'epigenetics_chipseq_hg19/ALL_h3k36me3_bamvals_merged.csv'
mre11_path = 'epigenetics_chipseq_hg19/ALL_mre11_hg19_1250_rc.csv'

gh2ax_df_raw, atac_df_raw, dnase_df_raw, h3k4me3_df_raw, h3k4me1_df_raw, h3k9me3_df_raw, h3k27ac_df_raw, h3k36me3_df_raw \
    = g.load_files(work_dir, gh2ax_path, atac_path, dnase_path, h3k4me3_path, h3k4me1_path, h3k9me3_path, h3k27ac_path,
                 h3k36me3_path)
                 
gh2ax_float, atac_float, dnase_float, h3k4me3_float, h3k4me1_float, h3k9me3_float, h3k27ac_float, h3k36me3_float = \
    g.process_files(gh2ax_df_raw, atac_df_raw, dnase_df_raw, h3k4me3_df_raw, h3k4me1_df_raw, h3k9me3_df_raw,
                  h3k27ac_df_raw, h3k36me3_df_raw)
                  
top10pct_ind, bot10pct_ind = g.gh2ax_enrichment(gh2ax_float)

avg_atac_array, avg_h3k4me1_array, avg_h3k4me3_array, avg_h3k27ac_array, avg_h3k36me3_array, avg_h3k9me3_array, avg_dnase_array = \
    g.avg_chipseq(gh2ax_float, atac_float, dnase_float, h3k4me3_float, h3k4me1_float, h3k9me3_float, h3k27ac_float,
                h3k36me3_float)

shapes_df, shape1_2_all_features = g.epigen(top10pct_ind, bot10pct_ind, avg_atac_array,
                                          avg_h3k4me1_array, avg_h3k4me3_array, avg_h3k27ac_array,
                                          avg_h3k36me3_array, avg_h3k9me3_array, avg_dnase_array)

g.boxplot(shapes_df, work_dir)
g.mre11_groups(top10pct_ind, bot10pct_ind, work_dir, mre11_path)
g.ttest(shape1_2_all_features)
g.to_bed(gh2ax_float, top10pct_ind, bot10pct_ind, work_dir, span_rad=50E3)


