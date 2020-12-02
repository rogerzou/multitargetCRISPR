"""
Script for generation of insulation scores
"""

cell_lines = ['GM12878', 'HMEC', 'HUVEC', 'IMR90', 'K562', 'KBM7', 'NHEK']

main_dir = "/Users/jayluo/HiC_analysis/"     # Containing directories of raw Hi-C matrices
wig_name = "all_chr_5kb"
gen_insu_scores(cell_lines, main_dir, wig_name, 250000, 5000)

