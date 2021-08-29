
import numpy as np
from . import msa as msa
from . import mtss as m
import subprocess as sp
import os


def get_primers_nested(gen, outfile, genome_str, savepath, ct_values, npr=3,
                       rad1=100, rad2=200, rad3=350):
    """
    :param gen: bowtie2 alignments for each putative protospacer sequence of format:
        - boolean: whether it is a new putative protospacer sequence
        - string: sequence of putative protospacer sequence
        - string: '+' or '-' that correspond to the sense of current protospacer sequence alignment
        - string: chromosome string (e.g. 'chr1') of current protospacer sequence alignment
        - int: coordinate of current protospacer sequence alignment
        - int: total number of alignments for current protospacer sequence
    :param outfile: string path to output file (extension omitted). Four fasta files are generated:
        - outfile + "_out_1/2.fa" for the outer PCR primers (PCR-1)
        - outfile + "_in_1/2.fa" for the inner PCR primers (PCR-2)
    :param genome_str: 'hg38', 'hg19', or 'mm10'
    :param savepath: save path to genome sequence that will be downloaded from NCBI if not already
    :param ct_values: restricts PCR primer generation only for protospacers with # of target sites
        included in this list of integers
    :param npr: # primer pairs that will be sought in each of the nested PCR, e.g. npr=3 will search
        for 3 for/rev primer pairs in PCR-1 and 3 for/rev primer pairs in PCR-2.
    :param rad1: max radius of inner amplicon range
    :param rad2: max radius of inner primer search range, bounded inside by rad1
    :param rad3: max radius of outer primer search range, bounded inside by rad2
    Inner nested primers are located rad1-rad2 from the cut site.
    Outer nested primers are located rad3-rad2 from the cut site.
    """
    proto_i, align_i = -1, -1
    msa.genome_initialized(savepath, genome_str)
    cter = 0
    with open(outfile + "_inn_1.fa", 'w') as f1inn, open(outfile + "_inn_2.fa", 'w') as f2inn, \
            open(outfile + "_out_1.fa", 'w') as f1out, open(outfile + "_out_2.fa", 'w') as f2out:
        for new_i, seq_i, sen_i, chr_i, coo_i, tct_i in gen:
            if tct_i in ct_values:
                proto_i = proto_i + 1 if new_i else proto_i
                align_i = 0 if new_i else align_i + 1
                if chr_i in msa.CHR:
                    seq = msa.get_target_sequence(chr_i, coo_i, sen_i, rad3)  # get sequence
                    rmin_out, rmax_out = rad3 + 1 - rad2, rad3 + 1 + rad2     # index of outer range
                    rmin_inn, rmax_inn = rad2 + 1 - rad1, rad2 + 1 + rad1     # index of inner range
                    r1out, r2out = get_primer3_primers(seq, num_primers=npr,
                                                       rng_min=rmin_out,
                                                       rng_max=rmax_out)
                    r1inn, r2inn = get_primer3_primers(seq[rmin_out:rmax_out], num_primers=npr,
                                                       rng_min=rmin_inn,
                                                       rng_max=rmax_inn)
                    for i, (r1out_i, r2out_i) in enumerate(zip(r1out, r2out)):
                        f1out.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" %
                                    (seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i, r1out_i))
                        f2out.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" %
                                    (seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i, r2out_i))
                    for i, (r1inn_i, r2inn_i) in enumerate(zip(r1inn, r2inn)):
                        f1inn.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" %
                                    (seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i, r1inn_i))
                        f2inn.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" %
                                    (seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i, r2inn_i))
            if cter % 10000 == 0:
                print("get_primers_nested(): processed %i samples" % cter)
            cter += 1


def get_primer3_primers(seq, num_primers, rng_min, rng_max, prod_size=None):
    """ Finds primers using primer3
    :param seq: template sequence to find primers
    :param num_primers:  # primer pairs that will be sought in each of the nested PCR
    :param rng_min: min range of sequence that is considered target region for PCR.
    :param rng_max: max range of sequence that is considered target region for PCR.
    :param prod_size: desired amplicon size.
    """
    inptprmr3, outprmr3 = "_tmp_primer3_input.dat", "_tmp_primer3_output.dat"
    _get_primer3_primers_helper(inptprmr3, seq, num_primers, rng_min, rng_max, prod_size)
    cmnd_prmr3 = "primer3_core " + inptprmr3 + "> " + outprmr3
    os.system(cmnd_prmr3)
    primers_list, primers_rt, primers_lt = [], [], []
    cter_l, cter_r = 0, 0
    with open(outprmr3, 'r') as prmrs_file:
        for line in prmrs_file:
            row = line.split('=')
            if row[0] == "PRIMER_LEFT_" + str(cter_l) + "_SEQUENCE":
                primers_lt.append(row[1].strip())
                cter_l += 1
            if row[0] == "PRIMER_RIGHT_" + str(cter_r) + "_SEQUENCE":
                primers_rt.append(row[1].strip())
                cter_r += 1
    # clean up!
    os.system("rm " + inptprmr3 + " " + outprmr3)
    return primers_lt, primers_rt


def _get_primer3_primers_helper(inptprmr3, seq, num_primers, rng_min, rng_max, ps=None):
    """ Write primer3 input file settings to "inptprmr3". Settings that can be changed are:
    :param num_primers:  # primer pairs that will be sought in each of the nested PCR
    :param rng_min: min range of sequence that is considered target region for PCR.
    :param rng_max: max range of sequence that is considered target region for PCR.
    :param ps: desired amplicon size.
    """
    # Write sequence target (primers will bind outside of this)
    with open(inptprmr3, 'w') as prim3file:
        prim3file.write("%s%s\n" % ("SEQUENCE_ID=Proto-", 0))
        prim3file.write("%s%s\n" % ("SEQUENCE_TEMPLATE=", seq))
        prim3_trgt = "SEQUENCE_TARGET=" + str(rng_min) + "," + str(rng_max-rng_min) + "\n"
        prim3file.write(prim3_trgt)
        prim3file.write("PRIMER_TASK=generic\n")            # generic: default option for primer3
        prim3file.write("PRIMER_PICK_LEFT_PRIMER=1\n")
        prim3_return = "PRIMER_NUM_RETURN=" + str(num_primers) + "\n"
        prim3file.write(prim3_return)                       # choose how many primer sets we want
        prim3file.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")   # we don't need internal oligos
        prim3file.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
        prim3file.write("PRIMER_OPT_SIZE=20\n")
        prim3file.write("PRIMER_MIN_SIZE=18\n")
        prim3file.write("PRIMER_MAX_SIZE=22\n")
        prim3file.write("PRIMER_MAX_GC=65\n")
        prim3file.write("PRIMER_MIN_GC=35\n")
        prim3file.write("PRIMER_MIN_TM=57\n")
        prim3file.write("PRIMER_MAX_TM=63\n")
        prim3file.write("PRIMER_PRODUCT_SIZE_RANGE=100-1000\n") if ps is None else prim3file.write(
            "PRIMER_PRODUCT_SIZE_RANGE=%s\n" % ps)
        prim3file.write("PRIMER_EXPLAIN_FLAG=1\n")
        prim3file.write("=\n")


def bowtie2_msa_primers(curfile, genome_path, k_count=10):
    """ Run bowtie2 to align paired-end reads in '-k' mode, followed by samtools to sort and index
        the alignment.

    :param curfile: path to base file name/path.
            The PE read files are of format: curfile + '_1.fa' | curfile + '_1.fa'.
            The output initial SAM file is of format: curfile + '_msa.sam'
            The output (unsorted) BAM file is of format: curfile + '_msa.bam'
            The output sorted BAM file is of format: curfile + '_msa_sorted.bam'
            The output BAM index file is of format: curfile + '_msa_sorted.bam.bai'
    :param genome_path: path to genome for bowtie2 to use
    :param k_count: number of distinct, valid alignments for each PE read
    """
    sp.run(['bowtie2', '-f', '-p', '8', '--local', '--score-min', 'L,0,1.5',
            '-k', str(k_count), '-X', '1000', '--no-mixed', '--no-discordant', '-L', '18',
            '-N', '1', '-x', genome_path[:-3],
            '-1', curfile + '_1.fa', '-2', curfile + '_2.fa', '-S', curfile + '_msa.sam'])
    sp.run(['samtools', 'view', '-h', '-S', '-b', '-o', curfile + '_msa.bam', curfile + '_msa.sam'])
    sp.run(['samtools', 'sort', '-o', curfile + '_msa_sorted.bam', curfile + '_msa.bam'])
    sp.run(['samtools', 'index', curfile + '_msa_sorted.bam'])


def get_stats_primers(outfile, min_prmrs=2):
    """ Given the output of 'parse_msa_sam_paired()' determine, for each gRNA:
        (1) # targets there are for that gRNA
        (2) # targets passed the primer3 filter, with at least min_prmrs number of primers
        (3) # targets passed the primer3 and bowtie2 filters, with at least min_prmrs # of primers
            when applying the bowtie2 filter, primers with more than one alignment are discarded,
            even if one of the alignments has much better score than the others
    :param outfile: path of output of 'parse_msa_sam_paired_*()', extension omitted
    :param min_prmrs: min # primer sets found in order to consider a target amplifiable

    Outputs csv file with following columns:
     0. original gRNA sequence
     1. total number of cut sites for gRNA sequence
     2. # cut sites w/ >= min_prmrs primers that pass primer3 filter
     3. # cut sites w/ >= min_prmrs primers that pass primer3 + bowtie2 (uniqueness mapping) filters
    """
    csv_out = open(outfile + '_stats.csv', 'w')
    msa_inn = m.load_nparray(outfile + "_inn_msa.csv")
    msa_out = m.load_nparray(outfile + "_out_msa.csv")
    msa_inn = np.hstack((msa_inn, np.ones((msa_inn.shape[0], 1))))
    msa_out = np.hstack((msa_out, np.zeros((msa_out.shape[0], 1))))
    msa_all = np.concatenate((msa_inn, msa_out), axis=0)
    msa_sorted = msa_all[np.lexsort((msa_all[:, 6].astype(int), msa_all[:, 5].astype(int),
                                     msa_all[:, 4].astype(int), -msa_all[:, 3].astype(int)))]
    proto_prev, tgt_prev, boo_prev, num_prev, prim_prev, seq_prev = None, None, None, None, None, None
    sumrow, numct = [], []
    max_num_tgts = int(msa_sorted[:, 3].max())
    num_prms_prim3, num_prms_final = [0]*max_num_tgts, [0]*max_num_tgts
    for msa_i in msa_sorted:
        seq_i, chr_i, coo_i, tnt_i = msa_i[:4]
        proto_i, tgt_i, prim_i, boo_i, num_i = msa_i[4:9].astype(int)
        if proto_i != proto_prev:
            if len(numct) > 0:
                sumrow[2] = str(len([i for i in num_prms_prim3 if i >= min_prmrs]))
                sumrow[3] = str(len([i for i in num_prms_final if i >= min_prmrs]))
                csv_out.write(','.join(sumrow) + '\n')
            sumrow = [seq_i, tnt_i, 0, 0]
            numct = []
            proto_prev = proto_i
            num_prms_prim3, num_prms_final = [0]*max_num_tgts, [0]*max_num_tgts
        if (prim_i == prim_prev) and (seq_i == seq_prev):
            num_prms_prim3[tgt_i] += 1
            if (int(boo_i) == int(boo_prev) == 1) and (int(num_i) == int(num_prev) == 1):
                num_prms_final[tgt_i] += 1
        numct.append(num_i)
        boo_prev, prim_prev, num_prev, seq_prev = boo_i, prim_i, num_i, seq_i
    if len(numct) > 0:
        sumrow[2] = str(len([i for i in num_prms_prim3 if i >= min_prmrs]))
        sumrow[3] = str(len([i for i in num_prms_final if i >= min_prmrs]))
        csv_out.write(','.join(sumrow) + '\n')
    csv_out.close()


def bowtie_parse_stats_wrapper(curfile, genome_path):
    """ Wrapper function that runs bowtie to align primer pairs to genome, parses SAM output,
        then obtains statistics on the best gRNAs.
        :param curfile: path to base file name/path. *** corresponds to either 'inn' or 'out'
            Initial SAM files are of format: curfile + '_***_msa.sam'
            Unsorted BAM files are of format: curfile + '_***_msa.bam'
            Sorted BAM files are of format: curfile + '_***_msa_sorted.bam'
            BAM index files are of format: curfile + '_***_msa_sorted.bam.bai'
            CSV files from parsing BAM are of format: curfile + '_***_msa.csv'
            CSV files of gRNA statistics are of format: curfile + '_stats.csv'
        :param genome_path: path to genome for bowtie2 to use.
    """
    # Run bowtie2 to align pairs of primers to the genome in order to check for uniqueness.
    bowtie2_msa_primers(curfile + "_out", genome_path)  # Check outer PCR primers (PCR-1)
    bowtie2_msa_primers(curfile + "_inn", genome_path)  # Check inner PCR primers (PCR-2)
    # Parse sam file to output info about PCR primers
    msa.parse_msa_sam_paired(curfile + "_out" + "_msa") # Check outer PCR primers (PCR-1)
    msa.parse_msa_sam_paired(curfile + "_inn" + "_msa") # Check inner PCR primers (PCR-2)
    # Get statistics for each gRNA
    get_stats_primers(curfile)
