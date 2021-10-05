
import numpy as np
import pickle
from . import msa as msa
from . import mtss as m
from . import chipseq as c
import subprocess as sp
import os
from Bio import SeqIO
from scipy import stats


def get_primers_nested(gen, outfile, genome_str, savepath, ct_values, npr=3,
                       rad1=75, rad2=150, rad3=250):
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
                        f1out.write(">%s_%s_%i_%i_%i_%i_%i_%s_%s\n%s\n" % (
                                    seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i,
                                    r1out_i, r2out_i, r1out_i))
                        f2out.write(">%s_%s_%i_%i_%i_%i_%i_%s_%s\n%s\n" % (
                                    seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i,
                                    r1out_i, r2out_i, r2out_i))
                    for i, (r1inn_i, r2inn_i) in enumerate(zip(r1inn, r2inn)):
                        f1inn.write(">%s_%s_%i_%i_%i_%i_%i_%s_%s\n%s\n" % (
                                    seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i,
                                    r1inn_i, r2inn_i, r1inn_i))
                        f2inn.write(">%s_%s_%i_%i_%i_%i_%i_%s_%s\n%s\n" % (
                                    seq_i, chr_i, coo_i, tct_i, proto_i, align_i, i,
                                    r1inn_i, r2inn_i, r2inn_i))
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


def _get_primer3_primers_helper(inptprmr3, seq, num_primers, rng_min, rng_max, ps=None, mltplx=0):
    """ Write primer3 input file settings to "inptprmr3". Settings that can be changed are:
    :param num_primers:  # primer pairs that will be sought in each of the nested PCR
    :param rng_min: min range of sequence that is considered target region for PCR.
    :param rng_max: max range of sequence that is considered target region for PCR.
    :param ps: desired amplicon size.
    :param mltplx: if multiplexed PCR is needed. Sets more stringent conditions on primers.
    """
    if mltplx == 1:
        min_gc, max_gc, min_tm, max_tm = 40, 60, 58, 60
    else:
        min_gc, max_gc, min_tm, max_tm = 35, 65, 57, 63
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
        prim3file.write("PRIMER_MAX_GC=" + str(max_gc) + "\n")
        prim3file.write("PRIMER_MIN_GC=" + str(min_gc) + "\n")
        prim3file.write("PRIMER_MAX_TM=" + str(max_tm) + "\n")
        prim3file.write("PRIMER_MIN_TM=" + str(min_tm) + "\n")
        prim3file.write("PRIMER_MAX_POLY_X=4\n") if mltplx == 1 else None
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


def get_stats_primers(outfile, min_prmrs=1):
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
     2. # cut sites w/ >= min_prmrs primers passing primer3 filter
     3. # cut sites w/ >= min_prmrs primers passing primer3 + bowtie2 (uniqueness mapping) filters
     4. # cut sites w/ >= min_prmrs primers passing primer3 + bowtie2 (single alignment) filters
    """
    csv_out = open(outfile + '_stats.csv', 'w')
    msa_inn = m.load_nparray(outfile + "_inn_msa.csv")
    msa_out = m.load_nparray(outfile + "_out_msa.csv")
    msa_inn = np.hstack((msa_inn, np.ones((msa_inn.shape[0], 1))))
    msa_out = np.hstack((msa_out, np.zeros((msa_out.shape[0], 1))))
    msa_all = np.concatenate((msa_inn, msa_out), axis=0)
    msa_sorted = msa_all[np.lexsort((msa_all[:, 6].astype(int), msa_all[:, 5].astype(int),
                                     msa_all[:, 4].astype(int), -msa_all[:, 3].astype(int)))]
    proto_prev, tgt_prev, boo_prev, num_prev, prim_prev, seq_prev = -1, -1, -1, -1, -1, ""
    sumrow, numct = [], []
    maxnumtgts = int(msa_sorted[:, 3].max())
    nprms_prim3, nprms_bowti, nprms_final = [0]*maxnumtgts, [0]*maxnumtgts, [0]*maxnumtgts
    for msa_i in msa_sorted:
        seq_i, chr_i, coo_i, tnt_i = msa_i[:4]
        proto_i, tgt_i, prim_i = msa_i[4:7].astype(int)
        boo_i, num_i = msa_i[9:11].astype(int)
        if proto_i != proto_prev:
            if len(numct) > 0:
                sumrow[2] = str(len([i for i in nprms_prim3 if i >= min_prmrs]))
                sumrow[3] = str(len([i for i in nprms_bowti if i >= min_prmrs]))
                sumrow[4] = str(len([i for i in nprms_final if i >= min_prmrs]))
                csv_out.write(','.join(sumrow) + '\n')
            sumrow = [seq_i, tnt_i, 0, 0, 0]
            numct = []
            proto_prev = proto_i
            nprms_prim3, nprms_bowti, nprms_final = [0]*maxnumtgts, [0]*maxnumtgts, [0]*maxnumtgts
        if (prim_i == prim_prev) and (seq_i == seq_prev):
            nprms_prim3[tgt_i] += 1
            if boo_i == boo_prev == 1:
                nprms_bowti[tgt_i] += 1
                if num_i == num_prev == 1:
                    nprms_final[tgt_i] += 1
        numct.append(num_i)
        boo_prev, prim_prev, num_prev, seq_prev = boo_i, prim_i, num_i, seq_i
    if len(numct) > 0:
        sumrow[2] = str(len([i for i in nprms_prim3 if i >= min_prmrs]))
        sumrow[3] = str(len([i for i in nprms_bowti if i >= min_prmrs]))
        sumrow[4] = str(len([i for i in nprms_final if i >= min_prmrs]))
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
    bowtie2_msa_primers(curfile + "_out", genome_path)      # Check outer PCR primers (PCR-1)
    bowtie2_msa_primers(curfile + "_inn", genome_path)      # Check inner PCR primers (PCR-2)
    # Parse sam file to output info about PCR primers
    msa.parse_msa_sam_paired(curfile + "_out" + "_msa")     # Check outer PCR primers (PCR-1)
    msa.parse_msa_sam_paired(curfile + "_inn" + "_msa")     # Check inner PCR primers (PCR-2)
    # Get statistics for each gRNA
    get_stats_primers(curfile)


def get_nested_primers(f_inn, f_out, outfile, protospacer):
    """
    """
    p_inn, p_out = [], []
    with open(f_inn + '.csv', 'r') as msa_inn:
        for msa_inn_i in msa_inn:
            split_i = msa_inn_i.strip().split(',')
            if split_i[0] == protospacer:
                p_inn.append(split_i)
    with open(f_out + '.csv', 'r') as msa_out:
        for msa_out_i in msa_out:
            split_i = msa_out_i.strip().split(',')
            if split_i[0] == protospacer:
                p_out.append(split_i)
    np.savetxt(outfile + '_inn_%s.csv' % protospacer, np.asarray(p_inn), fmt='%s', delimiter=',')
    np.savetxt(outfile + '_out_%s.csv' % protospacer, np.asarray(p_out), fmt='%s', delimiter=',')


def lineage_ngs_fq2sam(ngsfile, genome_path, outfile, align_length=75):
    """ Align a part of paired-end amplicon NGS to genome for indel/SNV detection.
        Assuming read1 is 100bp and read2 is 200bp, use 'align_length' of read1/read2 for bowtie2
        alignment, while storing read2 in header for downstream mutation detection.
    :param ngsfile: path to input fastq. read1/read2 as ngsfile_1.fastq/ngsfile_2.fastq respectively
    :param genome_path: path to genome for bowtie2 to use
    :param outfile: path to output file as outfile.sam
    :param align_length: truncated length of paired-end read to use for alignment. Default is 75,
                         that is, 75bp from read1 and read2 are used for bowtie2 alignment.
    :output bowtie2 output as 'outfile'.sam
    """
    reads1 = SeqIO.parse(open(ngsfile + "_1.fastq"), 'fastq')
    reads2 = SeqIO.parse(open(ngsfile + "_2.fastq"), 'fastq')
    out1, out2 = [], []
    for ngs_1, ngs_2 in zip(reads1, reads2):            # parse each paired-end read
        name_1, seq_1 = ngs_1.id, str(ngs_1.seq)
        name_2, seq_2 = ngs_2.id, str(ngs_2.seq)
        ngs_1 = ngs_1[:align_length]                    # crop read1 to 'align_length' bp
        ngs_2 = ngs_2[:align_length]                    # crop read2 to 'align_length' bp
        ngs_1.id = name_1 + "_" + seq_2                 # add read2 to read1 header/id
        ngs_2.id = name_2 + "_" + seq_2                 # add read2 to read2 header/id
        ngs_1.description = ""
        ngs_2.description = ""
        out1.append(ngs_1)
        out2.append(ngs_2)
    SeqIO.write(out1, outfile + "_1_trunc.fq", "fastq") # save truncated read1 as fastq
    SeqIO.write(out2, outfile + "_2_trunc.fq", "fastq") # save truncated read2 as fastq
    # bowtie2
    sp.run(['bowtie2', '-p', '8', '--local', '--no-discordant', '-x', genome_path[:-3],
            '-1', outfile + '_1_trunc.fq', '-2', outfile + '_2_trunc.fq', '-S', outfile + '.sam'])


def lineage_ngs_sam2dict(infile, targetfile, proto, rc, win=1000):
    """ Takes bowtie2 alignment of amplicon NGS from lineage_ngs_fq2sam(), and stores each alignment
        grouped by identity of gRNA target position in a dictionary structure.
        Also outputs statistics for reads that are intact or have indels/SNVs.
    :param infile: path to input file - output of lineage_ngs_fq2sam()
    :param targetfile: path to input file - output of msa.get_target_sequences()
    :param proto: 20bp protospacer sequence
    :param rc: boolean indicating if protospacer should be reverse complement
    :param win: window to tell whether an alignment falls within a target region
    """
    proto = c.get_reverse_complement(proto) if rc else proto
    chr_tgt, pos_tgt, seq_tgt = _lineage_ngs_gen_sequences(targetfile)  # get each target site
    num_tgt = len(chr_tgt)   # determine # of target sites
    num_align, num_intact, ratio_mutated = [0]*num_tgt, [0]*num_tgt, [0]*num_tgt
    dict_ind, dict_int = {}, {}
    with open(infile + '.sam', 'r') as sam_it:
        for read in sam_it:         # read every line of SAM file
            if read[0] == '@':      # skip SAM header lines
                continue
            row = read.strip().split('\t')
            tgt = 0
            for chr_i, pos_i, seq_i in zip(chr_tgt, pos_tgt, seq_tgt) :    # parse target sites
                seq_i = c.get_reverse_complement(seq_i) if rc else seq_i
                if row[2] == chr_i and int(row[3]) - win < int(pos_i) < int(row[3]) + win:
                    key_i = "%s-%s-%s" % (chr_i, pos_i, seq_i)
                    if key_i not in dict_ind.keys() and key_i not in dict_int.keys():   # new key
                        dict_ind[key_i] = []
                        dict_int[key_i] = []
                    num_align[tgt] += 1
                    seq_i = row[0].strip().split('_')[1]
                    if seq_i.find(proto) > 0:       # found protospacer in NGS read (i.e., intact)
                        num_intact[tgt] += 1
                        dict_int[key_i].append(seq_i)
                    else:                           # either indel or SNV
                        dict_ind[key_i].append(seq_i)
                tgt += 1
    # output statistics of intact vs indel/SNV reads
    for i in range(len(num_align)):
        ratio_mutated[i] = 1-num_intact[i]/num_align[i] if (num_align[i] > 0) else 0
    stat_results = [num_intact, num_align, ratio_mutated]
    np.savetxt(infile + '_stats.csv', np.asarray(stat_results), fmt='%s', delimiter=',')
    # save dictionary structures using pickle
    pickle.dump(dict_int, open(infile + '_intact.pickle', 'wb'))
    pickle.dump(dict_ind, open(infile + '_indels.pickle', 'wb'))


def _lineage_ngs_gen_sequences(targetfile):
    """ Parses output of msa.get_target_sequences() to yield the chromosomes, position, and local
        sequences of each target site.
    :param targetfile: csv output of get_target_sequences()
    :yield chromosome, coordinate, sequence of each target site, all as strings.
    """
    data = m.load_nparray(targetfile)
    chr_tgt, coo_tgt, seq_tgt = [], [], []
    for i in range(len(data)):
        chr_i, coo_i, seq_i, sen_i, s_i = data[i]
        chr_tgt.append(chr_i)
        coo_tgt.append(coo_i)
        seq_tgt.append(s_i)
    return chr_tgt, coo_tgt, seq_tgt


def lineage_ngs_dict2csv(infile, proto, rc):
    """ Takes amplicon NGS reads (as intact or indels/SNV) from lineage_ngs_sam2dict(), grouped by
        target position, and determines the # of reads for each unique indel/SNV.
    :param infile: path to input - output of lineage_ngs_sam2dict()
    :param proto: 20bp protospacer sequence
    :param rc: boolean indicating if protospacer should be reverse complement
    :output (1) fasta file listing the mutation type and trimmed amplicon sequence with SNV/indels
            (2) csv file listing # of reads for each unique SNV/indel
    """
    dict_int = pickle.load(open(infile + '_intact.pickle', 'rb'))
    dict_ind = pickle.load(open(infile + '_indels.pickle', 'rb'))
    fasta_out = open(infile + "_mut.fasta", 'w')
    mut_dict = {}
    for i, key_i in enumerate(dict_int.keys()):         # iterate over each target site (dict key)
        ki = key_i.strip().split('-')
        chr_i, pos_i, seq_i = ki[0], ki[1], ki[2]       # get coordinate, position, and sequence
        lt_seq, rt_seq, lt_len, rt_len = _lineage_ngs_define(seq_i, proto, rc)
        # determine consensus sequence from intact amplicon NGS reads, compare to genomic sequence
        seq_c = consensus_sequence(dict_int[key_i])
        mut_c, cro_c = _lineage_ngs_mutations(seq_i, seq_c, lt_seq, rt_seq, lt_len, rt_len)
        if not cro_c:   # if intact consensus sequence differs from genomic sequence, use the former
            seq_i = seq_c
            lt_seq, rt_seq, lt_len, rt_len = _lineage_ngs_define(seq_i, proto, rc)
        mut_list = []
        for j, seq_j in enumerate(dict_ind[key_i]):     # record indel/SNV type from mutated reads
            mut_j, cro_j = _lineage_ngs_mutations(seq_i, seq_j, lt_seq, rt_seq, lt_len, rt_len)
            if cro_j:   # if valid comparison between NGS read and reference, determine mutation
                fasta_out.write(">%03i_%s_%s_%08i_%s\n%s\n" % (i, chr_i, pos_i, j, mut_j, cro_j))
                mut_list.append(mut_j)
        for k, seq_k in enumerate(dict_int[key_i]):     # record intact reads
            mut_k, cro_k = _lineage_ngs_mutations(seq_i, seq_k, lt_seq, rt_seq, lt_len, rt_len)
            if cro_k:   # if valid comparison between NGS read and reference, determine mutation
                mut_list.append(mut_k)
        mut_dict["%s_%s" % (chr_i, pos_i)] = [[x, mut_list.count(x)] for x in set(mut_list)]
    np.savetxt(infile + '_mut.csv', _lineage_ngs_dict2np(mut_dict), fmt='%s', delimiter=',')
    fasta_out.close()


def _lineage_ngs_define(seq, proto, rc):
    """ Given local sequence and protospacer sequence, define a local region of indel/SNV analysis.
        Determines 20bp flanking sequences to the left and right of cut site.
        Helper function of lineage_ngs_dict2csv().
    :param seq: local genomic sequence
    :param proto: 20bp protospacer sequence
    :param rc: boolean indicating if protospacer should be reverse complement
    :return (left flanking 20bp sequence,
             right flanking 20bp sequence,
             bp distance from left flanking sequence alignment to cut site,
             bp distance from right flanking sequence alignment to cut site)
    """
    proto = c.get_reverse_complement(proto) if rc else proto
    proto_ind = seq.find(proto)
    if proto_ind == -1:
        raise ValueError("_lineage_ngs_define(): sequence not found!")
    proto_ind = proto_ind + 3 if rc else proto_ind + 16
    lt1, lt2 = max(0, proto_ind - 40), max(0, proto_ind - 19)
    rt1, rt2 = min(len(seq), proto_ind + 20), min(len(seq), proto_ind + 41)
    lt_seq = seq[lt1:lt2]
    rt_seq = seq[rt1:rt2]
    return lt_seq, rt_seq, proto_ind - lt1, rt1 - proto_ind


def _lineage_ngs_mutations(seq_e, seq_i, lt_seq, rt_seq, lt_len, rt_len):
    """ Determines presence of SNV, insertion, or deletion in amplicon NGS sequence compared to
        ground truth sequence characterized using _lineage_ngs_define().
        Helper function of lineage_ngs_dict2csv().
    :param seq_e: string ground truth sequence
    :param seq_i: string amplicon NGS sequence
    :param lt_seq: left flanking 20bp sequence
    :param rt_seq: right flanking 20bp sequence
    :param lt_len: bp distance from left flanking sequence alignment to cut site
    :param rt_len: bp distance from right flanking sequence alignment to cut site
    :return (1) string that characterizes the mutation type of amplicon NGS sequence compared to
                ground truth of format "mut_intact_snv_insertion_deletion",
            (2) Cropped sequence using [lt_seq, rt_seq]. If None, then 1 or both lt and rt flanking
                sequences were not found, so invalid comparison.
    """
    lt_e = seq_e.find(lt_seq)
    rt_e = seq_e.find(rt_seq)
    len_e = rt_e - lt_e
    lt_i = seq_i.find(lt_seq)
    rt_i = seq_i.find(rt_seq)
    len_i = rt_i - lt_i
    intact, snv, insertion, deletion, seq_cropped = "-", "-", "-", "-", None
    if lt_i != -1 and rt_i != -1:
        seq_cropped = seq_i[lt_i + 20:rt_i + 1]
        if len_i > len_e:                                           # insertion
            insertion = seq_i[lt_i + lt_len:rt_i - rt_len]
        elif len_e > len_i:                                         # deletion
            deletion = str(len_e - len_i)
        else:
            if seq_i[lt_i + lt_len] == seq_e[lt_e + lt_len]:        # intact
                intact = seq_i[lt_i + lt_len]
            else:                                                   # SNV
                snv = seq_i[lt_i + lt_len]
    # else:
    #     print("_lineage_ngs_mutations(): sequence not found: %s" % seq_i)
    return "mut_%s_%s_%s_%s" % (intact, snv, insertion, deletion), seq_cropped


def _lineage_ngs_dict2np(mut_dict):
    """ Convert dictionary that stores histogram of each unique mutation for each target site
        into a numpy array. Helper function of lineage_ngs_dict2csv().
    :param mut_dict: dictionary generated in lineage_ngs_dict2csv()
    :return dictionary converted to numpy array
    """
    dictkeys = list(mut_dict.keys())
    dictkeys.sort()
    clen = len(dictkeys)
    rlen = max([len(mut_dict[key_i]) for key_i in dictkeys])
    mut_np = np.full((rlen + 1, clen * 2), '', dtype=object)
    for c_i, key_i in enumerate(dictkeys):
        mut_np[0, c_i * 2:c_i * 2 + 2] = key_i.split("_")
        mut_list = mut_dict[key_i]
        mut_list.sort(reverse=True)
        for r_i, mut_i in enumerate(mut_list):
            mut_np[r_i + 1, c_i * 2:c_i * 2 + 2] = mut_i
    return mut_np


def consensus_sequence(seqs):
    """ Given a list of sequences, generate consensus sequence from (1) all sequences with the most
        common length and (2) the most common nucleotide in each position.
    :param seqs: list of sequences
    :return string of consensus sequence
    """
    lens = stats.mode([len(s) for s in seqs])[0][0]
    seq_cnt = np.zeros((4, lens))
    for seq_i in seqs:
        if len(seq_i) == lens:
            for i, s in enumerate(seq_i):
                if s == 'A':
                    seq_cnt[0, i] += 1
                if s == 'C':
                    seq_cnt[1, i] += 1
                if s == 'G':
                    seq_cnt[2, i] += 1
                if s == 'T':
                    seq_cnt[3, i] += 1
    seq_cnt = np.argmax(seq_cnt, axis=0)
    seq_out = ""
    for s in seq_cnt:
        if s == 0:
            seq_out += 'A'
        if s == 1:
            seq_out += 'C'
        if s == 2:
            seq_out += 'G'
        if s == 3:
            seq_out += 'T'
    return seq_out


def lineage_ngs_np2sum(csv_list, keystr):
    """ Given output of lineage_ngs_dict2np(), recreate the csv file by including, for each target
        site, all unique mutation types across all timepoints. The output of lineage_ngs_dict2np()
        may have different sets of mutations at different time points, which make different time
        points hard to compare.
    :param csv_list: list of file paths from lineage_ngs_dict2np() outputs
    :param keystr: a string to add to the output of this function to distinguish it
    :output csv files of similar structure to the output of lineage_ngs_dict2np(), except ensuring
            that each timepoint has the same list of mutations.
    """
    np_list = [m.load_nparray(f + "_mut.csv") for f in csv_list]    # list of time points (csv)
    n_points = len(csv_list)                                        # num of time points
    n_target = int(np_list[0].shape[1] / 2)                         # num of target sites
    for i_pts in range(n_points):               # iterate through each time point (n_points)
        m_d = {}
        for j_target in range(n_target):           # iterate through each target (n_target)
            j_sta, j_end = j_target * 2, j_target * 2 + 2
            key_i = "%s_%s" % (np_list[i_pts][0, j_sta], str(np_list[i_pts][0, j_sta + 1]))
            # get all mutation types across all time points (k) for specific target (j)
            mut_types = set()
            for k_pts in range(n_points):
                for t in list(np_list[k_pts][1:, j_sta]):
                    if t != '':
                        mut_types.add(t)
            # recreate dictionary with all mutation types; those from different time points set to 0
            m_d[key_i] = [[x[0], int(x[1])] for x in np_list[i_pts][1:, j_sta:j_end] if x[0] != '']
            m_d[key_i] += [[x, 0] for x in mut_types if x not in set(np_list[i_pts][1:, j_sta])]
        np.savetxt(csv_list[i_pts] + '_mut_%s.csv' % keystr, _lineage_ngs_dict2np(m_d),
                   fmt='%s', delimiter=',')


def lineage_ngs_aggregate(csv_list, keystr, outfile):
    """ Summarize the output of lineage_ngs_np2sum() across all time points into one file.
        Compute Shannon (base 2) entropy values of mutations for each time point and target site.
    :param csv_list: list of file paths from lineage_ngs_dict2np() outputs
    :param keystr: a string used to further distinguish output of lineage_ngs_np2sum()
    :param outfile: string output path of summary csv file
    """
    np_list2 = [m.load_nparray(f + "_mut_%s.csv" % keystr) for f in csv_list]
    n_points = len(csv_list)                                        # num of time points
    n_target = int(np_list2[0].shape[1] / 2)                         # num of target sites
    summary_np = np.full((np_list2[0].shape[0], 1), '', dtype=object)
    entropy_np = np.full((1, 1), '', dtype=object)
    for i_target in range(n_target):        # for each target site,
        for j_pts in range(n_points):       # get the change in mutation distribution over time
            i = i_target * 2
            app_i = np_list2[j_pts][:, i:i + 2] if j_pts == 0 else np_list2[j_pts][:, i + 1:i + 2]
            val_i = [int(x[1]) for x in np_list2[j_pts][1:, i:i + 2] if x[0] != '']
            ent_i = stats.entropy([x/sum(val_i) for x in val_i], base=2) if sum(val_i) != 0 else 0
            ent_i = ['', ent_i] if j_pts == 0 else [ent_i]
            summary_np = np.append(summary_np, app_i, axis=1)
            entropy_np = np.append(entropy_np, ent_i)
    summary_np = np.append(summary_np, entropy_np[None, :], axis=0)
    np.savetxt(outfile, summary_np, fmt='%s', delimiter=',')


def lineage_ngs1(ngsfile, genome_path, verbose=1):
    """
    Mutation analysis script #1.
        1) Truncates paired-end fastq files and aligns to the genome using bowtie in PE mode.
        2) Checks for alignments that are in proximity to the msa of the mgRNA targets.
        3) Then checks for the presence of the protospacer in the read. If present, the
            read is considered intact. Else, it's considered mutated.
    Parameters of the script:
        ngsfile: name of the input fastq.
        genome_path: path to genome in pickle format.
        verbose: displays progress in reading the sam file.
    """
    # Define the msa of the targets: chromosome + position:
    chr_tgt = ["chr16", "chr16", "chr6", "chr1", "chr19", "chr16", "chr11", "chr7"]
    pos_tgt = [11642579, 70421720, 86513528, 68153228, 1904845, 72558952, 93681278, 37733626]
    proto = "CCAGGCTGGAGTGCAGTGCT"
    win = 2000  # window to tell whether an alignment falls within a target region
    num_tgt = 10
    num_align = [0]*num_tgt
    num_intact = [0]*num_tgt
    ratio_mutated = [0]*num_tgt
    reads1 = SeqIO.parse(open(ngsfile + "_1.fastq"), 'fastq')
    reads2 = SeqIO.parse(open(ngsfile + "_2.fastq"), 'fastq')
    out1, out2 = [], []
    for ngs_1, ngs_2 in zip(reads1, reads2):
        name_1, seq_1 = ngs_1.id, str(ngs_1.seq)
        name_2, seq_2 = ngs_2.id, str(ngs_2.seq)
        ngs_1 = ngs_1[:50]
        ngs_1.id = name_1 + "_" + seq_1
        ngs_2.id = name_2 + "_" + seq_1
        ngs_1.description = ""
        ngs_2.description = ""
        out1.append(ngs_1)
        out2.append(ngs_2)
    SeqIO.write(out1, ngsfile + "_1_trunc.fq", "fastq")
    SeqIO.write(out2, ngsfile + "_2_trunc.fq", "fastq")
    sp.run(['bowtie2', '-p', '8', '--local', '--very-sensitive',
            '--no-discordant', '-x', genome_path[:-3],
            '-1', ngsfile + '_1_trunc.fq', '-2', ngsfile + '_2_trunc.fq', '-S', ngsfile + '.sam'])
    with open(ngsfile + '.sam', 'r') as sam_it:
        cter = 0
        for read in sam_it:     # read every line of SAM file
            # counter to track progress
            if verbose and cter % 10000 == 0:
                print("lineage_ngs(): Read %i lines of SAM file." % cter)
            cter += 1
            # skip SAM header lines
            if read[0] == '@':
                continue
            row = read.strip().split('\t')
            tgt = 0
            for chr_i, pos in zip(chr_tgt, pos_tgt):
                if row[2] == chr_i and int(row[3]) - win < pos < int(row[3]) + win:
                    num_align[tgt] += 1
                    sequence = row[0].strip().split('_')
                    if sequence[1].find(proto) > 0:
                        num_intact[tgt] += 1
                tgt += 1
    for i in range(len(num_align)):
        ratio_mutated[i] = 1-num_intact[i]/num_align[i] if (num_align[i] > 0) else 0
    # return num_align, ratio_mutated
    return ratio_mutated


def qual_filt(infile, outfile, qual=30, perc=95):
    # Performs quality filter using the command "fastq_quality_filter" from FASTX-Toolkit
    filtered_file = outfile + "_1.fastq"
    inputfile = infile + "_1.fastq"
    sp.run(['fastq_quality_filter', '-q', str(qual), '-p', str(perc), '-i', inputfile,
            '-Q', '33', '-o', filtered_file])


def lineage_ngs2(ngsfile):
    """
    Mutation analysis script #2.
        1) Checks first bases and demultiplexes according to the hamming distance to the second
           round PCR primers.
        2) Then checks for the presence of the protospacer in the read. If present, the
           read is considered intact. Else, it's considered mutated.
    Parameters of the script:
        ngsfile: name of the input fastq.
    """
    proto = "CCAGGCTGGAGTGCAGTGCT"
    threshold = 2
    # 2nd round PCR primers for each of the targets are defined for demultiplexing.
    prmrs_tgts = {0: "TGGAGGTGTCATTTGCCCAG", 1: "ACAGCTTGGGGTTACACAGG",
                  2: "TGTTCAAGCAAAAGGTTCAGCT", 3: "TCCTCGGGGTTCTACCTCTC",
                  4: "CTGGGCCAAACAAGCACATT", 5: "GCTGCCTTGCCAGAGTTTTT",
                  6: "TCAGCTGTTTGATCTCAGGCA", 7: "GGGTCTAGCAAGTGGGCAAT",
                  8: "CCAACGTTGTTCAGGCACAC", 9: "CTCCTTGGCCCTTGGTTCAT"}
    num_tgt = 10
    num_align = [0]*num_tgt
    num_intact = [0]*num_tgt
    ratio_mutated = [0]*num_tgt
    # Perform quality filter.
    ngs_filtered = ngsfile + "_filt"
    qual_filt(ngsfile, ngs_filtered, qual=30, perc=95)
    # Go through filtered file to demultiplex and find the mutations in protospacer.
    reads1 = SeqIO.parse(open(ngs_filtered + "_1.fastq"), 'fastq')
    for ngs_1 in reads1:
        for i in prmrs_tgts:
            check_primer = prmrs_tgts[i]
            if hamming_distance(ngs_1[0:len(check_primer)], check_primer) < threshold:
                num_align[i] = num_align[i] + 1
                if ngs_1.seq.find(proto) > 1:
                    num_intact[i] = num_intact[i] + 1
    for i in range(len(num_align)):
        ratio_mutated[i] = 1-num_intact[i]/num_align[i] if (num_align[i] > 0) else 0
    return ratio_mutated


def lineage_ngs3(ngsfile, threshold=2):
    """
    Mutation analysis script #3.
        1) Checks first bases and demultiplexes according to the hamming distance to the second
           round PCR primers.
        2) Checks for "key" sequences: sets of 8 bases that are adjacent to the protospacer on
           either side. If the distance between them is different than 28
           (20 for proto + 8 for key), then an indel is called.
    Parameters of the script:
        ngsfile: name of the input fastq.
    """
    prmrs_tgts = {0: "TGGAGGTGTCATTTGCCCAG", 1: "ACAGCTTGGGGTTACACAGG",
                  2: "TGTTCAAGCAAAAGGTTCAGCT", 3: "TCCTCGGGGTTCTACCTCTC",
                  4: "CTGGGCCAAACAAGCACATT", 5: "GCTGCCTTGCCAGAGTTTTT",
                  6: "TCAGCTGTTTGATCTCAGGCA", 7: "GGGTCTAGCAAGTGGGCAAT",
                  8: "CCAACGTTGTTCAGGCACAC", 9: "CTCCTTGGCCCTTGGTTCAT"}
    key1 = {0: "TCTGTTGG", 1: "ACTCCAGC", 2: "TGTGTTAT", 3: "TCTGTTGT", 4: "TCTGTCAC",
            5: "TCCGTCAC", 6: "TCTGTCGC", 7: "TCTGTCAC", 8: "TCTGTTGC", 9: "TCTGTTGC"}
    key2 = {0: "GGGATCTT", 1: "AGGATCAT", 2: "TGGCTTAC", 3: "GGGATCTC", 4: "AGGCTCAC",
            5: "GGGATCTT", 6: "GGGATCTT", 7: "GGGATCTC", 8: "CGGCTCAC", 9: "GGGGGATC"}
    num_tgt = 10
    num_align = [0] * num_tgt
    num_intact = [0] * num_tgt
    ratio_mutated = [0] * num_tgt
    reads_merged = SeqIO.parse(open(ngsfile + "_1.fastq"), 'fastq')
    for ngs_1 in reads_merged:
        for i in prmrs_tgts:
            check_primer = prmrs_tgts[i]
            if hamming_distance(ngs_1[0:len(check_primer)], check_primer) < threshold:
                if ngs_1.seq.find(key1[i]) > 0 and ngs_1.seq.find(key2[i]) > 0:
                    num_align[i] = num_align[i] + 1
                    if (ngs_1.seq.find(key2[i])-ngs_1.seq.find(key1[i])) == 28:
                        num_intact[i] = num_intact[i] + 1
    for i in range(len(num_align)):
        ratio_mutated[i] = 1-num_intact[i]/num_align[i] if (num_align[i] > 0) else 0
    return ratio_mutated


def hamming_distance(string1, string2):
    if len(string1) != len(string2):
        raise ValueError('Strings should be the same length')
    else:
        return sum(c1 != c2 for c1, c2 in zip(string1, string2))
