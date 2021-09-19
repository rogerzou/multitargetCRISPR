
import numpy as np
from . import msa as msa
from . import mtss as m
import subprocess as sp
import os
from Bio import SeqIO


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


def lineage_ngs(ngsfile, genome_path):
    """
    Assuming read1 is 250bp, read2 is 50bp.
    """
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

    sp.run(['bowtie2', '-p', '8', '--local', '--no-discordant', '-x', genome_path[:-3],
            '-1', ngsfile + '_1_trunc.fq', '-2', ngsfile + '_2_trunc.fq', '-S', ngsfile + '.sam'])


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
