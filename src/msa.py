
import pysam
import numpy as np
from Bio import Entrez, SeqIO
from . import chipseq as c
from . import mtss as m
import pickle
import subprocess as sp
import statistics
import re
import os
import csv

hg38size = None
hg38id = None
hg38seq = None
CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
       'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
       'chr22', 'chrX', 'chrY']
CHROMHMM = ['1_TssA', '2_TssAFlnk', '3_TxFlnk', '4_Tx', '5_TxWk', '6_EnhG', '7_Enh', '8_ZNF/Rpts',
            '9_Het', '10_TssBiv', '11_BivFlnk', '12_EnhBiv', '13_ReprPC', '14_ReprPCWk', '15_Quies']
fullflags_paired = ['83', '163', '99', '147', '339', '419', '355', '403', '77', '141']
initflags_paired = ['83', '99', '339', '355']
primflags_paired = ['83', '99']
fullflags_single = ['0', '16', '256', '272', '4']
initflags_single = ['0', '16', '256', '272']
primflags_single = ['0', '16']


def hg38_seq(savep=None):
    """ Return dict that holds entire sequence for each hg38 chromosome. Stores data with pickle to
        avoid multiple retrievals from NCBI.
    :param savep: save path of hg38 sequence using pickle. If not provided, will be saved to the
                  'lib' folder of this project
    :return: dict with keys as chromosomes, values as SeqIO sequence of each chromosome key
    """
    dirname, filename = os.path.split(os.path.abspath(__file__))
    mfile = savep + "/hg38.pickle" if savep else os.path.dirname(dirname) + "/lib/hg38.pickle"
    if os.path.isfile(mfile):
        print("hg38_seq(): Loading existing hg38 sequence dict from %s." % mfile)
        return pickle.load(open(mfile, 'rb'))
    else:
        print("hg38_seq(): Downloading hg38 sequence from NCBI, saving to %s." % mfile)
        d = {}
        hg38ids = hg38_ids()
        Entrez.email = "bob@example.org"
        for chr_i in hg38ids.keys():
            handle = Entrez.efetch(db="nucleotide", id=hg38ids[chr_i], rettype="fasta", strand=1)
            read_i = SeqIO.read(handle, "fasta")
            d[chr_i] = read_i
        pickle.dump(d, open(mfile, 'wb'))
        return d


def hg38_ids():
    """ Return dict that holds the NCBI GI number for each hg38 chromosome
    :return: dict with keys as chromosomes, values as NCBI GI number of each chromosome key
    """
    d = {}
    dirname, filename = os.path.split(os.path.abspath(__file__))
    with open(os.path.dirname(dirname) + "/lib/hg38_entrez.csv", 'r') as f:
        for row in csv.reader(f, delimiter=','):
            d[row[0]] = row[2].split(":")[1]
    return d


def hg38_initialized(savepath):
    """ Initialize downloading of hg38 GI indices, full genome, and sizes as global variables
    :param savepath: path to save the hg38 genome
    """
    global hg38seq, hg38id, hg38size
    if not hg38seq:
        hg38seq = hg38_seq(savepath)
    if not hg38id:
        hg38id = hg38_ids()
    if not hg38size:
        hg38size = c.hg_dict('hg38')


def get_targets_fasta(outfile, seqstr, numbases):
    """ Local exhaustive search for all possible protospacer sequences, up to 3 mismatches
    at the PAM-proximal side, from an input sequence.

    :param outfile: name of output fasta file (with .fa extension)
    :param seqstr: input sequence as a string
    :param numbases: number of bases counting from PAM-proximal side available to mutate

    """
    seqstr = seqstr.upper()
    unique_ele = set()
    with open(outfile + ".fa", 'w') as f:       # only write unique sequences to file
        # check all 20bp substrings from input sequence
        for i in range(len(seqstr)-19):
            seq_i = seqstr[i:i+20]
            for seqmod in _get_targets_fasta_helper(seq_i, numbases):
                if seqmod not in unique_ele:
                    f.write(">%s\n%s\n" % (seqmod, seqmod))
                    unique_ele.add(seqmod)
        # check all 20bp substrings from reverse complement of input sequence
        seqstr = c.get_reverse_complement(seqstr)
        for i in range(len(seqstr)-19):
            seq_i = seqstr[i:i+20]
            for seqmod in _get_targets_fasta_helper(seq_i, numbases):
                if seqmod not in unique_ele:
                    f.write(">%s\n%s\n" % (seqmod, seqmod))
                    unique_ele.add(seqmod)


def _get_targets_fasta_helper(init_str, numbases):
    """ Helper function to yield all sequences with up to 3 mismatches from template.
    Uniqueness is not checked, but the GC content is restricted to between 0.4 and 0.7

    :param init_str: template sequence string
    :param numbases: number of bases counting from PAM-proximal side available to mutate

    """
    init_list = list(init_str)
    for i in range(1, numbases - 1):
        for j in range(i + 1, numbases):
            for k in range(j + 1, numbases + 1):
                for nb_i in {'A', 'C', 'G', 'T'}:
                    for nb_j in {'A', 'C', 'G', 'T'}:
                        for nb_k in {'A', 'C', 'G', 'T'}:
                            mod_ijk = list(init_list)
                            mod_ijk[-i] = nb_i
                            mod_ijk[-j] = nb_j
                            mod_ijk[-k] = nb_k
                            seq_str = "".join(mod_ijk)
                            if 0.4 < get_gc(seq_str) < 0.7:      # restrict to 40-70% GC content
                                yield seq_str + "NGG"


def get_gc(seq_str):
    """ Return GC content from a sequence string. """
    return float(seq_str.count('G') + seq_str.count('C')) / len(seq_str)


def get_targets_bowtie2(inf, hg38):
    """ Run bowtie2 to align 'inf' fasta file to hg38, result in 'outfile' sam file """
    sp.run(['bowtie2', '-k', '1000', '-f', '-x', hg38[:-3], '-U', inf + ".fa", '-S', inf + ".sam"])


def gen_putative(samfile):
    """ Generator of bowtie2 alignments for each putative protospacer sequence of format:
        - boolean: whether it is a new putative protospacer sequence
        - string: sequence of putative protospacer sequence
        - string: '+' or '-' that correspond to the sense of current protospacer sequence alignment
        - string: chromosome string (e.g. 'chr1') of current protospacer sequence alignment
        - int: coordinate of current protospacer sequence alignment
        - int: total number of alignments for current protospacer sequence

    :param samfile: path to SAM file that contains bowtie2 alignment of each protospacer sequence
    """
    with open(samfile, 'r') as sam_it:
        cter = 0
        seq_list = []
        for read in sam_it:     # read every line of SAM file
            # counter to track progress
            if cter % 10000 == 0:
                print("gen_putative(): Read %i lines of SAM file." % cter)
            cter += 1
            # skip SAM header lines
            if read[0] == '@':
                continue
            # read each bowtie2 alignment in SAM format
            row = read.strip().split('\t')
            if row[1] == '0' or row[1] == '16':         # first alignment for new sequence
                sense = '+' if row[1] == '0' else '-'
                seq_list_len = len(seq_list)
                if seq_list_len > 0:
                    for seq in seq_list:
                        yield seq[0], seq[1], seq[2], seq[3], seq[4], seq_list_len
                seq_list = [(True, row[0], sense, row[2], int(row[3]))]
            elif row[1] == '256' or row[1] == '272':    # more alignments for current sequence
                sense = '+' if row[1] == '256' else '-'
                seq_list.append((False, row[0], sense, row[2], int(row[3])))
            elif row[1] == '4':
                pass
            else:
                raise ValueError("gen_putative(): unexpected flag value of %s." % row[1])
        seq_list_len = len(seq_list)
        if seq_list_len > 0:
            for seq in seq_list:
                yield seq[0], seq[1], seq[2], seq[3], seq[4], seq_list_len


def get_targets_stats(generator, outfile):
    """ Given a bowtie2 samfile output containing multiple sequence alignments for each Cas9 target
        sequence, outputs:
    - (*_counts.csv) For each protospacer sequence, tally number of alignments and epigenetic state
      of each alignment position
    - (*_align.csv) For each protospacer sequence, list each alignment with the following columns:
        protospacer + PAM | sequence # | alignment # for specific sequence |
        total # of alignments for specific sequence | chr | coord | gene name | gene orientation |
        protospacer orientation

    :param generator: bowtie2 alignments for each putative protospacer sequence of format:
        - boolean: whether it is a new putative protospacer sequence
        - string: sequence of putative protospacer sequence
        - string: '+' or '-' that correspond to the sense of current protospacer sequence alignment
        - string: chromosome string (e.g. 'chr1') of current protospacer sequence alignment
        - int: coordinate of current protospacer sequence alignment
        - int: total number of alignments for current protospacer sequence
    :param outfile: string path to output file (extension omitted)

    """
    cnt_list, cnt_set, sam_list,  = [], set(), []
    p_read, p_coun, p_gene, p_active, p_chmm = "", 0, 0, 0, [0] * len(CHROMHMM)
    new_ct = 0
    for new_i, seq_i, sen_i, chr_i, cor_i, tnt_i in generator:
        if new_i:                                               # new sequence
            # record previous sequence alignment counts
            if p_coun != 0:
                if p_read in cnt_set:
                    raise ValueError("get_targets_stats(): duplicate and separate alignments")
                cnt_list.append([p_read, p_coun, p_gene, p_active] + p_chmm)
                cnt_set.add(p_read)
            # reset alignment counts for new sequence
            new_ct += 1
            ig = c.is_gene_refseq(chr_i, cor_i)                 # get gene status
            anno_i = c.get_chromhmm_annotation(chr_i, cor_i)    # get ChromHMM annotation
            active_i = c.get_chromhmm_active(anno_i)            # check if epigenetically active
            p_read, p_coun, = seq_i, 1
            p_gene = 1 if ig else 0
            p_active = 1 if active_i == 'active' else 0
            p_chmm = [0] * len(CHROMHMM)
            if anno_i:
                p_chmm[CHROMHMM.index(anno_i)] += 1
            # record new alignment for new sequence
            nam_g, sen_g = (ig[0], ig[1]) if ig else ("", "")
            sam_list.append([seq_i, new_ct, p_coun, tnt_i, chr_i, cor_i, sen_i, nam_g, sen_g,
                             anno_i, active_i])
        else:                                                   # old sequence
            # update alignment counts for current sequence
            ig = c.is_gene_refseq(chr_i, cor_i)                 # get gene status
            anno_i = c.get_chromhmm_annotation(chr_i, cor_i)    # get ChromHMM annotation
            active_i = c.get_chromhmm_active(anno_i)            # check if epigenetically active
            p_coun += 1
            p_gene = p_gene + 1 if ig else p_gene
            p_active = p_active + 1 if active_i == 'active' else p_active
            if anno_i:
                p_chmm[CHROMHMM.index(anno_i)] += 1
            # record new alignment for current sequence
            nam_g, sen_g = (ig[0], ig[1]) if ig else ("", "")
            sam_list.append([seq_i, new_ct, p_coun, tnt_i, chr_i, cor_i, sen_i, nam_g, sen_g,
                             anno_i, active_i])
    # record last sequence alignment counts
    if p_coun != 0:
        if p_read in cnt_set:
            raise ValueError("get_targets_stats(): duplicate and separate alignments")
        cnt_list.append([p_read, p_coun, p_gene, p_active] + p_chmm)
        cnt_set.add(p_read)
    # finalize and save to file
    cnt_np = np.asarray(cnt_list)
    cnt_np = cnt_np[cnt_np[:, 1].astype(int).argsort()[::-1], :]   # sort descending by counts
    cnt_h = ",".join(["protospacer", "# of alignments", "# in genes", "active region"] + CHROMHMM)
    sam_h = ",".join(["protospacer", "sequence #", "alignment #",  "total alignments", "chr",
                      "position", "chr orientation", "gene", "gene orientation",
                      "ChromHMM annotation", "ChromHMM active"])
    np.savetxt(outfile + '_count.csv', cnt_np, fmt='%s', delimiter=',', header=cnt_h)
    np.savetxt(outfile + '_align.csv', np.asarray(sam_list), fmt='%s', delimiter=',', header=sam_h)


def get_targets_dist(alignfile, fileout):
    """ Given get_targets_stats() *_align.csv output, outputs:
    - average distance between adjacent putative on-target sites (*_dist.csv)

    :param alignfile: *_align.csv output from get_targets_stats()
    :param fileout: output file name/path (no extension)

    """
    aln_np = m.load_nparray(alignfile)
    curseq, curind, dist_list = None, [], []
    for i in range(aln_np.shape[0]):       # iterate over alignlist
        row = aln_np[i, :]
        if curseq != row[0]:               # a new sequence
            if curseq is not None:         # and not the very first one in list
                avdist, numalign = _get_targets_dist_helper(aln_np[curind, :])
                dist_list.append([curseq, numalign, avdist])   # calculate avg dist for previous seq
            curind, curseq = [i], row[0]   # for the new sequence, record its indices from alignlist
        else:
            curind.append(i)               # if same sequence as previous row, add its index
    avdist, numalign = _get_targets_dist_helper(aln_np[curind, :])
    dist_list.append([curseq, numalign, avdist])               # calculate avg dist for the last seq
    distnp = np.asarray(dist_list)
    distnp = distnp[distnp[:, 1].astype(int).argsort()[::-1], :]   # sort descending by counts
    np.savetxt(fileout + '_dist.csv', distnp, fmt='%s', delimiter=',')


def _get_targets_dist_helper(aln):
    """ Helper function to calculate average distance between adjacent putative on-target sites

    :param aln: a subset of alignfile (*_align.csv) rows that correspond to all putative on-target
    sites for a particular target sequence

    aln[4] = chr
    aln[5] = position
    """
    aln_sorted = aln[np.lexsort((aln[:, 5].astype(int), aln[:, 4]))]
    numalign, avdist = aln_sorted.shape[0], []
    if numalign > 1:
        for i in range(aln_sorted.shape[0]-1):
            aln1, aln2 = aln_sorted[i, :], aln_sorted[i+1, :]
            if aln1[4] == aln2[4]:
                avdist.append(int(aln2[5]) - int(aln1[5]))
    if len(avdist) > 0:
        return statistics.mean(avdist), numalign
    else:
        return None, numalign


def get_artifical_pe_reads(generator, outfile, hg38savepath, ct_min=100, ct_max=300):
    """
    :param generator: bowtie2 alignments for each putative protospacer sequence of format:
        - boolean: whether it is a new putative protospacer sequence
        - string: sequence of putative protospacer sequence
        - string: '+' or '-' that correspond to the sense of current protospacer sequence alignment
        - string: chromosome string (e.g. 'chr1') of current protospacer sequence alignment
        - int: coordinate of current protospacer sequence alignment
        - int: total number of alignments for current protospacer sequence
    :param outfile: string path to output file (extension omitted)
    :param hg38savepath: save path to hg38 sequence that will be downloaded from NCBI if not already
    :param ct_min: minimum number of target sites for current protospacer sequence
    :param ct_max: maximum number of target sites for current protospacer sequence
    """
    proto_i, align_i = -1, -1
    hg38_initialized(hg38savepath)
    cter = 0
    with open(outfile + "_1.fa", 'w') as f1, open(outfile + "_2.fa", 'w') as f2:
        for new_i, seq_i, sen_i, chr_i, coo_i, tct_i in generator:
            # only get mock PE reads for putative protospacers with 100-300 expected target sites
            if ct_min <= tct_i <= ct_max and chr_i in CHR:
                proto_i = proto_i + 1 if new_i else proto_i
                align_i = 0 if new_i else align_i + 1
                for i, (r1, r2) in enumerate(_get_artificial_pe_reads_helper(chr_i, coo_i)):
                    f1.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                              align_i, i, r1.seq))
                    f2.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
                                                              align_i, i, r2.seq))
                if cter % 10000 == 0:
                    print("get_artifical_pe_reads(): processed %i samples" % cter)
                cter += 1


def _get_artificial_pe_reads_helper(chrom, coord, width_min=200, width_max=600, count=100, rlen=36):
    """ Given chromosome and coordinate position, generate artificial paired-end ChIP-seq reads that
        mimics those from MRE11 ChIP-seq. Those include PE reads that abut the cut site on either
        side, and reads that span the cut site.

    :param chrom: chromosome of expected cut site
    :param coord: coordinate of expected cut site
    :param width_min: minimum width of artificial DNA fragment
    :param width_max: maximum width of artificial DNA fragment
    :param count: number of PE reads to generate
    :param rlen: number of bases to read from the ends of each PE read (DNA fragment)

    :return list of tuples of format (read1, read2). read1 is taken from the sense strand, while
            read2 is the reverse complement of the sense strand (--fr orientation for bowtie2
            that is consistent with Illumina sequencing).
    """
    global hg38id, hg38size, hg38seq
    chr_i = hg38seq[chrom]
    chr_i_len = len(chr_i)
    outlist = []
    rand_steps = np.random.exponential(4, count).astype(int)
    rand_width = np.random.uniform(width_min, width_max, count).astype(int)
    rand_devia = np.random.uniform(-int((width_max-width_min) * 3 / 4), 0, count).astype(int)
    rand_rtype = np.random.randint(0, 3, count)      # left of cut, right of cut, or span cut
    for i in range(count):
        steps, width, devia, rtype = rand_steps[i], rand_width[i], rand_devia[i], rand_rtype[i]
        if rtype == 0:                   # left of cut
            pe_rt = coord + 1 - steps
            pe_lt = pe_rt - width
        elif rtype == 1:                 # right of cut
            pe_lt = coord + steps
            pe_rt = pe_lt + width
        else:                           # span cut
            pe_lt = coord + 1 + devia
            pe_rt = pe_lt + width
        if pe_lt >= 0 and pe_rt < chr_i_len:
            outlist.append((chr_i[pe_lt:pe_lt + rlen + 1],
                            chr_i[pe_rt - rlen:pe_rt + 1].reverse_complement()))
    return outlist


def bowtie2_msa_single(curfile, hg38, k_count=10):
    """ Run bowtie2 to align SINGLE-end reads in '-k' mode, followed by samtools to sort and index
        the alignment.

    :param curfile: path to base file name/path.
                    The SE read files are of format: curfile + '.fa'
                    The output initial SAM file is of format: curfile + '_msa.sam'
                    The output (unsorted) BAM file is of format: curfile + '_msa.bam'
                    The output sorted BAM file is of format: curfile + '_msa_sorted.bam'
                    The output BAM index file is of format: curfile + '_msa_sorted.bam.bai'
    :param hg38: path to hg38 for bowtie2 to use
    :param k_count: number of distinct, valid alignments for each PE read
    """
    sp.run(['bowtie2', '-f', '-p', '12', '--local', '-k', str(k_count), '--no-mixed',
            '--no-discordant', '-x', hg38[:-3], '-U', curfile + '.fa', '-S', curfile + '_msa.sam'])
    sp.run(['samtools', 'view', '-h', '-S', '-b', '-o', curfile + '_msa.bam', curfile + '_msa.sam'])
    sp.run(['samtools', 'sort', '-o', curfile + '_msa_sorted.bam', curfile + '_msa.bam'])
    sp.run(['samtools', 'index', curfile + '_msa_sorted.bam'])


def bowtie2_msa_paired(curfile, hg38, k_count=10):
    """ Run bowtie2 to align paired-end reads in '-k' mode, followed by samtools to sort and index
        the alignment.

    :param curfile: path to base file name/path.
                    The PE read files are of format: curfile + '_1.fa' | curfile + '_1.fa'.
                    The output initial SAM file is of format: curfile + '_msa.sam'
                    The output (unsorted) BAM file is of format: curfile + '_msa.bam'
                    The output sorted BAM file is of format: curfile + '_msa_sorted.bam'
                    The output BAM index file is of format: curfile + '_msa_sorted.bam.bai'
    :param hg38: path to hg38 for bowtie2 to use
    :param k_count: number of distinct, valid alignments for each PE read
    """
    sp.run(['bowtie2', '-f', '-p', '12', '--local', '-k', str(k_count), '-X', '1000', '--no-mixed',
            '--no-discordant', '-x', hg38[:-3], '-1', curfile + '_1.fa', '-2', curfile + '_2.fa',
            '-S', curfile + '_msa.sam'])
    sp.run(['samtools', 'view', '-h', '-S', '-b', '-o', curfile + '_msa.bam', curfile + '_msa.sam'])
    sp.run(['samtools', 'sort', '-o', curfile + '_msa_sorted.bam', curfile + '_msa.bam'])
    sp.run(['samtools', 'index', curfile + '_msa_sorted.bam'])


def parse_msa_sam_single(outfile):
    """ Parse the SAM file output from bowtie2 alignment of paired-end ChIP-seq reads to record the
        information about each PE read, number of alignments, uniqueness, and alignment score.
    :param outfile: path to output csv file (extension omitted)

    Outputs csv file with following columns:
     0. original protospacer sequence
     1-2. chromosome and coordinate of cut site
     3. total number of cut sites for specific protospacer sequence
     4. index of protospacer sequence
     5. index of cut site for a particular protospacer sequence
     6. index of PE ChIP-seq read for a particular cut site
     7. '1' if primary alignment of PE read is to the expected cut site, and there is only one
            alignment or multiple alignments but the primary alignment has the highest score
    """
    global fullflags_single, initflags_single, primflags_single
    """ Parse bowtie2 samfile output, determine all primary/secondary alignments """
    csv_out = open(outfile + '.csv', 'w')
    outrow, outstr, scores = None, "", []
    sam_it = open(outfile + '.sam', 'r')
    for read in sam_it:         # read every alignment of SAM file
        if read[0] == '@':      # skip SAM header lines
            continue
        stp = read.strip()
        row = stp.split('\t')           # read each bowtie2 alignment in SAM format
        n_str = row[0].split('_')       # split the information of each aligned PE read
        seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i = n_str[:7]
        # properly mapped alignment of a new SE-read
        if row[1] in initflags_single:
            # get quality scores and location of alignment
            sco_c = int(re.findall(r'AS:i:[-0-9]+', stp)[0].split(':')[2])  # quality score
            chr_c, coo_c = row[2], int(row[3])              # chromosome and position of alignment
            str_c = "%s_%i_%0.3f" % (chr_c, coo_c, sco_c)   # summary string of alignment
            # primary alignment - the alignment of a new SE-read designated 'primary'
            if row[1] in primflags_single:
                if outrow is not None:  # save previous alignment set
                    len_scores = len(scores)
                    outrow[-1] = '1' if outrow[-1] and (
                            len_scores == 1 or scores[0] > max(scores[1:])) else '0'
                    outrow.append(str(len_scores))
                    outrow.append(outstr)
                    csv_out.write(','.join(outrow) + '\n')
                # determine if primary alignment matches previously assigned location
                int_c = int(coo_i)
                boo = True if chr_c == chr_i and int_c - 2e3 <= coo_c <= int_c + 2e3 else False
                outrow = [seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i, boo]
                outstr, scores = "", []
            outstr += str_c + "|"
            scores.append(sco_c)
        elif row[1] not in fullflags_single:
            raise ValueError("parse_msa_sam_single(): unexpected SAM flag: %s" % row[1])
    if outrow is not None:                  # save final alignment set
        len_scores = len(scores)
        outrow[-1] = '1' if outrow[-1] and (len_scores == 1 or scores[0] > max(scores[1:])) else '0'
        outrow.append(str(len_scores))
        outrow.append(outstr)
        csv_out.write(','.join(outrow) + '\n')
    sam_it.close()
    csv_out.close()


def parse_msa_sam_paired(outfile):
    """ Parse the SAM file output from bowtie2 alignment of paired-end ChIP-seq reads to record the
        information about each PE read, number of alignments, uniqueness, and alignment score.
    :param outfile: path to output csv file (extension omitted)

    Outputs csv file with following columns:
     0. original protospacer sequence
     1-2. chromosome and coordinate of cut site
     3. total number of cut sites for specific protospacer sequence
     4. index of protospacer sequence
     5. index of cut site for a particular protospacer sequence
     6. index of PE ChIP-seq read for a particular cut site
     7. '1' if primary alignment of PE read is to the expected cut site, and there is only one
            alignment or multiple alignments but the primary alignment has the highest score
    """
    global fullflags_paired, initflags_paired, primflags_paired
    """ Parse bowtie2 samfile output, determine all primary/secondary alignments """
    csv_out = open(outfile + '.csv', 'w')
    outrow, outstr, scores = None, "", []
    sam_it = open(outfile + '.sam', 'r')
    for read in sam_it:         # read every alignment of SAM file
        if read[0] == '@':      # skip SAM header lines
            continue
        stp = read.strip()
        row = stp.split('\t')           # read each bowtie2 alignment in SAM format
        n_str = row[0].split('_')       # split the information of each aligned PE read
        seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i = n_str[:7]
        # first in pair, properly mapped - one (first) alignment of a new PE-read.
        if row[1] in initflags_paired:
            # get quality scores and location of alignment
            AS = int(re.findall(r'AS:i:[-0-9]+', stp)[0].split(':')[2])
            YS = int(re.findall(r'YS:i:[-0-9]+', stp)[0].split(':')[2])
            sco_c = (AS + YS) / 2.0                         # quality score of alignment
            chr_c, coo_c = row[2], int(row[3])              # chromosome and position of alignment
            str_c = "%s_%i_%0.3f" % (chr_c, coo_c, sco_c)   # summary string of alignment
            # first in pair, primary alignment - the alignment of a new PE-read designated 'primary'
            if row[1] in primflags_paired:
                if outrow is not None:  # save previous alignment set
                    len_scores = len(scores)
                    outrow[-1] = '1' if outrow[-1] and (
                            len_scores == 1 or scores[0] > max(scores[1:])) else '0'
                    outrow.append(str(len_scores))
                    outrow.append(outstr)
                    csv_out.write(','.join(outrow) + '\n')
                # determine if primary alignment matches previously assigned location
                int_c = int(coo_i)
                boo = True if chr_c == chr_i and int_c - 2e3 <= coo_c <= int_c + 2e3 else False
                outrow = [seq_i, chr_i, coo_i, tnt_i, proto_i, align_i, pread_i, boo]
                outstr, scores = "", []
            outstr += str_c + "|"
            scores.append(sco_c)
        elif row[1] not in fullflags_paired:
            raise ValueError("parse_msa_sam_paired(): unexpected SAM flag: %s" % row[1])
    if outrow is not None:                  # save final alignment set
        len_scores = len(scores)
        outrow[-1] = '1' if outrow[-1] and (len_scores == 1 or scores[0] > max(scores[1:])) else '0'
        outrow.append(str(len_scores))
        outrow.append(outstr)
        csv_out.write(','.join(outrow) + '\n')
    sam_it.close()
    csv_out.close()


def get_msa_stats(outfile):
    """ Given the output of 'parse_msa_sam_paired()' or 'parse_msa_sam_single()', determine
        the number of reads that:
        (1) have more than 1 equivalently good alignment [bad outcome]
        (2) have more than 1 alignment, but one unique 'best' alignment, i.e. highest alignment
            score [good outcome]
        (3) only have 1 alignment [best outcome]
    :param outfile: path of output of 'parse_msa_sam_*()', extension omitted

    Outputs csv file with following columns:
     0. original protospacer sequence
     1. total number of cut sites for specific protospacer sequence
     2. # of PE reads with >1 equivalently good alignment [bad outcome]
     3. # of PE reads with >1 alignment, but one unique 'best' alignment [good outcome]
     4. total number of alignments for each protospacer sequence
    """
    csv_out = open(outfile + '_stats.csv', 'w')
    msa = m.load_nparray(outfile + ".csv")
    msa_sorted = msa[np.lexsort((msa[:, 6].astype(int), msa[:, 5].astype(int),
                                 msa[:, 4].astype(int), -msa[:, 3].astype(int)))]
    proto_prev = None
    sumrow, numct = [], []
    for msa_i in msa_sorted:
        seq_i, chr_i, coo_i, tnt_i = msa_i[:4]
        proto_i, boo_i, num_i = int(msa_i[4]), msa_i[7], int(msa_i[8])
        if proto_i != proto_prev:
            if len(numct) > 0:
                sumrow[2], sumrow[3], sumrow[4] = str(sumrow[2]), str(sumrow[3]), str(sumrow[4])
                csv_out.write(','.join(sumrow) + '\n')
            sumrow = [seq_i, tnt_i, 0, 0, 0]
            numct = []
            proto_prev = proto_i
        sumrow[2] = sumrow[2] + 1 if boo_i == '0' else sumrow[2]
        sumrow[3] = sumrow[3] + 1 if boo_i == '1' and num_i > 1 else sumrow[3]
        sumrow[4] += 1
        numct.append(num_i)
    if len(numct) > 0:
        sumrow[2], sumrow[3], sumrow[4] = str(sumrow[2]), str(sumrow[3]), str(sumrow[4])
        csv_out.write(','.join(sumrow) + '\n')
    csv_out.close()


def target_gen(alignfile, span_r, guide):
    """ Generator to yield all putative on-target sites for a given protospacer

    :param alignfile: CSV file generated from get_targets_stats() output
    :param span_r: radius of window from peak center for analysis of associated epigenetic info
    :param guide: on-target protospacer sequence (no PAM)
    :yield: ( span_rs, cut_i, sen_i, pam_i, gui_i, mis_i, guide )
            ( region string, cut site, sense/antisense, PAM, discovered protospacer,
              # mismatches, non-mismatched protospacer )
            - The format is designed to be consistent with other generators that may yield
              mismatched sequences. Note that no mismatched sequences will be outputted by design.

    """
    aln = m.load_nparray(alignfile)
    hg38dict = c.hg_dict('hg38')
    pam_i = 'NGG'
    outlist = []
    for i in range(aln.shape[0]):
        row = aln[i, :]
        if row[0] == guide + pam_i:
            chr_i = row[4]
            sen_i = row[6]
            if sen_i == '+':
                cut_i = int(row[5]) + 16
            else:
                cut_i = int(row[5]) + 6
            span_sta = max(1, cut_i - span_r)
            span_end = min(hg38dict[chr_i], cut_i + span_r)
            span_rs = "%s:%i-%i" % (chr_i, span_sta, span_end)
            mis_i = 0
            outlist.append([chr_i, span_rs, cut_i, sen_i, pam_i, guide, mis_i, guide])
    outlist = sorted(outlist, key=lambda x: (x[2]))
    outlist = sorted(outlist, key=lambda x: (x[0]))
    for out in outlist:
        yield tuple(out[1:])


def get_bamfile_pe_reads(generator, bamfile, outfile):
    """

    """
    bamin = pysam.AlignmentFile(bamfile, 'rb')
    list_f1, list_f2 = [], []
    num_sites = 0
    for i, (rs, cut, sen, pam, gui, mis, tar) in enumerate(generator):  # iterate over each cut site
        num_sites += 1
        # get window centered at Cas9 cut site
        chr_i = re.split('[:-]', rs)[0]
        for j, (read1, read2) in enumerate(c.read_pair_generator(bamin, rs)):
            s1, q1, s2, q2 = _get_bamfile_pe_reads_helper(read1, read2)
            list_f1.append((gui + "NGG", chr_i, cut, -1, 1, i, j, s1))
            list_f2.append((gui + "NGG", chr_i, cut, -1, 1, i, j, s2))
    # write all reads in target site window to paired fastq files in original '-fr' format
    with open(outfile + '_1.fa', 'w') as f1, open(outfile + '_2.fa', 'w') as f2:
        for l1, l2 in zip(list_f1, list_f2):
            f1.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (
                l1[0], l1[1], l1[2], num_sites, l1[4], l1[5], l1[6], l1[7]))
            f2.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (
                l2[0], l2[1], l2[2], num_sites, l2[4], l2[5], l2[6], l2[7]))


def _get_bamfile_pe_reads_helper(read1, read2):
    """ Extract read pair sequence and quality from original paired-end ChIP-seq data

    :param read1: read #1 of pair in pysam AlignedSegment format
    :param read2: read #2 of pair in pysam AlignedSegment format
    :return read1 sequence, read1 quality, read2 sequence, read2 quality from original ChIP-seq data

    """
    if read1.is_reverse and not read2.is_reverse:
        return c.get_reverse_complement(read1.seq), read1.qual[::-1], read2.seq, read2.qual
    elif read2.is_reverse and not read1.is_reverse:
        return read1.seq, read1.qual, c.get_reverse_complement(read2.seq), read2.qual[::-1]
    else:
        raise ValueError("_get_bamfile_pe_reads_helper(): Unexpected paired-end read conditions.")


# def site_uniqueness(chromosome, coordinate):
#
#     for i, (r1, r2) in enumerate(_get_artificial_pe_reads_helper(chromosome, coordinate)):
#         f1.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
#                                                   align_i, i, r1.seq))
#         f2.write(">%s_%s_%i_%i_%i_%i_%i\n%s\n" % (seq_i, chr_i, coo_i, tct_i, proto_i,
#                                                   align_i, i, r2.seq))
