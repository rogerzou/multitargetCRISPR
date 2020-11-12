

filelist = [
"/mnt/d/201022_chipseq/A36",
"/mnt/d/201022_chipseq/A37",
"/mnt/d/201022_chipseq/A38",
"/mnt/d/201022_chipseq/A39",
"/mnt/d/201022_chipseq/A40",
"/mnt/d/201022_chipseq/A41",
"/mnt/d/201022_chipseq/A42",
"/mnt/d/201022_chipseq/A43",
"/mnt/d/201022_chipseq/A44",
"/mnt/d/201022_chipseq/A45",
"/mnt/d/201022_chipseq/A46"
]
matchstr = "CGCCATCACGCCT"


def process_raw(file_r1):
    print("process_bliss.py: Processing %s" % file_r1)
    total_cter, dup_cter = 0, 0
    umi_set = set()
    with open(file_r1 + "_mod_1.fastq", 'w') as w1, open(file_r1 + "_mod_2.fastq", 'w') as w2, \
            open(file_r1 + "_1.fastq", 'r') as r1, open(file_r1 + "_2.fastq", 'r') as r2:
        while True:
            r1_l1, r2_l1 = r1.readline().replace('\n', ''), r2.readline().replace('\n', '')
            r1_l2, r2_l2 = r1.readline().replace('\n', ''), r2.readline().replace('\n', '')
            r1_l3, r2_l3 = r1.readline().replace('\n', ''), r2.readline().replace('\n', '')
            r1_l4, r2_l4 = r1.readline().replace('\n', ''), r2.readline().replace('\n', '')
            if not r1_l1:
                break
            matchind = r1_l2.find(matchstr)
            umi_i = r1_l2[:matchind]
            if matchind == 12:
                total_cter += 1
                if umi_i not in umi_set:
                    umi_set.add(umi_i)
                    r1_l1, r2_l1 = r1_l1 + "+" + umi_i, r2_l1 + "+" + umi_i
                    r1_l2, r1_l4 = r1_l2[matchind + 13:], r1_l4[matchind + 13:]
                    w1.write("\n".join([r1_l1, r1_l2, r1_l3, r1_l4, ""]))
                    w2.write("\n".join([r2_l1, r2_l2, r2_l3, r2_l4, ""]))
                else:
                    dup_cter += 1
    print("process_bliss.py: %i out of %i reads had duplicate UMIs (PCR duplicates)." %
          (dup_cter, total_cter))


for file in filelist:
    process_raw(file)
