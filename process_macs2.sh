#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a dirlist=(\
"/mnt/d/200206_chipseq_AluGG/" \
"/mnt/d/200206_chipseq_AluGG/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200206_chipseq_AluGG/" \
"/mnt/d/200206_chipseq_AluGG/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200316_chipseq/" \
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
"/mnt/d/200804_chipseq/"
)

declare -a samplist=(\
"AluGG-Cas9_hg19_final" \
"AluGG-MRE11_hg19_final" \
"AluCT-cas9-rep1_hg19_final" \
"AluCT-mre11-rep1_hg19_final" \
"AluTA-cas9-rep1_hg19_final" \
"AluTA-mre11-rep1_hg19_final" \
"AluGG-Cas9_hg38_final" \
"AluGG-MRE11_hg38_final" \
"AluCT-cas9-rep1_hg38_final" \
"AluCT-mre11-rep1_hg38_final" \
"AluTA-cas9-rep1_hg38_final" \
"AluTA-mre11-rep1_hg38_final" \
"A03_hg19_final" \
"A06_hg19_final" \
"A09_hg19_final" \
"A12_hg19_final" \
"A03_hg38_final" \
"A06_hg38_final" \
"A09_hg38_final" \
"A12_hg38_final" \
)

declare -a ctrllist=(\
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg19_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A16_cas9_hg38_final" \
"/mnt/d/200212_chipseq_WT1/A17_mre11_hg38_final" \
)
# setup commands
cd /mnt/c/users/Roger/bioinformatics
source /mnt/c/Users/Roger/bioinformatics/python3biovenv/bin/activate
##########################################


# processing arguments, proceed with macs2
main() {
  numthreads=1                      # default number of parallel processes
  while getopts 'p:' opt; do       # pass number of threads with -p
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  for index in ${!samplist[*]}; do
    ((i=i%numthreads)); ((i++==0)) && wait
    macs2_call ${dirlist[$index]} ${samplist[$index]} ${ctrllist[$index]} &
  done
}


macs2_call() {
  macs2 callpeak -t "$1$2.bam" -c "$3.bam" --outdir "$1/macs/" --name $2 -f BAMPE -g hs
}

main "$@"; exit
