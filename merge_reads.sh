#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
"/mnt/d/210225_atac/N16" \
"/mnt/d/210225_atac/N18" \
"/mnt/d/210225_atac/N19" \
"/mnt/d/210225_atac/N01" \
"/mnt/d/210225_atac/N08" \
"/mnt/d/210225_atac/N09" \
"/mnt/d/210225_atac/N10" \
"/mnt/d/210225_atac/N04" \
"/mnt/d/210225_atac/N05" \
"/mnt/d/210225_atac/N06" \
)

##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  genome="hg38"
  while getopts 'p:g:' opt; do       # pass number of threads with -p | which genome with -g
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      g) genome="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads] [-g genome for alignment (hg38, hg19, mm10)]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  echo "Genome file previously for alignment: $genome"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    mergebam ${file} ${genome} &
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
mergebam() {

  # Remove potential PCR duplicates
  samtools merge \
  "$1_$2_merged.bam" "$1_L1_$2_final.bam" "$1_L2_$2_final.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$1_$2_merged.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$1_$2_merged.bam" > "$1_$2_flagstats.txt"

}

main "$@"; exit
