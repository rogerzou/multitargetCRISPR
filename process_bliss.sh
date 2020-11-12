#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
"/mnt/d/201022_chipseq/A36_mod" \
"/mnt/d/201022_chipseq/A37_mod" \
"/mnt/d/201022_chipseq/A38_mod" \
"/mnt/d/201022_chipseq/A39_mod" \
"/mnt/d/201022_chipseq/A40_mod" \
"/mnt/d/201022_chipseq/A41_mod" \
"/mnt/d/201022_chipseq/A42_mod" \
"/mnt/d/201022_chipseq/A43_mod" \
"/mnt/d/201022_chipseq/A44_mod" \
"/mnt/d/201022_chipseq/A45_mod" \
"/mnt/d/201022_chipseq/A46_mod" \
)
# Enter path to indexed genome
genomepath="/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38"
##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  while getopts 'p:' opt; do
    case "$opt" in
      p) numthreads="$OPTARG"
        ;;
      \?) echo "Usage: $(basename $0) [-p number of threads]"
      exit 1
        ;;
    esac
  done
  echo "Number of parallel processes: $numthreads"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    align2bam ${genomepath} ${file} &
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam() {

  # Align reads to either hg38 or mm10 using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2.sam" ;
  # bowtie2 -p 6 -q --local -x $1 -U "$2_2.fastq" -S "$2.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -q 25 \
  "$2.sam" > "$2_unsorted.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_unsorted.bam" > "$2_sorted.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_sorted.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$2_sorted.bam" > "$2_flagstats.txt"

}

main "$@"; exit
