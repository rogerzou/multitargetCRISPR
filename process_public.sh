#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
#"/mnt/d/public/RNAseq_HEK293_SRR5627161" \
"/mnt/d/public/MNaseseq_HEK293_ERR2403161" \
"/mnt/d/public/ATACseq_HEK293_SRR6418075" \
"/mnt/d/public/POLR2A_HEK293_SRR502470" \
"/mnt/d/public/POLR2A_HEK293_SRR442119" \
)
declare -a typelist=(\
# "pe" \
 "pe" \
 "pe" \
"se" \
"se" \
)
declare -a genolist=(\
# "hg38" \
 "hg38" \
 "hg38" \
"hg38" \
"hg38" \
)

# Enter path to indexed genome
hg38path="/home/roger/bioinformatics/hg38/hg38"
hg19path="/home/roger/bioinformatics/hg19/hg19"
mm10path="/home/roger/bioinformatics/mm10/mm10"
##########################################

# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  genomepath=""
  genome=""
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
  for index in ${!filelist[*]}; do
    ((i=i%numthreads)); ((i++==0)) && wait
    if [ ${genolist[$index]} = "hg19" ]; then
      genomepath="$hg19path"
      genome="hg19"
    elif [ ${genolist[$index]} = "mm10" ]; then
      genomepath="$mm10path"
      genome="mm10"
    else
      genomepath="$hg38path"
      genome="hg38"
    fi
    if [ ${typelist[$index]} = "pe" ]; then
      align2bam_pe ${genomepath} ${filelist[$index]} ${genome} &
    else
      align2bam_se ${genomepath} ${filelist[$index]} ${genome} &
    fi
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam_se() {

  # Align reads using bowtie2
  bowtie2 -p 6 -q --local -x $1 -U "$2.fastq" -S "$2_$3.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25
  samtools view -h -S -b -q 25 \
  "$2_$3.sam" > "$2_$3_unsorted.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_$3_unsorted.bam" > "$2_$3_final.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_$3_final.bam"
}

align2bam_pe() {

  # Align reads using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2_$3.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -F 0x08 -q 25 \
  "$2_$3.sam" > "$2_$3_unsorted.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_$3_unsorted.bam" > "$2_$3_final.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_$3_final.bam"

}

main "$@"; exit
