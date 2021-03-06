#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
declare -a filelist=(\
# "/mnt/d/200206_chipseq_AluGG/AluGG-53BP1" \
# "/mnt/d/200206_chipseq_AluGG/AluGG-Cas9" \
# "/mnt/d/200206_chipseq_AluGG/AluGG-gH2AX" \
# "/mnt/d/200206_chipseq_AluGG/AluGG-MRE11" \
# "/mnt/d/200316_chipseq/AluCT-53bp1-rep1" \
# "/mnt/d/200316_chipseq/AluCT-cas9-rep1" \
# "/mnt/d/200316_chipseq/AluCT-gh2ax-rep1" \
# "/mnt/d/200316_chipseq/AluCT-mre11-rep1" \
# "/mnt/d/200316_chipseq/AluTA-53bp1-rep1" \
# "/mnt/d/200316_chipseq/AluTA-cas9-rep1" \
# "/mnt/d/200316_chipseq/AluTA-gh2ax-rep1" \
# "/mnt/d/200316_chipseq/AluTA-mre11-rep1" \
# "/mnt/d/200212_chipseq_WT1/A15_53bp1" \
# "/mnt/d/200212_chipseq_WT1/A16_cas9" \
# "/mnt/d/200212_chipseq_WT1/A17_mre11" \
# "/mnt/d/200212_chipseq_WT1/A18_gh2ax" \
# "/mnt/d/200316_chipseq/AluGG-mre11-noD-rep1" \
# "/mnt/d/200316_chipseq/AluGG-mre11-PKi-rep1" \
# "/mnt/d/200316_chipseq/AluGG-53bp1-noD-rep1" \
# "/mnt/d/200316_chipseq/AluGG-53bp1-PKi-rep1" \
# "/mnt/d/200804_chipseq/A01" \
# "/mnt/d/200804_chipseq/A02" \
# "/mnt/d/200804_chipseq/A03" \
# "/mnt/d/200804_chipseq/A04" \
# "/mnt/d/200804_chipseq/A05" \
# "/mnt/d/200804_chipseq/A06" \
# "/mnt/d/200804_chipseq/A07" \
# "/mnt/d/200804_chipseq/A08" \
# "/mnt/d/200804_chipseq/A09" \
# "/mnt/d/200804_chipseq/A10" \
# "/mnt/d/200804_chipseq/A11" \
# "/mnt/d/200804_chipseq/A12" \
# "/mnt/d/201012_chipseq/A01" \
# "/mnt/d/201012_chipseq/A02" \
# "/mnt/d/201012_chipseq/A03" \
# "/mnt/d/201012_chipseq/A04" \
# "/mnt/d/201012_chipseq/A14" \
# "/mnt/d/201012_chipseq/A15" \
# "/mnt/d/201012_chipseq/A16" \
# "/mnt/d/201012_chipseq/A17" \
# "/mnt/d/201012_chipseq/A18" \
# "/mnt/d/201207_atac/N701" \
# "/mnt/d/201207_atac/N703" \
# "/mnt/d/201207_atac/N704" \
# "/mnt/d/201207_atac/N705" \
# "/mnt/d/201207_atac/N706" \
# "/mnt/d/210225_atac/N16_L1" \
# "/mnt/d/210225_atac/N18_L1" \
# "/mnt/d/210225_atac/N19_L1" \
# "/mnt/d/210225_atac/N01_L1" \
# "/mnt/d/210225_atac/N08_L1" \
# "/mnt/d/210225_atac/N09_L1" \
# "/mnt/d/210225_atac/N10_L1" \
# "/mnt/d/210225_atac/N04_L1" \
# "/mnt/d/210225_atac/N05_L1" \
# "/mnt/d/210225_atac/N06_L1" \
# "/mnt/d/210225_atac/N16_L2" \
# "/mnt/d/210225_atac/N18_L2" \
# "/mnt/d/210225_atac/N19_L2" \
# "/mnt/d/210225_atac/N01_L2" \
# "/mnt/d/210225_atac/N08_L2" \
# "/mnt/d/210225_atac/N09_L2" \
# "/mnt/d/210225_atac/N10_L2" \
# "/mnt/d/210225_atac/N04_L2" \
# "/mnt/d/210225_atac/N05_L2" \
# "/mnt/d/210225_atac/N06_L2" \
)

# Enter path to indexed genome
hg38path="/mnt/c/Users/Roger/bioinformatics/hg38_bowtie2/hg38"
hg19path="/mnt/c/Users/Roger/bioinformatics/hg19_bowtie2/hg19"
mm10path="/mnt/c/Users/Roger/bioinformatics/mm10_bowtie2/mm10"
##########################################


# processing arguments, proceed with bioinformatics pipeline
main() {
  numthreads=1                      # default number of parallel processes
  genome=""
  genomepath=""
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
  if [ "$genome" = "hg19" ]; then
    genomepath="$hg19path"
  elif [ "$genome" = "mm10" ]; then
    genomepath="$mm10path"
  else
    genomepath="$hg38path"
    genome="hg38"
  fi
  echo "Number of parallel processes: $numthreads"
  echo "Genome file for alignment: $genome | $genomepath"
  for file in "${filelist[@]}"; do  # process each file
    ((i=i%numthreads)); ((i++==0)) && wait
    align2bam ${genomepath} ${file} ${genome} &
  done
}


# main bioinformatics pipeline (alignment to indexing and read statistics)
align2bam() {

  # Align reads to either hg38 or mm10 using bowtie2
  bowtie2 -p 6 -q --local -X 1000 -x $1 \
  -1 "$2_1.fastq" -2 "$2_2.fastq" -S "$2_$3.sam" ;

  # Convert from SAM to BAM, filter for mapping quality >=25 and singleton reads
  samtools view -h -S -b -F 0x08 -q 25 \
  "$2_$3.sam" > "$2_$3_unsorted.bam" ;

  # Add mate score tags to ensure that paired-end reads contain correct information about the mate reads
  samtools fixmate -m \
  "$2_$3_unsorted.bam" "$2_$3_fixmate.bam" ;

  # Sort BAM file entries by genomic position
  samtools sort \
  "$2_$3_fixmate.bam" > "$2_$3_sorted.bam" ;

  # Remove potential PCR duplicates
  samtools markdup -r \
  "$2_$3_sorted.bam" "$2_$3_final.bam" ;

  # Index the sorted BAM file
  samtools index \
  "$2_$3_final.bam" ;

  # Retrieve read count statistics
  samtools flagstat \
  "$2_$3_final.bam" > "$2_$3_flagstats.txt"

}

main "$@"; exit
