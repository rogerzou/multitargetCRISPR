#!/usr/bin/env bash

########### USER ENTRY SECTION ###########
# Enter paths to all samples for processing into this array
download_path=/mnt/c/Users/rzou4/Downloads/
# List of all SRA samples with their labels. Comment out samples that are not needed.
declare -a filelist=(\
"SRR502470" "HEK_POLR2A_1"
"SRR442119" "HEK_POLR2A_2"
)
##########################################

# processing arguments, proceed with bioinformatics pipeline
main() {
  downloadfastq "${filelist[@]}"
}

# pipeline for downloading FASTQ reads from SRA
downloadfastq() {

  cd $download_path

  arraylength=${#filelist[@]}
  for (( i=0; i<${arraylength}; i+=2 ));
  do
    var1="${filelist[$i]}"
    var2="${filelist[$i+1]}"
    cd $download_path
    prefetch $var1 -O $var1
    cd $var1/$var1
    fasterq-dump "${var1}.sra"
  done

}

main "$@"; exit
