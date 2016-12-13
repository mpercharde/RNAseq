#!/bin/bash

###################################################
#       RNAseq analysis pipeline                  #
#       Michelle Percharde, PhD 2016              #
#                                                 #
#                v1.1                             #
###################################################

#~~~~~~~~~~EDIT BELOW THIS LINE ~~~~~~~~~~~~~~~~~~#

#README~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Script that takes input with name "sample.fq" or "sample.fq.gz" and generates sorted bams
# Bam files can then be DLed and fed into FeatureCounts in R

#usage: ./runRNAseq.sh [options] [-i path/to/folder/]

#FOLDERS NEEDED IN ROOT:
  #raw/ (where raw files are)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#things to add - print out help messages
#              - error check folders created

flagcheck=0

while getopts ':gchi:' flag; do
    case ${flag} in
      i) dir=$OPTARG
        flagcheck=1 ;;
        g) gz='true' ;;
        c) clontech='true' ;;
        h) echo ""
           echo "Usage: $0 [-h] [-g] [-c] [-i <path/to/files/>]"
           echo ""
  	       echo "    -h        Help mode, prints Usage"
  	       echo "    -g        Input FASTQ is .gz compressed"
           echo "    -c        Libraries are CLONTECH"
  	       echo "    -i        input file directory/. use "./" for current dir (not rec)"
  	       echo ""
           flagcheck=1
           exit 1 ;;
        \?)
          echo ""
          echo "Invalid option, type -h for help"
          echo ""
          flagcheck=1
          exit 1 ;;
    esac
done

if [ "$flagcheck" == 0 ] ; then
  echo ""
  echo "Incorrect or no -i specified, please read help -h and try again!"
  echo ""
  exit 1
fi

if [ "$gz" == "true" ]; then
  echo ""
  echo "your file is compressed, trim and align will be run on .gz files"
fi

mkdir -p trimmed/fastqc/
mkdir -p sorted_bam/

for file in "$dir"* ; do
    echo ""
    if [ "$gz" == "true" ]; then
      name=$(basename $file .fq.gz)
    else
      name=$(basename $file .fq)
    fi
    echo "analysing file: $name"

    if [ "$clontech" == "true" ]; then
      echo ""
      echo "analyzing clontech libraries, adaptor trimming will include clipping 3bp"
      trim=$(echo "$name is clontech")
      echo ""
    else
      echo ""
      trim=$(echo "$name is nebnext")
      echo ""
    fi

    # if [ "$clontech" == "true" ]; then
      # trim = "trim_galore --clip_R 3" --fastqc --fastqc_args " --outdir trimmed/fastqc/" --illumina $file -o trimmed/
      # trim = "trim_galore" trim_galore --clip_R 3 --fastqc --fastqc_args " --outdir trimmed/fastqc/" --illumina $file -o trimmed/
    # else
    #   trim = "trim_galore --fastqc --fastqc_args " --outdir trimmed/fastqc/" --illumina $file -o trimmed/"
    # fi

    echo "1. trimming $name"
    echo ""
    echo $trim
    echo ""
    # trim_galore --fastqc --fastqc_args " --outdir trimmed/fastqc/" --illumina $file -o trimmed/ #check what to happen if not gz.

    echo "2. aligning $name to mm10 plus ERCCs"
    echo ""
    mkdir ${name}_aligned/
  # tophat -o ${name}_aligned/ -p 4 -g 20 -G /data/refs/mm10/BTmm10_ercc/genes_ercc.gtf --no-coverage-search --library-type fr-firststrand \
  # --no-novel-indels /data/refs/mm10/BTmm10_ercc/mm10_ercc trimmed/${name}_trimmed.fq

    echo "3. sorting $name bam file ready for DL"
    echo ""
    # samtools sort -o sorted_bam/${name}.sorted.bam ${name}_aligned/accepted_hits.bam
    echo "$name DONE!"
    echo ""
done
