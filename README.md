###################################################
#       RNAseq analysis pipeline                  #
#       Michelle Percharde, PhD 2016              #
#        michelle.percharde@ucsf.edu              #
#                v1.2                             #
###################################################


Script that takes input with name "sample.fq" or "sample.fq.gz" and generates sorted bams
Bam files can then be DLed and fed into FeatureCounts in R

UPDATE 12/18/16: Tophat working with .gz or .fq files.
                 Currently, alignment is set to mm10 + ERCCs. 


FOLDERS NEEDED IN ROOT:

raw/ (where raw files are)

NB: test-fq is a folder of two small, test .fq files that can be used as input folder to check pipeline.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ### Usage: ./runRNAseq2.sh [-h] [-g] [-c] [-i <path/to/files/>] ###

    -h        Help mode, prints Usage
    -g        Input FASTQ is .gz compressed
    -c        Libraries are CLONTECH
    -i        input file directory/. use ./ for current dir (not recommended)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PIPELINE STEPS

1. trim_galore, fastQC analysis on trimmed files in trimmed/fastqc

2. tophat2, currently set to align to mm10 (+ERCCs)

3. samtools sort

4. exit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KEY NOTES:

1. Create a sep. dir where raw files are kept. Raw files should end .fq or .gz NOT .fastq

2. If your file is compressed as .fq.gz, use -g option for correct file nomenclature

3. If your libraries are clontech, use -c to trim off first 3 repetitive bp

   Good luck!
