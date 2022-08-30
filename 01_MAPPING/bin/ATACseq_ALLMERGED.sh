#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ATACseq MERGE ALL CALL
#PBS -l select=1:ncpus=4:mem=12g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V

source activate EnvATAC

macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/*.5Shift.nochrM.bam --name=ALL_MERGED.nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAM -B -g mm --nomodel --extsize 146 --nolambda --keep-dup 'all' --call-summits 

source deactivate
                