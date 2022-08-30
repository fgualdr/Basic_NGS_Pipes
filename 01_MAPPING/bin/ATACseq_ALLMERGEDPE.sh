#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ATACseq MERGE ALL CALL
#PBS -l select=1:ncpus=4:mem=50g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V

source activate EnvATAC

macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/*R[0-9].bam --name=ALL_MERGED.nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g mm --call-summits
macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/*R*.rmdup.bam --name=ALL_MERGED.rmdup.nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g mm --call-summits

source deactivate
                

