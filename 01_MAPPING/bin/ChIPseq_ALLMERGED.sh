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

if [ "$INPUT_CHIP" = "NA" ]; then
    echo NO INPUT
    macs2 callpeak -t $OUT_FOLDER/2_BAM/*.clean.bam  --name=ALL_MERGED.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g mm --extsize 250 --keep-dup 'all' --call-summits
    #macs2 callpeak -t $OUT_FOLDER/2_BAM/*.clean.filterdup.bam  --name=ALL_MERGED.filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ -B -g hs --extsize 250 --call-summits
else
    macs2 callpeak -t $OUT_FOLDER/2_BAM/*.clean.bam  -c $INPUT_CHIP/INPUT_BMDM.bam --name=ALL_MERGED.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g mm --extsize 250 --keep-dup 'all' --call-summits
    #macs2 callpeak -t $OUT_FOLDER/2_BAM/*.clean.filterdup.bam  -c $INPUT_CHIP/INPUT.filterdup.bam --name=ALL_MERGED.filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ -B -g hs --extsize 250 --call-summits
fi                  

source deactivate
                
