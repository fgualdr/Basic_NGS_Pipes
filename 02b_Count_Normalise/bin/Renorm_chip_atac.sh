#!/bin/bash

#PBS -l select=1:ncpus=8:mem=20g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V

echo $BAM_DIR
echo $OUT_DIR
mkdir $OUT_DIR
echo $NORM_DAT # 
echo $col_ID # 

source activate EnvAlignRNAseq

for f in $BAM_DIR/*.bam
do

    filename=$(basename "$f")
    identifier=${filename%%.*}
    col=$(awk -F'\t' -v var="$col_ID" -v id="^$identifier" '{
            for(i=1;i<=NF;i++) {
                if($i == var) {
                    printf(i)
                }
            }
            exit 0
            }' $NORM_DAT)
    scale_factor=$(awk -v var="^$identifier" -v col="$col" '$0~var{print $col}' $NORM_DAT)

    bamCoverage --bam $f -o $OUT_DIR/${identifier}.scaled.bw --binSize 1 --extendReads --scaleFactor $scale_factor --numberOfProcessors 8 --ignoreDuplicates 

done

conda deactivate
