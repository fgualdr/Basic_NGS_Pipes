#!/bin/bash

#PBS -l select=1:ncpus=8:mem=20g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V

##copyrights(c) Francesco Gualdrini

## This tool generate normalised BigWig based on a Normalisation file generated with the 03_RNAseq_Read_counts_normalise.R pipe
## To have the script running you need the  EnvToolkit environment (bedtools and ucsc-wigToBigWig)

echo "Below the variables you provided:"
echo " "
echo $BAM_FILES
echo $NORM_FILE
echo $OUT_FOLDER
echo $ChromSizes_path
echo $STRANDNESS

ls -l $BAM_FILES

source activate EnvAlignRNAseq

cd $OUT_FOLDER
mkdir BIGWIG_NORM
cd BIGWIG_NORM

for f in $BAM_FILES/*.clean.bam
do
    if [ "$STRANDNESS" = "NO" ]; then
    	filename=$(basename "$f" .bam)
    	identifier=$(basename "$f" .clean.bam)
     	echo $f
    	echo $filename
    	echo $identifier
    	scale_factor=$(awk -v var="^$identifier" '$0~var{print $NF}' $NORM_FILE)
       	echo $scale_factor
        bamCoverage -p 8 -bs 1 --scaleFactor $scale_factor -b $f -o $OUT_FOLDER/BIGWIG_NORM/"$identifier".normalised.allreads.bw
    else
        filename=$(basename "$f" .bam)
        identifier=$(basename "$f" .clean.bam)
        echo $f
        echo $filename
        echo $identifier
        scale_factor=$(awk -v var="^$identifier" '$0~var{print $NF}' $NORM_FILE)
        echo $scale_factor
        bamCoverage -p 8 -bs 1 --scaleFactor $scale_factor --filterRNAstrand forward -b $f -o $OUT_FOLDER/BIGWIG_NORM/"$identifier".normalised.allreads.reverse.bw
        bamCoverage -p 8 -bs 1  --scaleFactor $scale_factor --filterRNAstrand reverse -b $f -o $OUT_FOLDER/BIGWIG_NORM/"$identifier".normalised.allreads.forward.bw
	fi
done

cd $OUT_FOLDER
rm -rf BIGWIG

conda deactivate
