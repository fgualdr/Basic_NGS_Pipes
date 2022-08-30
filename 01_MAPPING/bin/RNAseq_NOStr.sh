#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by RNAseq_NOStrand
# Reads are going to be mapped using tophat and big wig non-strand specific created

#PBS -l select=1:ncpus=4:mem=16g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -l maxarray_5=1

source activate EnvAlignRNAseq

Scale_to=10000000

echo "Below the variables you provided:" 
echo " " 
echo $TYPE_AL 
echo $OUT_FOLDER 
echo $BOWTIE2_INDEX 
echo $BOWTIE_INDEX 
echo $REFSEQ 
echo $BL 
echo $ChromSizes_path 
echo $Path_to_ShortStack 
echo $INPUT_CHIP 
echo $STRANDNESS 
echo $BROAD_PEAKS 
echo $DUPLICATES 
echo $FILE 

filename=$(basename "$FILE" .fastq)

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/BAM_TOPHAT

echo "Aligning" 
echo $filename 

#---Align using TOPHAT ------
tophat -p 4 -G $REFSEQ -o $scratch/$tmp/BAM_TOPHAT/$filename --max-multihits 2 --b2-very-sensitive $BOWTIE2_INDEX/genome $FILE 
gzip $FILE
# Select those with single mapping
samtools view -h -F 4 $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam | grep -v "NH:i:2" | samtools view -b -S - > $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam  
# RM black list
echo "    Remove black list" 
bedtools intersect -a $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam -b $BL -v > $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam 

# Indexing
echo "    Indexing" 
samtools index -b $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam 
# We create the bigwig
echo "    Create scaled BigWig" 
echo "    Scaling" 
factor=$(samtools view -c -F 260 $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam)
scale_f=$(echo "scale=7 ; $Scale_to / $factor" | bc)
# Make big_wig norm to depth
bamCoverage -p 20 -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam -o $OUT_FOLDER/BIGWIG/"$filename".scale_depth.bw 

# Before deleting all we create 
# MAPPING ALL READS IN INITIAL FQ
echo "Reads in fq" 
echo $(cat $FILE | wc -l)/4|bc 

# MAPPING number of reads uniquely mapped
echo "Uniquelty mapping" 
samtools view -h $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam | grep -v "NH:i:1" | wc -l  

# MAPPING multimappers
echo "Multi mapping" 
samtools view -h $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam | grep "NH:i:1" | wc -l  

# MAPPING Cleaned
echo "Multi mapping" 
samtools view -h $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam | wc -l  

# we keep only the final clean file
rm -rf $scratch/$tmp/BAM_TOPHAT/$filename*

source deactivate