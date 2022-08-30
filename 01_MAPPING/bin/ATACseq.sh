#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ATACseq
#PBS -l select=1:ncpus=4:mem=30g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -l maxarray_5=1

source activate EnvATAC

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
echo $ORG


tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/1_StartDataTrim
mkdir $scratch/$tmp/2_ShortStack_align

filename=$(basename "$FILE" .fastq)
echo $filename
echo "ATACseq processing"

echo "Trimming"

##------trim------##
trimmomatic SE -threads 4 -phred33 $FILE $scratch/$tmp/1_StartDataTrim/"$filename".trim.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:38
gzip $FILE
##-----Align-----##
echo "Aligning"
$Path_to_ShortStack/ShortStack --mismatches 2 --bowtie_m 50 --outdir $scratch/$tmp/2_ShortStack_align/$filename --nohp --bowtie_cores 4 --sort_mem 7G --align_only --mmap u --readfile $scratch/$tmp/1_StartDataTrim/"$filename".trim.fastq --genomefile $BOWTIE_INDEX/genome.fa

##------Filtering Unmap as we use ShortStuck we keep mapped uniquely XY:Z:U and XY:Z:P ##
echo FILTERING
samtools view -h $scratch/$tmp/2_ShortStack_align/$filename/"$filename".trim.bam | grep -v -w "XY:Z:N\|XY:Z:M\|XY:Z:O\|XY:Z:R" > $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sam
samtools view -h -bS $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sam > $scratch/$tmp/2_ShortStack_align/"$filename".mapped.bam
samtools sort -m 7G -@ 4 $scratch/$tmp/2_ShortStack_align/"$filename".mapped.bam > $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sort.bam

# BL remove
bedtools intersect -a $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sort.bam -b $BL -v > $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sort.clean.bam

##------5' Shift------##
## Converts to BED then
echo SHIFTING
bam2bed < $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sort.clean.bam > $scratch/$tmp/2_ShortStack_align/"$filename".temp.bed
# -- shift reads +4 on the plus strand and -5 on the minus strand (as BED are 0 based) -- WARNING issues with reads at the ends of chromosomes
bedtools shift -i $scratch/$tmp/2_ShortStack_align/"$filename".temp.bed -p 4 -m -5 -g $ChromSizes_path > $scratch/$tmp/2_ShortStack_align/"$filename".temp.5Shift0.bed
sort -k 1,1 -k2,2n $scratch/$tmp/2_ShortStack_align/"$filename".temp.5Shift0.bed > $scratch/$tmp/2_ShortStack_align/"$filename".temp.5Shift.bed
bedToBam -i $scratch/$tmp/2_ShortStack_align/"$filename".temp.5Shift.bed -g $ChromSizes_path > $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam

# Make it without chrM
samtools index $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam
samtools idxstats $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam | cut -f 1 | grep -v 'chrM' | xargs samtools view -b $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam > $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.nochrM.bam

##------Check parameters  ------##
echo CALLING PEAKS

if [ "$ORG" = "mm" ]; then
    macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.nochrM.bam --name="$filename".nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAM -B -g mm --nomodel --extsize 146 --nolambda --call-summits --keep-dup 'all'
else
    if [ "$ORG" = "hg" ]; then
        macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.nochrM.bam --name="$filename".nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAM -B -g hs --nomodel --extsize 146 --nolambda --call-summits --keep-dup 'all'
    fi
fi

echo Scaling
factor=$(samtools view -c -F 260 $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.nochrM.bam)
scale_f=$(echo "scale=7 ; $Scale_to / $factor" | bc)
awk -v var=$scale_f '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.bdg > $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg
wigToBigWig -clip $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg $ChromSizes_path $OUT_FOLDER/6_bw/"$filename".norm.depth.bw
rm $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg

echo $filename
echo NUMBER OF READS IN .fastq
echo $(cat $FILE |wc -l)/4|bc
gzip $FILE
echo NUMBER OF READS MAPPED
samtools view -c $scratch/$tmp/2_ShortStack_align/"$filename".mapped.sort.bam

echo NUMBER OF READS MAPPED w/o BL
samtools view -c $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam

echo NUMBER OF READS IN BAM w/o ChrM
samtools view -c $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.nochrM.bam

rm -rf $scratch/$tmp/1_StartDataTrim/"$filename"* $scratch/$tmp/2_ShortStack_align/$filename*
rm $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam
rm $OUT_FOLDER/3_UM5ShiftData/"$filename".5Shift.bam.bai

source deactivate