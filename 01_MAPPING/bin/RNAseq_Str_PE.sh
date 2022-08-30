#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by RNAseq_Strand
# Reads are going to be mapped using tophat and big wig strand specific created

#PBS -l select=1:ncpus=7:mem=40g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V
#PBS -l maxarray_5=1

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
echo $FILE_1
echo $FILE_2

filename=$(basename "$FILE_1" .R1.fastq)
echo $filename 

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/BAM_TOPHAT

source activate EnvAlignRNAseq
#---Align using TOPHAT ------
# Now map with tophat using the reference
tophat2 -p 7 --no-coverage-search --no-discordant -G $REFSEQ -o $scratch/$tmp/BAM_TOPHAT/$filename --max-multihits 2 --b2-very-sensitive -r 200 $BOWTIE2_INDEX/genome $FILE_1 $FILE_2
gzip $FILE_2
gzip $FILE_1

echo REDS IN  $filename AFTER MAPPING &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

# Select Clean the bam file
samtools view -h -Sb -f 2 -F 524 $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam > $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam
rm $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits.bam

echo REDS IN  $filename AFTER FILTERING MATES &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

# sort BAM ---> index and flagstat as required
samtools sort -m 4G -@ 5 -n -o $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.sort.bam $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam
rm $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.bam

# Remove Black list regions
pairToBed -abam $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.sort.bam -b $BL -type neither > $scratch/$tmp/"$filename".temp.clean.bam

echo REDS IN  $filename AFTER FILTERING BLACKLIST &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/"$filename".temp.clean.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
rm $scratch/$tmp/BAM_TOPHAT/$filename/accepted_hits_UM.sort.bam

# Add ms and MC tags for markdup to use later - this is essential to assign the full read
samtools fixmate -m $scratch/$tmp/"$filename".temp.clean.bam $scratch/$tmp/"$filename".temp.clean.fixmate.bam
rm $scratch/$tmp/"$filename".temp.clean.bam
# Markdup needs position order
samtools sort -m 4G -@ 5 -o $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.bam $scratch/$tmp/"$filename".temp.clean.fixmate.bam
rm $scratch/$tmp/"$filename".temp.clean.fixmate.bam
# Finally mark duplicates - do not remove them
samtools markdup -@ 5 $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.bam $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.rmdup.bam

echo REDS IN  $filename AFTER DUP REMOVAL &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.rmdup.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
rm $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.bam

mv $scratch/$tmp/"$filename".temp.clean.fixmate.positionsort.rmdup.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam

samtools index $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam

# We create the bigwig
echo "    Create scaled BigWig" 
echo "    Scaling" 
factor=$(samtools view -c $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam)
factor_pe=$(echo "$factor / 2" | bc)
scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)

# include reads that are 2nd in a pair (128);
# exclude reads that are mapped to the reverse strand (16)
 samtools view -b -f 128 -F 16 $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam > $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd1.bam
# exclude reads that are mapped to the reverse strand (16) and
# first in a pair (64): 64 + 16 = 80
 samtools view -b -f 80 $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam > $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd2.bam
# combine the temporary files
 samtools merge -f $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fw.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd2.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd1.bam
# index the filtered BAM file
 samtools index $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fw.bam
# run bamCoverage
 bamCoverage -p 7 -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fw.bam -o $OUT_FOLDER/BIGWIG/"$filename".scale_depth.fwd.bw
# remove the temporary files
 rm $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd1.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.fwd2.bam

# include reads that map to the reverse strand (128)
# and are second in a pair (16): 128 + 16 = 144
 samtools view -b -f 144 $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam > $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev1.bam
# include reads that are first in a pair (64), but
# exclude those ones that map to the reverse strand (16)
 samtools view -b -f 64 -F 16 $OUT_FOLDER/CLEAN_BAM/"$filename".clean.bam > $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev2.bam
# merge the temporary files
 samtools merge -f $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev2.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev1.bam
# index the merged, filtered BAM file
 samtools index $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev.bam
# run bamCoverage
bamCoverage -p 7 -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev.bam -o $OUT_FOLDER/BIGWIG/"$filename".scale_depth.rev.bw
# remove temporary files
rm $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev1.bam $OUT_FOLDER/CLEAN_BAM/"$filename".clean.rev2.bam

rm -rf $scratch/$tmp/BAM_TOPHAT/"$filename"* $scratch/$tmp/"$filename"*


source deactivate
