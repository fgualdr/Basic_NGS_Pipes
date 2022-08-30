#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ChIPseq
#PBS -l select=1:ncpus=5:mem=20g
#PBS -S /bin/bash
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
echo $CALL_PEAKS

filename=$(basename "$FILE" .fastq)

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/1_ShortStack_align
mkdir $scratch/$tmp/2_BAM

echo $filename 
echo "Alagning ChIPseq" 

#- shortstuck mapping
$Path_to_ShortStack/ShortStack --mismatches 2 --bowtie_m 50 --outdir $scratch/$tmp/1_ShortStack_align/$filename --nohp --bowtie_cores 4 --sort_mem 4G --align_only --mmap u --readfile $FILE --genomefile $BOWTIE_INDEX/genome.fa
#gzip $FILE
##------Filtering Unmap as we use ShortStuck we keep mapped uniquely XY:Z:U and XY:Z:P ##
echo FILTERING
samtools view -h $scratch/$tmp/1_ShortStack_align/$filename/"$filename".bam | grep -v -w "XY:Z:N\|XY:Z:M\|XY:Z:O\|XY:Z:R" > $scratch/$tmp/2_BAM/"$filename".mapped.sam
samtools view -h -bS $scratch/$tmp/2_BAM/"$filename".mapped.sam > $scratch/$tmp/2_BAM/"$filename".mapped.bam
echo "Sort"
samtools sort -m 3G -@ 3 $scratch/$tmp/2_BAM/"$filename".mapped.bam > $scratch/$tmp/2_BAM/"$filename".mapped.sort.bam

##------Cleaning--------##
echo "rm BL"
bedtools intersect -a $scratch/$tmp/2_BAM/"$filename".mapped.sort.bam -b $BL -v > $scratch/$tmp/2_BAM/"$filename".mapped.sort.clean.bam
##------MarkDup--------##
# The first sort can be omitted if the file is already name ordered
samtools sort -m 3G -@ 3 -n -o $scratch/$tmp/2_BAM/"$filename".namesort.bam $scratch/$tmp/2_BAM/"$filename".mapped.sort.clean.bam
# Add ms and MC tags for markdup to use later
samtools fixmate -m $scratch/$tmp/2_BAM/"$filename".namesort.bam $scratch/$tmp/2_BAM/"$filename".fixmate.bam
# Markdup needs position order
samtools sort -m 3G -@ 3 -o $scratch/$tmp/2_BAM/"$filename".positionsort.bam $scratch/$tmp/2_BAM/"$filename".fixmate.bam
# Finally mark duplicates
samtools markdup -r $scratch/$tmp/2_BAM/"$filename".positionsort.bam $OUT_FOLDER/2_BAM/"$filename".clean.bam # ONLY MARKED

echo REDS IN  $filename AFTER BLACK LIST RM &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $OUT_FOLDER/2_BAM/"$filename".clean.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt   
          

# Predict D
echo PREDICTD
macs2 predictd -i $OUT_FOLDER/2_BAM/"$filename".clean.bam --format=BAM -g mm --outdir $OUT_FOLDER/4.1_PREDICTD/ --rfile $filename.R

echo "    Indexing" 
samtools index -b $OUT_FOLDER/2_BAM/"$filename".clean.bam

echo REDS IN  $filename AFTER RM DUP &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $OUT_FOLDER/2_BAM/"$filename".clean.bam&>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

if [ "$CALL_PEAKS" = "YES" ]; then
    # Calling Peaks                                         
    echo CALLING PEAKS NARROW W/O DUP
    if [ "$INPUT_CHIP" = "NA" ]; then
        echo NO INPUT
        if [ "$ORG" = "mm" ]; then
            macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam --name="$filename".filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g mm --extsize 250 --call-summits
        else
            if [ "$ORG" = "hg" ]; then
                macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam --name="$filename".filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g hs --extsize 250 --call-summits
            fi
        fi
    else
        if [ "$ORG" = "mm" ]; then
            macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam -c $INPUT_CHIP/*.clean.filterdup.bam --name="$filename".filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g mm --extsize 250 --call-summits
            rm $OUT_FOLDER/4_Called_Peaks/"$filename".filterdup.nol_control*.bdg  
        else
            if [ "$ORG" = "hg" ]; then
                macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam -c $INPUT_CHIP/*.clean.filterdup.bam --name="$filename".filterdup.nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAM -B -g hs --extsize 250 --call-summits
                rm $OUT_FOLDER/4_Called_Peaks/"$filename".filterdup.nol_control*.bdg 
            fi
        fi
    fi
fi

# Make BW                                               
echo Scaling
factor=$(samtools view -c -F 260 $OUT_FOLDER/2_BAM/"$filename".clean.bam)
scale_f=$(echo "scale=7 ; $Scale_to / $factor" | bc)
bamCoverage -b  $OUT_FOLDER/2_BAM/"$filename".clean.bam -o $OUT_FOLDER/5_bw/"$filename".filterdup.norm.depth.bw --scaleFactor $scale_f --extendReads 300

source deactivate
cd $OUT_FOLDER/4.1_PREDICTD/
source activate EnvR
R_FILE=$OUT_FOLDER/4.1_PREDICTD/$filename.R
Rscript $R_FILE
source deactivate

rm -rf $scratch/$tmp/1_ShortStack_align/$filename* $scratch/$tmp/2_BAM/"$filename"*