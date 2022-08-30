#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ChIPseq
#PBS -l select=1:ncpus=5:mem=20g
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
echo $CALL_PEAKS


filename=$(basename "$FILE" .fastq)

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/$filename

echo $filename 
echo "Alagning ChIPseq" 

#-  mapping
bowtie2 -k 2 -t --phred33 -p 4 -q -x $BOWTIE2_INDEX/genome -U $FILE -S $scratch/$tmp/$filename/$filename.sam 
gzip $FILE

##----- SAM to BAM -------##
samtools view -h -S -b -o $scratch/$tmp/$filename/$filename.bam $scratch/$tmp/$filename/$filename.sam 
echo REDS IN  $filename AFTER mapping &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/$filename/$filename.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

# -------- Remove multi mappers ----------- ##
samtools view -h $scratch/$tmp/$filename/$filename.bam  | grep -v "XS:i" > $scratch/$tmp/$filename/$filename.sam
samtools view -h -bS $scratch/$tmp/$filename/$filename.sam > $scratch/$tmp/$filename/$filename.rmmm.bam

echo REDS IN  $filename AFTER removing multi mappers &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/$filename/$filename.rmmm.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

# -------- Remove poor quality mappers and unmapped ----------- ##
samtools view -h -b -F 4 -q 30 -@ 5 -q 1 -o $scratch/$tmp/$filename/$filename.filter.bam $scratch/$tmp/$filename/$filename.rmmm.bam
echo REDS IN  $filename AFTER  Remove poor quality mappers and unmapped &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/$filename/$filename.filter.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 
#---------------------------------------
# Sorting, Fixmate and Duplicate removal
#---------------------------------------
samtools sort -n -O BAM -@ 5 $scratch/$tmp/$filename/$filename.filter.bam | samtools fixmate -m -@ 5 - - | samtools sort -O BAM -@ 5 - | samtools markdup -r -S -@ 5 - $scratch/$tmp/$filename/$filename.filter.sorted_duprmv.bam

echo REDS IN  $filename AFTER  Remove poor quality mappers and unmapped &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/$filename/$filename.filter.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

##------BL remove--------##
echo "rm BL"
bedtools intersect -a $scratch/$tmp/$filename/$filename.filter.sorted_duprmv.bam -b $BL -v > $OUT_FOLDER/2_BAM/"$filename".clean.bam
#----------
# Indexing
#----------
samtools index $OUT_FOLDER/2_BAM/"$filename".clean.bam $OUT_FOLDER/2_BAM/"$filename".clean.bam.bai

echo REDS IN  $filename AFTER  BL remove  &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $OUT_FOLDER/2_BAM/"$filename".clean.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

# Predict D
echo PREDICTD
macs2 predictd -i $OUT_FOLDER/2_BAM/"$filename".clean.bam --format=BAM -g mm --outdir $OUT_FOLDER/4.1_PREDICTD/ --rfile $filename.R


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
bamCoverage -b  $OUT_FOLDER/2_BAM/"$filename".clean.bam -o $OUT_FOLDER/5_bw/"$filename".filterdup.norm.depth.bw --scaleFactor $scale_f --extendReads 250


source deactivate
cd $OUT_FOLDER/4.1_PREDICTD/
source activate EnvR
R_FILE=$OUT_FOLDER/4.1_PREDICTD/$filename.R
Rscript $R_FILE
source deactivate

rm -rf $scratch/$tmp/$filename/