#!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ChIPseq
#PBS -l select=1:ncpus=8:mem=20g
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
echo $BOWTIE2_INDEX=
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

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/2_BAM

echo $filename 
echo "Alagning ChIPseq" 

source activate EnvATAC 

#- shortstuck mapping
bowtie2 -X 1000 -k 2 -t --phred33 -p 8 -q -x $BOWTIE2_INDEX/genome -1 $FILE_1 -2 $FILE_2 -S $scratch/$tmp/2_BAM/"$filename".temp.sam
gzip $FILE_1
gzip $FILE_2

##------Filtering
# SAM to BAM
samtools view -bS $scratch/$tmp/2_BAM/"$filename".temp.sam > $scratch/$tmp/2_BAM/"$filename".temp.bam

echo REDS IN  $filename AFTER MAPPING &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

# keep reads mapped in proper pairs (-f 2) and we exclude reads unmapped, mate unmapped, and failing quality (-F 524)
samtools view -b -f 2 -F 524 $scratch/$tmp/2_BAM/"$filename".temp.bam > $scratch/$tmp/2_BAM/"$filename".temp.pp.bam

echo REDS IN  $filename AFTER FILTERING MATES &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.pp.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

# -------- Remove multi mappers ----------- ##
samtools view -h $scratch/$tmp/2_BAM/"$filename".temp.pp.bam | grep -v "XS:i" > $scratch/$tmp/2_BAM/"$filename".temp.pp.sam
rm -rf $scratch/$tmp/2_BAM/"$filename".temp.pp.bam
samtools view -h -bS $scratch/$tmp/2_BAM/"$filename".temp.pp.sam > $scratch/$tmp/2_BAM/"$filename".temp.pp.bam

echo REDS IN  $filename AFTER removing multi mappers &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.pp.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt 

# sort BAM ---> index and flagstat as required
samtools sort -m 2G -@ 8 -n -o $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.bam $scratch/$tmp/2_BAM/"$filename".temp.pp.bam

# Remove Black list regions
pairToBed -abam $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.bam -b $BL -type neither > $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.bam

echo REDS IN  $filename AFTER FILTERING BLACKLIST &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
rm $OUT_FOLDER/2_BAM/"$filename".temp.pp.sort.bam

# Add ms and MC tags for markdup to use later - this is essential to assign the full read
samtools fixmate -m $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.bam $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.bam

# Markdup needs position order
samtools sort -o $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.positionsort.bam $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.bam

# Finally mark duplicates - and remove them with the -r option
samtools markdup -r $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.positionsort.bam $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.positionsort.rmdup.bam

echo REDS IN  $filename AFTER DUP REMOVAL &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.positionsort.rmdup.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

mv $scratch/$tmp/2_BAM/"$filename".temp.pp.sort.clean.fixmate.positionsort.rmdup.bam $OUT_FOLDER/2_BAM/"$filename".clean.bam

# Index
samtools index $OUT_FOLDER/2_BAM/"$filename".clean.bam


#############
############# To the Peak calling work out the INPUT
#############

# Predict-D in macs - get the D value and use it for peak calling Specified Mouse mm
macs2 predictd -f BAMPE -i $OUT_FOLDER/2_BAM/"$filename".clean.bam --format=BAMPE -g mm --outdir $OUT_FOLDER/4.1_PREDICTD/ --rfile $filename.R
cd $OUT_FOLDER/4.1_PREDICTD/
Rscript $OUT_FOLDER/4.1_PREDICTD/$filename.R 


if [ "$CALL_PEAKS" = "YES" ]; then
    # Calling Peaks                                         
    echo CALLING PEAKS NARROW W/O DUP
    if [ "$INPUT_CHIP" = "NA" ]; then
        echo NO INPUT
        if [ "$ORG" = "mm" ]; then
            macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam --name="$filename".nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAMPE -B -g mm --broad --broad-cutoff 0.1
        else
            if [ "$ORG" = "hg" ]; then
                macs2 callpeak -t $OUT_FOLDER/2_BAM/"$filename".clean.bam --name="$filename".nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAMPE -B -g hs --broad --broad-cutoff 0.1
            fi
        fi
    else
        if [ "$ORG" = "mm" ]; then
            macs2 callpeak -t $$OUT_FOLDER/2_BAM/"$filename".clean.bam -c $INPUT_CHIP --name="$filename".nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAMPE -B -g mm --broad --broad-cutoff 0.1
            rm $OUT_FOLDER/4_Called_Peaks/"$filename".nol_control*.bdg 
        else
            if [ "$ORG" = "hg" ]; then
                macs2 callpeak -t $$OUT_FOLDER/2_BAM/"$filename".clean.bam -c $INPUT_CHIP --name="$filename".nol --outdir $OUT_FOLDER/4_Called_Peaks/ --format=BAMPE -B -g hs --broad --broad-cutoff 0.1
                rm $OUT_FOLDER/4_Called_Peaks/"$filename".nol_control*.bdg
            fi
        fi
    fi
fi

# Make BW                                        
echo Scaling
factor=$(samtools view -c $OUT_FOLDER/2_BAM/"$filename".clean.bam)
factor_pe=$(echo "$factor / 2" | bc)
scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)
awk -v var=$scale_f '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.bdg > $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg
sort -k1,1 -k2,2n $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg > $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg
rm $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg
bedGraphToBigWig $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg $ChromSizes_path $OUT_FOLDER/5_bw/"$filename".norm.depth.bw
rm $OUT_FOLDER/4_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg

rm -rf $scratch/$tmp/2_BAM/"$filename"*

