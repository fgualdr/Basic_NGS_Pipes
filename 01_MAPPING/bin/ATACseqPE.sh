    #!/bin/bash
# copyrights(c) Francesco Gualdrini
# List of commands run by ATACseq
#PBS -l select=1:ncpus=4:mem=16g
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
echo $FILE_1

rm -rf $OUT_FOLDER/2_ShortStack_align
rm -rf $OUT_FOLDER/4_tagAlign

mkdir $OUT_FOLDER/2_BAM/

tmp=$(basename "$OUT_FOLDER")
scratch=/hpcnfs/scratch/GN/fgualdrini/
mkdir $scratch/$tmp
mkdir $scratch/$tmp/2_BAM
mkdir $scratch/$tmp/1_StartDataTrim

filename=$(basename "$FILE_1" .R1.fastq)
echo $filename 
echo "Alagning ATACseq" 

##-----Trimming---##

trimmomatic PE -threads 4 -phred33 $FILE_1 $FILE_2 $scratch/$tmp/1_StartDataTrim/"$filename".R1.trimPair.fastq $scratch/$tmp/1_StartDataTrim/"$filename".R1.trimUnPair.fastq $scratch/$tmp/1_StartDataTrim/"$filename".R2.trimPair.fastq $scratch/$tmp/1_StartDataTrim/"$filename".R2.trimUnPair.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
rm $scratch/$tmp/1_StartDataTrim/"$filename".R1.trimUnPair.fastq $scratch/$tmp/1_StartDataTrim/"$filename".R2.trimUnPair.fastq

gzip $FILE_1
gzip $FILE_2

##-----Align-----##
echo "Aligning"
echo BOWTIE2 MAPPING LOG $filename &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
bowtie2 --very-sensitive -k 2 -t --phred33 -p 4 -q -x $BOWTIE2_INDEX/genome  -1 $scratch/$tmp/1_StartDataTrim/"$filename".R1.trimPair.fastq  -2 $scratch/$tmp/1_StartDataTrim/"$filename".R2.trimPair.fastq -S $scratch/$tmp/2_BAM/"$filename".temp.sam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

samtools view -bS $scratch/$tmp/2_BAM/"$filename".temp.sam > $scratch/$tmp/2_BAM/"$filename".temp.bam

echo REDS IN  $filename AFTER MAPPING &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

##------- SET FOR HMMRATAC --- need re-testing

 samtools sort $scratch/$tmp/2_BAM/"$filename".temp.bam -o $scratch/$tmp/2_BAM/"$filename".temp.sort.bam 
 samtools index $scratch/$tmp/2_BAM/"$filename".temp.sort.bam  $scratch/$tmp/2_BAM/"$filename".temp.sort.bam.bai
# samtools view -H $scratch/$tmp/2_BAM/"$filename".temp.sort.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > $scratch/$tmp/2_BAM/"$filename".genome.info
# HMMRATAC -b $scratch/$tmp/2_BAM/"$filename".temp.sort.bam -i $scratch/$tmp/2_BAM/"$filename".temp.sort.bam.bai -g $scratch/$tmp/2_BAM/"$filename".genome.info -o $OUT_FOLDER/5_Called_Peaks/"$filename".HMMRATAC -p TRUE --bedgraph TRUE --blacklist $BL

##----- keep reads mapped in proper pairs (-f 2) and we exclude reads unmapped, mate unmapped, and failing quality (-F 524)
samtools view -b -f 3 -F 524 $scratch/$tmp/2_BAM/"$filename".temp.sort.bam > $scratch/$tmp/2_BAM/"$filename".temp.pp.bam

echo REDS IN  $filename AFTER FILTERING MATES &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".temp.pp.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

##----- RM chrM reads-----##
samtools view -h  $scratch/$tmp/2_BAM/"$filename".temp.pp.bam  |  grep -v 'chrM'  |  samtools view -b -  >  $scratch/$tmp/2_BAM/"$filename".nochrM.bam

echo REDS IN  $filename AFTER RM CHRM &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".nochrM.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

# sort BAM ---> index and flagstat as required
samtools sort -m 4G -@ 3 -n -o $scratch/$tmp/2_BAM/"$filename".nochrM.sort.bam $scratch/$tmp/2_BAM/"$filename".nochrM.bam

##----- Remove Black list regions
pairToBed -abam $scratch/$tmp/2_BAM/"$filename".nochrM.sort.bam -b $BL -type neither > $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.bam

echo REDS IN  $filename AFTER FILTERING BLACKLIST &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

##------Cleaning DUP --------##
# Add ms and MC tags for markdup to use later - this is essential to assign the full read
samtools fixmate -m $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.bam $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.bam
# Markdup needs position order
samtools sort -o $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.bam $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.bam
# Finally mark duplicates - and remove them with the -r option
samtools markdup -r $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.bam $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.rmdup.bam

echo REDS IN  $filename AFTER DUP REMOVAL &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt
samtools flagstat $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.rmdup.bam &>> $OUT_FOLDER/LOG/"$filename".MAPPED_COUNTS.txt

mv $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.rmdup.bam $scratch/$tmp/2_BAM/"$filename".rmdup.bam
mv $scratch/$tmp/2_BAM/"$filename".nochrM.sort.clean.fixmate.positionsort.bam $scratch/$tmp/2_BAM/"$filename".bam

##------5' Shift------##
samtools index $scratch/$tmp/2_BAM/"$filename".rmdup.bam
samtools index $scratch/$tmp/2_BAM/"$filename".bam
alignmentSieve -p 4 --ATACshift -b $scratch/$tmp/2_BAM/"$filename".rmdup.bam -o $OUT_FOLDER/3_UM5ShiftData/"$filename".rmdup.bam
alignmentSieve -p 4 --ATACshift -b $scratch/$tmp/2_BAM/"$filename".bam -o $OUT_FOLDER/3_UM5ShiftData/"$filename".bam

##------Check parameters  ------##
# Calling Peaks       
if [ "$ORG" = "mm" ]; then
    macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".rmdup.bam --name="$filename".rmdup.nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g mm --call-summits --keep-dup all
    macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".bam --name="$filename".nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g mm --call-summits --keep-dup all
else
    if [ "$ORG" = "hg" ]; then
        macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".rmdup.bam --name="$filename".rmdup.nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g hs --call-summits --keep-dup all
        macs2 callpeak -t $OUT_FOLDER/3_UM5ShiftData/"$filename".bam --name="$filename".nol --outdir $OUT_FOLDER/5_Called_Peaks/ --format=BAMPE -B -g hs --call-summits --keep-dup all
    fi
fi

# Make BW - with dups                                  
echo Scaling
factor=$(samtools view -c $OUT_FOLDER/3_UM5ShiftData/"$filename".bam)
factor_pe=$(echo "$factor / 2" | bc)
scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)

awk -v var=$scale_f '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.bdg > $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg
sort -k1,1 -k2,2n $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg > $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg
rm $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.bdg
bedGraphToBigWig $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg $ChromSizes_path $OUT_FOLDER/6_bw/"$filename".norm.depth.bw
rm $OUT_FOLDER/5_Called_Peaks/"$filename".nol_treat_pileup.norm.depth.sort.bdg

# Make BW - w/o                               
echo Scaling
factor=$(samtools view -c $OUT_FOLDER/3_UM5ShiftData/"$filename".rmdup.bam)
factor_pe=$(echo "$factor / 2" | bc)
scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)

awk -v var=$scale_f '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.bdg > $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.bdg
sort -k1,1 -k2,2n $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.bdg > $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.sort.bdg
rm $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.bdg
bedGraphToBigWig $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.sort.bdg $ChromSizes_path $OUT_FOLDER/6_bw/"$filename".rmdup.norm.depth.bw
rm $OUT_FOLDER/5_Called_Peaks/"$filename".rmdup.nol_treat_pileup.norm.depth.sort.bdg

source deactivate

rm -rf $scratch/$tmp/1_StartDataTrim/"$filename"* $scratch/$tmp/2_BAM/"$filename"*