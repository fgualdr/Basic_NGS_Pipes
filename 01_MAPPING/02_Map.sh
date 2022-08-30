#!/bin/bash

##copyrights(c) Francesco Gualdrini
##This tool will map RNAseq / ChromRNAseq / Chip-seq / ATAC-seq samples
##RNA-seq are going to be aligned with TopHat
##ChIP-seq and ATAC-seq are going to be aligned with bowtie2
##To work three conda environments are necessary:
##EnvAlignRNAseq = containing tophat2, samtools, bedtools,  ucsc-wigtobigwig , ucsc-liftover , ucsc-bedtobigbed, bowtie, bowtie2,bc, deeptools
##EnvATAC = ucsc-wigtobigwig , ucsc-liftover , ucsc-bedtobigbed  , trimmomatic , samtools , macs2 , bowtie , bedtools , bedops (ChIP-seq will be processed with same tools as ATACseq) also bc, deeptools
## Add Picard to all
##The user has to provide the location of the Bowtie2 and Bowtie1 index
##the reference transcriptome in GFF format (eg mm10, mm9, hg19 etc..)
##black_list file in bed format

echo -e $INITIAL_MESS

while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "02_Map tool to align .fastq to the desire genome producing BAM files."
                        echo " REMEMBER!! The BIGWIG files are going to be normalised by the depth of sequencing"
                        echo " "
                        echo "To work two conda environments are necessary:"
                        echo "1) EnvAlignRNAseq = containing tophat2, samtools, bedtools"
                        echo "2) EnvATAC = ucsc-wigtobigwig , ucsc-liftover , ucsc-bedtobigbed  , trimmomatic , samtools , macs2 , bowtie2 , bedtools , bedops (ChIP-seq will be processed with same tools as ATACseq)"
                        echo " "
                        echo "options:"
                        echo "-h, --help                show brief help"
                        echo "-t, --type=TYPE_AL       specify the type of alignment e.g. RNAseq, ChromRNAseq, ATACseq, ChIP-seq"
                        echo "-o, --output-dir=OUT_FOLDER      specify a directory to store output in must contain the Merged_Fastq folder (Run 01_Merge_fq.R first)"
                        echo "-i, --index-bowtie2=BOWTIE2_INDEX      specify the bowtie2 index full path (folder containing the genome files) this is mandatory for RNAseq and ChromRNAseq"
                        echo "-x, --index-bowtie=BOWTIE_INDEX      specify the bowtie index full path mandatory for ChromRNAseq, ChIP-seq and ATACseq(folder containing the genome files)"
                        echo "-r, --refseq-reference=REFSEQ      specify the Refseq reference file (GFF,GFF3 format)"
                        echo "-b, --black-list=BL      specify the Black list file"
                        echo "-c, --ChromSizes=ChromSizes_path      specify the ChromSizes file"
                        echo "-s, --ShortStack-path=Path_to_ShortStack      specify the path to ShortStuck tool"
                        echo "-n, --chip-input=chip-input      specify the location of the input file"
                        echo "-d, --Starndness=library-strandness   specify if the library has to be treated in a strand specific manner YES or NO"
                        echo "-k, --Broad-peaks=macs2 broad peaks YES NO"
                        echo "-u, --Duplicates=Use duplicates (DEPRECATED - RNAseq are treated with all read counts is impossible if SE to assess PCR dup - Chip/ATAC both files are produced)"
                        echo "-Sm, --Sequencing-mode=Tell the sequencing mode either PE or SE"
                        exit 0
                        ;;
                -t)
                        shift
                        if test $# -gt 0; then
                                export TYPE_AL=$1
                        else
                                echo "no Alignment type specified"
                                exit 1
                        fi
                        shift
                        ;;
                --type*)
                        export TYPE_AL=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -d)
                        shift
                        if test $# -gt 0; then
                                export STRANDNESS=$1
                        else
                                echo "no strandness specified"
                                exit 1
                        fi
                        shift
                        ;;
                --Starndness*)
                        export STRANDNESS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -o)
                        shift
                        if test $# -gt 0; then
                                export OUT_FOLDER=$1
                        else
                                echo "no output dir specified"
                                exit 1
                        fi
                        shift
                        ;;
                --output-dir*)
                        export OUT_FOLDER=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -i)
                        shift
                        if test $# -gt 0; then
                                export BOWTIE2_INDEX=$1
                        else
                                echo "no bowtie2 Index specified; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --index-bowtie2*)
                        export BOWTIE2_INDEX=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -x)
                        shift
                        if test $# -gt 0; then
                                export BOWTIE_INDEX=$1
                        else
                                echo "no bowtie2 Index specified; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --index-bowtie*)
                        export BOWTIE_INDEX=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -r)
                        shift
                        if test $# -gt 0; then
                                export REFSEQ=$1
                        else
                                echo "no transcriptome reference specified (Refseq in GFF GFF3 format); if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --refseq-reference*)
                        export REFSEQ=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -b)
                        shift
                        if test $# -gt 0; then
                                export BL=$1
                        else
                                echo "no Black list file specified visit https://www.encodeproject.org/search/?type=Annotation; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --black-list*)
                        export BL=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -c)
                        shift
                        if test $# -gt 0; then
                                export ChromSizes_path=$1
                        else
                                echo "no chromosome sizes specified; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --ChromSizes*)
                        export ChromSizes_path=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -s)
                        shift
                        if test $# -gt 0; then
                                export Path_to_ShortStack=$1
                        else
                                echo "no shortstuck path specified; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --ShortStack-path*)
                        export Path_to_ShortStack=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                 -n)
                        shift
                        if test $# -gt 0; then
                                export INPUT_CHIP=$1
                        else
                                echo "no input control; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --chip-input*)
                        export INPUT_CHIP=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -k)
                        shift
                        if test $# -gt 0; then
                                export BROAD_PEAKS=$1
                        else
                                echo "YES for broad peak macs2 or NO if narrow required only for ChIP-seq alignment and peak calling"
                                exit 1
                        fi
                        shift
                        ;;
                --Broad-peaks*)
                        export BROAD_PEAKS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

                -Pa)
                        shift
                        if test $# -gt 0; then
                                export PPath=$1
                        else
                                echo "Provide path for mapping script"
                                exit 1
                        fi
                        shift
                        ;;
                --Path*)
                        export PPath=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;


                -u)
                        shift
                        if test $# -gt 0; then
                                export DUPLICATES=$1
                        else
                                echo "YES the peak calling will be done with DUPLICATES usefull for conditions with repeats elements"
                                exit 1
                        fi
                        shift
                        ;;
                --Duplicates*)
                        export DUPLICATES=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

                -Sm)
                        shift
                        if test $# -gt 0; then
                                export SEQ_MODE=$1
                        else
                                echo "YES the peak calling will be done with DUPLICATES usefull for conditions with repeats elements"
                                exit 1
                        fi
                        shift
                        ;;
                --Sequencing-mode*)
                        export SEQ_MODE=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

                -Cp)
                        shift
                        if test $# -gt 0; then
                                export CALL_PEAKS=$1
                        else
                                echo "YES the peak calling will be done with DUPLICATES usefull for conditions with repeats elements"
                                exit 1
                        fi
                        shift
                        ;;
                --Call-peaks*)
                        export CALL_PEAKS=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                -Org)
                        shift
                        if test $# -gt 0; then
                                export ORG=$1
                        else
                                echo "Organism can either be mm for mus musculus or hs for Human. This flag determine the genome size for MCAS2 peak calling only for the Chipseq and ATACseq mode"
                                exit 1
                        fi
                        shift
                        ;;
                --Organism*)
                        export ORG=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;

                *)
                        break
                        ;;
        esac
done

echo CHECK PATH
echo $PATH

## VARS
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
echo BROAD_PEAKS $BROAD_PEAKS 
echo $DUPLICATES 
echo $PPath
echo $SEQ_MODE
echo $CALL_PEAKS
echo $ORG

## SORT by type

SAMP_PATH=${OUT_FOLDER}/Merged_Fastq

if [ "$TYPE_AL" = "RNAseq" ]; then

        echo "RNAseq processing"
        #mkdir $OUT_FOLDER/BAM_TOPHAT
        mkdir $OUT_FOLDER/CLEAN_BAM
        mkdir $OUT_FOLDER/BIGWIG
        mkdir $OUT_FOLDER/LOG
        # Loop through each sample and lounch the pipe
        if [ "$SEQ_MODE" = "SE" ]; then
                echo "SE mode"
                for f in $SAMP_PATH/*.fastq
                do
                        FILE=$f
                        filename=$(basename "$f" .fastq)
                        VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE="$FILE"
                        if [ "$STRANDNESS" = "NO" ]; then 
                                echo SEND JOB $filename.RNAseqNStr.map.fgualdrini
                                qsub -N $filename.RNAseqNStr.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/RNAseq_NOStr.sh  
                        else         
                                echo SEND JOB $filename.RNAseqStr.map.fgualdrini
                                qsub -N $filename.RNAseqStr.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/RNAseq_Str.sh
                        fi
                done
        else
                if [ "$SEQ_MODE" = "PE" ]; then
                        echo "PE mode"
                        for f in $SAMP_PATH/*.R1.fastq
                        do
                                filename=$(basename "$f" .R1.fastq)
                                FILE_1=$SAMP_PATH/"$filename".R1.fastq
                                FILE_2=$SAMP_PATH/"$filename".R2.fastq
                                VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE_1="$FILE_1",FILE_2="$FILE_2"
                                if [ "$STRANDNESS" = "NO" ]; then 
                                        echo SEND JOB $filename.RNAseqNStrPE.map.fgualdrini
                                        qsub -N $filename.RNAseqNStrPE.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/RNAseq_NOStr_PE.sh  
                                else         
                                        echo SEND JOB $filename.RNAseqStrPE.map.fgualdrini
                                        qsub -N $filename.RNAseqStrPE.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/RNAseq_Str_PE.sh
                                fi
                        done
                fi

        fi
else
        if [ "$TYPE_AL" = "ChromRNAseq" ]; then
                echo "ChromRNAseq processing"
                #mkdir $OUT_FOLDER/BAM_TOPHAT
                #mkdir $OUT_FOLDER/SHORT_STUCK_UNMAPPED
                mkdir $OUT_FOLDER/CLEAN_BAM
                mkdir $OUT_FOLDER/BIGWIG
                mkdir $OUT_FOLDER/LOG
                if [ "$SEQ_MODE" = "SE" ]; then
                        for f in $SAMP_PATH/*.fastq
                        do  
                                FILE=$f
                                filename=$(basename "$f" .fastq)
                                VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE="$FILE"
                                if [ "$STRANDNESS" = "NO" ]; then
                                        echo SEND JOB $filename.ChARNStr.map.fgualdrini
                                        qsub -N $filename.ChARNStr.map.fgualdrini -v $VAR  -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChARseq_NOStr.sh  
                                else
                                        echo SEND JOB $filename.ChARStr.map.fgualdrini
                                        qsub -N $filename.ChARStr.map.fgualdrini -v $VAR  -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChARseq_Str.sh
                                fi
                        done
                else
                        if [ "$SEQ_MODE" = "PE" ]; then
                                for f in $SAMP_PATH/*.R1.fastq
                                do
                                        filename=$(basename "$f" .R1.fastq)
                                        FILE_1=$SAMP_PATH/"$filename".R1.fastq
                                        FILE_2=$SAMP_PATH/"$filename".R2.fastq
                                        VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE_1="$FILE_1",FILE_2="$FILE_2"
                                        if [ "$STRANDNESS" = "NO" ]; then 
                                                echo SEND JOB $filename.ChARNStrPE.map.fgualdrini
                                                qsub -N $filename.ChARNStrPE.map.fgualdrini -v $VAR  -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChARseq_NOStrPE.sh 
                                        else         
                                                echo SEND JOB $filename.ChARStr.map.fgualdrini
                                                qsub -N $filename.ChARStr.map.fgualdrini -v $VAR  -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChARseq_StrPE.sh
                                        fi
                                done
                        fi
                fi

        else
        	if [ "$TYPE_AL" = "ATACseq" ]; then
                echo "ATAC processing"
                        #mkdir $OUT_FOLDER/1_StartDataTrim
                        #mkdir $OUT_FOLDER/2_ShortStack_align
        		mkdir $OUT_FOLDER/3_UM5ShiftData
        		mkdir $OUT_FOLDER/4_tagAlign
                        mkdir $OUT_FOLDER/4.1_PREDICTD
        		mkdir $OUT_FOLDER/5_Called_Peaks
        		mkdir $OUT_FOLDER/6_bw
                        mkdir $OUT_FOLDER/LOG                        
                        args=()
                        if [ "$SEQ_MODE" = "SE" ]; then
                                for f in $SAMP_PATH/*.fastq
                                do      
                                        echo SEND JOB $filename.ATAC.map.fgualdrini
                                        FILE=$f
                                        filename=$(basename "$f" .fastq)
                                        VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE="$FILE"
                                        qsub -N $filename.ATAC.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ATACseq.sh
                                done
                        else
                                if [ "$SEQ_MODE" = "PE" ]; then
                                        for f in $SAMP_PATH/*.R1.fastq
                                        do
                                                echo SEND JOB $filename.ATACPE.map.fgualdrini
                                                filename=$(basename "$f" .R1.fastq)
                                                FILE_1=$SAMP_PATH/"$filename".R1.fastq
                                                FILE_2=$SAMP_PATH/"$filename".R2.fastq
                                                VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE_1="$FILE_1",FILE_2="$FILE_2"
                                                qsub -N $filename.ATACPE.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ATACseqPE.sh
                                        done
                                fi
                        fi
        	else
                        if  [ "$TYPE_AL" = "Chipseq" ]; then
                                echo "Chipseq processing"    
                                #mkdir $OUT_FOLDER/1_ShortStack_align
                                mkdir $OUT_FOLDER/2_BAM
                                mkdir $OUT_FOLDER/4_Called_Peaks
                                mkdir $OUT_FOLDER/5_bw
                                mkdir $OUT_FOLDER/4.1_PREDICTD
                                mkdir $OUT_FOLDER/LOG
                                if [ "$SEQ_MODE" = "SE" ]; then
                                        for f in $SAMP_PATH/*.fastq
                                        do  
                                                FILE=$f
                                                filename=$(basename "$f" .fastq)
                                                VAR=Path_to_ShortStack="$Path_to_ShortStack",ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE="$FILE",CALL_PEAKS="$CALL_PEAKS"
                                                if [ "$BROAD_PEAKS" != "YES" ]; then
                                                        echo SEND JOB $filename.ChIPN.map.fgualdrini
                                                        qsub -N $filename.ChIPN.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChIPseqNarrow.sh
                                                else
                                                        if [ "$BROAD_PEAKS" == "YES" ]; then
                                                                echo SEND JOB $filename.ChIPB.map.fgualdrini
                                                                qsub -N $filename.ChIPB.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.log -e $OUT_FOLDER/LOG/$filename.Error.log $PPath/ChIPseqBroad.sh
                                                        fi
                                                fi
                                        done
                                else
                                        if [ "$SEQ_MODE" = "PE" ]; then
                                                for f in $SAMP_PATH/*.R1.fastq
                                                do  
                                                        filename=$(basename "$f" .R1.fastq)
                                                        FILE_1=$SAMP_PATH/"$filename".R1.fastq
                                                        FILE_2=$SAMP_PATH/"$filename".R2.fastq
                                                        VAR=ORG="$ORG",TYPE_AL="$TYPE_AL",OUT_FOLDER="$OUT_FOLDER",BOWTIE2_INDEX="$BOWTIE2_INDEX",BOWTIE_INDEX="$BOWTIE_INDEX",REFSEQ="$REFSEQ",BL="$BL",ChromSizes_path="$ChromSizes_path",INPUT_CHIP="$INPUT_CHIP",STRANDNESS="$STRANDNESS",BROAD_PEAKS="$BROAD_PEAKS",DUPLICATES="$DUPLICATES",FILE_1="$FILE_1",FILE_2="$FILE_2",CALL_PEAKS="$CALL_PEAKS"
                                                        if [ "$BROAD_PEAKS" = "NO" ]; then
                                                                echo SEND JOB $filename.ChIPN.map.fgualdrini
                                                                args+=$(qsub -N $filename.ChIPNPE.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.PE.log -e $OUT_FOLDER/LOG/$filename.ErrorPE.log $PPath/ChIPseqNarrowPE.sh)
                                                        else
                                                                echo SEND JOB $filename.ChIPB.map.fgualdrini
                                                                args+=$(qsub -N $filename.ChIPBPE.map.fgualdrini -v $VAR -o $OUT_FOLDER/LOG/$filename.PE.log -e $OUT_FOLDER/LOG/$filename.ErrorPE.log $PPath/ChIPseqBroadPE.sh)
                                                                
                                                        fi
                                                done
                                        fi
                                fi
                        fi
                fi
        fi
fi