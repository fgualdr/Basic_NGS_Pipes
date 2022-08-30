#!/bin/bash

echo -e $INITIAL_MESS

while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "Trim adapters in SE mode or PE mode"
                        echo "To work two conda environments are necessary:"
                        echo "1) EnvCutadap"
                        echo "2) EnvATAC = ucsc-wigtobigwig , ucsc-liftover , ucsc-bedtobigbed  , trimmomatic , samtools , macs2 , bowtie2 , bedtools , bedops (ChIP-seq will be processed with same tools as ATACseq)"
                        echo " "
                        echo "options:"
                        echo "-o, --output-dir=OUT_FOLDER      specify a directory to store output in must contain the Merged_Fastq folder (Run 01_Merge_fq.R first)"
                        echo "-s, --seq-mode=SEQ_MODE      specify SE or PE"
                        exit 0
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
                 -s)
                        shift
                        if test $# -gt 0; then
                                export SEQ_MODE=$1
                        else
                                echo "no bowtie2 Index specified; if you are running other methods yes type -n NO"
                                exit 1
                        fi
                        shift
                        ;;
                --seq-mode*)
                        export SEQ_MODE=`echo $1 | sed -e 's/^[^=]*=//g'`
                        shift
                        ;;
                *)
                        break
                        ;;
        esac
done

echo $OUT_FOLDER
echo $SEQ_MODE

#Trim adapters 
source activate EnvCutadap
SAMP_PATH=$OUT_FOLDER/Merged_Fastq
resource=/hpcnfs/data/GN2/fgualdrini/UsefullData/DB/

if [ "$SEQ_MODE" = "SE" ]; then

    for f in $SAMP_PATH/*.fastq;
    do
        filename=$(basename "$f" .fastq)
        echo $filename
        #Adapter Trimming
        bbduk.sh -Xmx1g threads=7 in=$f out=$SAMP_PATH/"$filename".clean.fastq ref=$resource/adapters.fa ktrim=r k=28 mink=13 hdist=1 tbo tpe
        rm $f
        mv $SAMP_PATH/"$filename".clean.fastq $SAMP_PATH/"$filename".fastq
    done

    else

        if [ "$SEQ_MODE" = "PE" ]; then

            for f in $SAMP_PATH/*_R1.fastq;
            do
                filename=$(basename "$f" _R1.fastq)
                echo $filename
                #Adapter Trimming
                bbduk.sh -Xmx1g threads=7 in1=$SAMP_PATH/"$filename"_R1.fastq in2=$SAMP_PATH/"$filename"_R2.fastq out1=$SAMP_PATH/"$filename".clean.R1.fastq out2=$SAMP_PATH/"$filename".clean.R2.fastq ref=$resource/adapters.fa ktrim=r k=28 mink=13 hdist=1 tbo tpe
                rm $SAMP_PATH/"$filename"_R1.fastq $SAMP_PATH/"$filename"_R2.fastq
                mv $SAMP_PATH/"$filename".clean.R1.fastq $SAMP_PATH/"$filename".R1.fastq 
                mv $SAMP_PATH/"$filename".clean.R2.fastq $SAMP_PATH/"$filename".R2.fastq
            done
        

        fi

fi


conda deactivate
