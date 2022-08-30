#!/bin/bash

#PBS -l select=1:ncpus=8:mem=20g
#PBS -S /bin/bash
#PBS -M francesco.gualdrini@ieo.it
#PBS -m abe
#PBS -j oe
#PBS -V


echo "Below the variables you provided:"
echo " "
echo $PEAK_FILE
echo $NORM_FILE
echo $OUT_FOLDER
echo $ChromSizes_path
echo $DUP

source activate EnvATAC
mkdir $OUT_FOLDER
cd $OUT_FOLDER


count=`ls -1 $PEAK_FILE/*.rmdup.nol_treat_pileup.bdg 2>/dev/null | wc -l`
if [ $count != 0 ]
then 
    DUP_RM=$(ls $PEAK_FILE/*.rmdup.nol_treat_pileup.bdg  | grep -v "ALL_MERGED" )
    ALL_READS=$(ls $PEAK_FILE/*.nol_treat_pileup.bdg | grep -v ".rmdup." | grep -v "ALL_MERGED" )
else
    DUP_RM=$(ls $PEAK_FILE/*.nol_treat_pileup.bdg  | grep -v "ALL_MERGED" )
    ALL_READS=$(ls $PEAK_FILE/*.nol_treat_pileup.bdg | grep -v ".rmdup." | grep -v "ALL_MERGED" )
fi

if [ "$DUP" = "YES" ]; then
    for f in $ALL_READS
    do
        filename=$(basename "$f" .bdg)
        identifier=$(basename -s .rmdup.nol_treat_pileup.bdg $(basename -s .nol_treat_pileup.bdg "$f"))
        echo $f
        echo $filename
        echo $identifier
        scale_factor=$(awk -v var="^$identifier" '$0~var{print $NF}' $NORM_FILE)
        echo $scale_factor
        if [ -z "$scale_factor" ]
        then
            echo SKIP
        else
            awk -v var=$scale_factor '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $f > $OUT_FOLDER/"$filename".scaled.bdg
            wigToBigWig -clip $OUT_FOLDER/"$filename".scaled.bdg $ChromSizes_path $OUT_FOLDER/"$filename".scaled.bw
            rm $OUT_FOLDER/"$filename".scaled.bdg $f
        fi
    done
else
    for f in $DUP_RM
    do
        filename=$(basename "$f" .bdg)
        identifier=$(basename -s .rmdup.nol_treat_pileup.bdg $(basename -s .nol_treat_pileup.bdg "$f"))
        echo $f
        echo $filename
        echo $identifier
        scale_factor=$(awk -v var="^$identifier" '$0~var{print $NF}' $NORM_FILE)
        echo $scale_factor
        if [ -z "$scale_factor" ]
        then
            echo SKIP
        else
            awk -v var=$scale_factor '{print $1,$2,$3, $4*var+0}' OFMT="%.0f" $f > $OUT_FOLDER/"$filename".scaled.rmdup.bdg
            wigToBigWig -clip $OUT_FOLDER/"$filename".scaled.rmdup.bdg $ChromSizes_path $OUT_FOLDER/"$filename".scaled.rmdup.bw
            rm $OUT_FOLDER/"$filename".scaled.rmdup.bdg
        fi
    done
fi

source deactivate

