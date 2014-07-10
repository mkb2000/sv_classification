#!/bin/bash
# run preprocess on output of Socrates, and merge the results for training but left one for test.

if [ $# -ne 2 ]; then
 echo "Usage: <dir of Socrates results> <output dir for processed data for ML>"
 exit
fi

if [ ! -d $2 ]; then
mkdir $2
fi


for i in $1/*_paired* ; do
#echo $i
#sf="simulated_data/results_Socrates_paired_simulated_chr12_"$i"_simseq_s_long_sc_l25_q5_m5_i95.txt"
    echo $2 $i
    python ./ML_preprocess_for_test_data.py $i $2
done

numOfFiles=`ls ${2}/*.arff | wc -l`
#numOneLeft=`expr ${numOfFiles} - 0`
python ./mergeARFF.py `ls $2/*.arff`

