#!/bin/bash
if [ $# -ne 2 ]; then
 echo "Usage: <total batch no> <outdir>"
 exit
fi
basedir=$2"/"
for i in `seq 1 $1` ; do
 f=${basedir}"simulated_chr12_"$i"_simseq_s.bam"
 basename="./"
${basename}/Socrates-master/Socrates all --bowtie2_db ${basename}/data/chr12 ${f}
done

dropbox="/home/kmo/Desktop/Dropbox/"$2
if [ ! -d ${dropbox} ]; then
mkdir ${dropbox}
fi
cp ./results_Socrates_*.txt ${dropbox}
cp ${basedir}/simulated_chr12_*.fa_breakpoints.txt ${dropbox}

mv ./results_Socrates_*.txt ${basedir}
mv ./simulated_*.bam ${basedir}
mv ./simulated_*.bai ${basedir}
mv ./simulated_*.metrics ${basedir}
mv ./simulated_*.fastq.gz ${basedir}

