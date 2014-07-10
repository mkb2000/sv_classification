#!/bin/bash
if [ $# -ne 2 ]; then
 echo "Usage: <total batch no> <outdir>"
 exit
fi
basedir=$2"/"
for i in `seq 1 $1` ; do
 in=${basedir}"simulated_chr12_"$i".fa_reference.fa"
 out=${basedir}"simulated_chr12_"$i"_simseq.sam"
 java  -jar -Xmx2048m ./src/SimSeq-master/SimSeqNBProject/store/SimSeq.jar -1 100 -2 100 --error ./src/SimSeq-master/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt --error2 ./src/SimSeq-master/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt -l 300 -s 30 -n 30000000 -r $in -o $out -u 0.01
 base=${basedir}"simulated_chr12_"$i"_simseq"
 java -Xmx3g -XX:-UseGCOverheadLimit -jar ./src/picard-tools-1.77/SamToFastq.jar I=$out FASTQ=$base"_R1.fastq" F2=$base"_R2.fastq" VALIDATION_STRINGENCY=SILENT
bowtie2 --local -p 8 -x ./data/chr12 -1 $base"_R1.fastq" -2 $base"_R2.fastq" -S $out
samtools view -bS -o $base".bam" $out; samtools sort $base".bam" $base"_s"; samtools index $base"_s.bam"
 rm  $out $base".bam" $base"_R1.fastq" $base"_R2.fastq"
done
