#!/bin/bash

if [ $# -ne 2 ]; then
 echo "Usage: <batch no> <outdir>"
 exit
fi
basedir=$2"/"
out=${basedir}"simulated_chr12_"$1
java -cp src VariationRandomization data/chr12.fa.fai 1500000 $out hom

java -cp src/commons-lang3-3.2.1.jar:src/picard-tools-1.77/picard-1.77.jar:src/picard-tools-1.77/sam-1.77.jar:src ReferenceSimulation data/chr12.fa $out".fa"
