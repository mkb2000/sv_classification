#!/bin/bash
if [ $# -ne 2 ]; then
 echo "Usage: <total batch no> <outdir>"
 exit
fi


if [ ! -d $2 ]; then
mkdir $2
fi


for i in `seq 1 $1`;do
./script1.sh $i $2
done
./script2_with_ill_err.sh $1 $2
./script3.sh $1 $2

