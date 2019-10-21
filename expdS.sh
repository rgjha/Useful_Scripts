#!/bin/bash
cd ~/Desktop/Runs
for d in SLNC* ; do
if ls $d/results/exp_dS.dat 1> /dev/null 2>&1 ; then
echo -ne "$d "
temp=`awk '{print int(100 * $1)}' $d/results/accP.dat`
echo -ne "had $temp% acceptance with <e^(-dS)> = "
temp=`awk '{print $1, "+/-", $2}' $d/results/exp_dS.dat`
echo -ne "$temp "
temp=`awk 'OFMT = "%.2f" {print sqrt((1.0-$1)^2) / $2}' $d/results/exp_dS.dat`
echo "(${temp}sigma from 1)"
fi
done



