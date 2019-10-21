#!/bin/bash
# Rename all files by moving them back by specified count 

if [[ $# -ne 3 ]] ; then
echo '*** Forgetting arguments !! Order : "FIRST"  "LAST"  "SKIP" '
echo 'First is where you want to start the new numbering from, i.e if 100 becomes 0 then FIRST=100, LAST=5000 something, SKIP=10'
exit 0
fi

first=$1
last=$2
skip=$3
back=$1

for(( i=$first ; $i<$last ; i+=$skip )); do
iter=$[$iter + 1]
next=$[$i + $skip]
new=$[$i - $back]
new_next=$[$i - $back + $skip]
mv ./Out/out.$i-$next ./Out/out.$new-$new_next
mv ./Out/eig.$i ./Out/eig.$new
mv ./Out/corr.$i ./Out/corr.$new
mv ./Configs/gauge.$i ./Configs/gauge.$new
mv ./Configs/gauge.$i.info ./Configs/gauge.$new.info
done



