#!/bin/bash

# Find file with some name 'abc'
find . -type f -name "abc*"

# Change space separated to comma separated
sed 's/ /,/g' filename

# Change tab separated to comma separated
sed 's/\t/,/g' filename

# Change comma separated by tab separated and add each to new line 
# Useful for Mathematica sometimes
cat test.txt | tr '[,]' '[\t]' | tr " " "\n"
or, 
cat matrix.txt | tr '[,]' '[\t]' > mat2.txt

# Print specific column based on matching other columns, print appropriate decimals, sort the output by 1st colum ascending
awk '{ if($4==45 && $5==50 && $7==0) print $1, $3}' FILENAME | awk '{printf "%.6f %.4f\n", $1, $2}' |  sort -k1 -n

# Round off a text file 
awk '{ printf("%.3g %.3g\n", $1, $2) }' file

# Merge line 'n' and 'n+1' together on one line 
awk 'NR%2{printf "%s ",$0;next;}1' FILENAME

# Check the magnetization for my tensor script 
awk 'NR%2{printf "%s ",$0;next;}1' input >> tmp1
awk '{print $1, ($9-$3)/($12-$6), 0.50*($12+$6),$4, $5}' tmp1

# Delete files below some size 2700 (usually in units of ls -lrt output)
find . -type f -size -2700c -exec rm '{}' \;

# Find files whose size is less than 100 bytes
find . -type f -size -100c

# Untar ABC.tar.gz
gunzip -c ABC.tar.gz | tar xopf -

# Line with no keyword
grep -w -v -e "done"  FILE

# Restrict decimal of 'awk' output
Some output | awk '{printf "%.2f %.2f\n", $1, $2}'

# Check internet speed etc. (also prints IP address)
curl -s  https://raw.githubusercontent.com/sivel/speedtest-cli/master/speedtest.py | python -

# Print 1 to N (integer step) in bash
printf "%d\n" {1..N} >> MDTU

# Cut Line 10 to 100 from some file
sed -n -e '10,100p' input.txt > output.txt

# Merge many eps files into one pdf #
convert -density 600 *.eps -resize 100% new.pdf

# Sort file in ascending order
sort -n $file

# Cut first 50 lines of a file
tail -n +50  infile > outfile

# RMS of a file
awk '{ sum += $1*$1; n++ } END { if (n > 0) print sqrt(sum / n); }' FILE

# Average of a column in a file
awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' FILE

# Maximum entry of a file
cat FILE | awk '{if ($1 > max) max=$1}END{print max}'

# Minimum entry
cat FILE  | awk 'BEGIN{min=512} {if ($1 < min) min=$1}END{print min}'

# If there is empty/blank lines in a file, it will remove them #
sed '/^\s*$/d' OLDFILE > NEWFILE

# Prints number of files with specific name in a given folder #
exts=( out.*); printf "There are ${#exts[@]} of them " extensions;

# Number of lines in a folder
find . -name '*.h' | xargs wc -l

# Merge files with one column each side by side #
pr -m -t OLDFILE1 OLDFILE2 > NEWFILE

# Alternative to pr -m -t above 
paste N* | column -s $'\t' -t

# Takes all the files with specific name in a folder and averages line by line #
awk '{a[FNR]+=$0;b[FNR]++;}END{for(i=1;i<=FNR;i++)print a[i]/b[i];}' FILES >> FILE

# Run on Condor without interruption #
nohup time bash ./script.sh > OUTFILE 2>&1

# Kill processes
lsof +D .
kill -9 PID#

# Delete files whose name is mentioned in FILE
xargs rm < FILE

# Sort a file
sort FILE -o FILE

# Run through files named *.txt in ascending order
for d in `ls *.txt | sort -V`; do

# Split SB.csv comma separated to file
awk -F',' '{print $2}' SB.csv >> sb.txt

# Cut a file between two line numbers passed
sed -n "${args[0]},${args[1]} p" FINAL_4K_m0125 > 4K_0125

# Change all extensions in directory
for d in *.c ; do
filename="${d%.*}"
mv $filename.c $filename.cpp
done

# Change SSH key to not being asked to enter password each time
ssh-keygen
ssh-copy-id -i ~/.ssh/id_rsa.pub XXX@comet.sdsc.edu # Replace appropraitely
ssh XXX@comet.sdsc.edu # Check that all is fine

# Python 2 to 3
for d in *.py ; do
2to3 -w $d  >> /dev/null 2>&1
sed 's/python/python3/g' $d >> tmp
mv tmp $d
done

# Submit jobs on Symmetry 
#!/bin/bash 

if [[ $# -ne 2 ]] ; then
    echo '*** Forgetting arguments !! Order : "TEMP"  "FIELD" ' 
    exit 0
fi

args=("$@")

echo "#!/bin/bash" >> sub.sh
echo "#SBATCH --job-name=SMA" >> sub.sh
echo "#SBATCH --nodes=1" >> sub.sh
echo "#SBATCH --output=/home/rjha1/log${args[0]}_${args[1]}.txt" >> sub.sh
echo "#SBATCH --time=10:00:00" >> sub.sh
echo "" >> sub.sh

echo "set -euxo pipefail" >> sub.sh
echo "" >> sub.sh
echo "pwd" >> sub.sh
echo "echo \"SLURM_JOB_ID=\$SLURM_JOB_ID\"" >> sub.sh
echo "env OMP_NUM_THREADS=4 python 2d_XYv2.py $1 $2" >> sub.sh
sbatch sub.sh
rm sub.sh


# Move back configuration & files etc. 

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








