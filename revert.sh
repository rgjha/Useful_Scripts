# Make sure we're calling in from right place and with correct number of args.

if [ $# -lt 2 ]; then
echo "Usage: $0 <filename>"
exit
fi

start=$1
end=$2

if [ -d Out ] ; then
echo "Okay" > /dev/null 2>&1
else
echo "Directory not found"
exit
fi

# Check whether this output file was run with usual definition or convoluted(my) definition
# Make changes if grep returns success

for(( i=$start ; $i<end ; i+=10 )); do
next=$[$i + 10]
out=./Out/out.$i-$next

kappa=$(awk '{if(/kappa=/) print $3}' < $out | cut -d')' -f 2 | cut -d'=' -f 2)
N=$(head -n4 $out | grep "Nc" | awk '{print $5}' | sed 's/,/ /g')
lambda_lat=$(bc <<< "scale=5;($N/(2.0*$kappa))")

if [  -n "`grep \"slnc_code_2\" $out`"  and lambda_lat=1/rt=t ]; then
    string1=$(grep -m 1 "lambda"  $out)
    sed -i "/$string1/c\lambda  $lambda_lat" $out (# Variable file name pass to sed ?)
    sed -i "/lambda=/c\lambda=$lambda_lat --> kappa=Nc/(2lambda)=$kappa" # ?
else
   echo "Nothing to change"
exit 1
fi
done
