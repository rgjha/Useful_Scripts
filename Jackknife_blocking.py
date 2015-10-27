'''This code calculates the jackknife error
   First check/edit : September 12, 2015
   Implemented blocking on October 27, 2015 '''

import sys
from math import *
data = []; data_tot = 0. ; Data = [] ; data_jack = []
filename = sys.argv[1]      # Run this as : python knife.py filename.dat

file = open(filename, "r")
for line in file:
    line = line.split()
    data_i = float(line[0])
    data.append(data_i)
    data_tot += data_i
n = len(data)               # The number of lines of data
print "===THIS IS BLOCKED JACKKNIFE CODE==="
print "Data Size is :: ", n
nb = input('Blocking/Binning Size : ')
blocksize = int(nb)

n_b = n / blocksize
B = 0.

for k in range(n_b):
    for w in range((k*blocksize)+1,(k*blocksize)+blocksize+1):
        B += data[w-1]
    
    Data.insert(k,B)
    B = 0

'''print "Data[0] is" , Data[0]
print "Data[1] is" , Data[1]
print "Data[2] is" , Data[2]
print "Data[3] is" , Data[3]
print "Data[4] is" , Data[4]
print "Data TOT", data_tot'''

''' Do the jackknife estimates '''


for i in xrange(n_b-1):
    data_jack.append((data_tot - Data[i]) / (n - blocksize))
    data_av = data_tot / n   # Do the overall averages
    data_av = data_av
    data_jack_av = 0.; data_jack_err = 0.
for i in xrange(n_b-1):
    dR = data_jack[i]
    data_jack_av += dR
    data_jack_err += dR**2

data_jack_av /= n_b-1
data_jack_err /= n_b-1

data_jack_err = sqrt((n_b - 2) * abs(data_jack_err - data_jack_av**2))
print " JK Average : %.2E "  " JK Error : %.2E" % (data_jack_av, data_jack_err)



