import fileinput
from math import *
N2 = []; N2_tot = 0.
N4 = []; N4_tot = 0.
filename = 'B.dat'
file = open(filename, "r")
for line in file:
    line = line.split()
    N2_i = float(line[0])
    N4_i = float(line[1])
    N2.append(N2_i) # put x2_i as the i-th element in an array x2
    N4.append(N4_i)
    N2_tot += N2_i
    N4_tot += N4_i
n = len(N2) # The number of lines in N2
# Do the jackknife estimates
N2_jack = []
N4_jack = []
for i in xrange(n):
    N2_jack.append((N2_tot - N2[i]) / (n - 1))
    N4_jack.append((N4_tot - N4[i]) / (n - 1))
    N2_av = N2_tot / n # do the overall averages
    N4_av = N4_tot / n
    R_av = ((N2_av / N4_av)* 0.4766 - 1)
    R_jack_av = 0.; R_jack_err = 0.
for i in xrange(n):
    dR = ((N2_jack[i] / N4_jack[i])* 0.4766 - 1)
    R_jack_av += dR
    #print 'Error',R_jack_av
    R_jack_err += dR**2
R_jack_av /= n
R_jack_err /= n
print 'R', R_jack_err
R_jack_err = sqrt((n - 1) * abs(R_jack_err - R_jack_av**2))
print " Overall average of Regge Curvature is %8.4f" % R_av
print " Jackknife average of Regge Curvature is %8.4f /pm %6.4f" % (R_jack_av, R_jack_err)
