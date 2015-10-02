# -*- coding: utf-8 -*-

# This code fits the data to the function defined below #
# October 1, 2015 #
# =============================================================

from pylab import *
from scipy.optimize import curve_fit
from scipy.stats import chi2



fname = sys.argv[1]
x, y, err = np.loadtxt(fname, unpack=True) # Read the data in #
n = len(x) # No. of data points #
p0 = [20, 0, 0] # Initial values of parameters #


f = lambda x, a, b, c: a * x**-b * np.exp(-c * x)  # Fitting function #
#f = lambda x, a, b, c: a*x**2 + b*x + c


# Read about lambda notation here : http://www.secnetix.de/olli/Python/lambda_functions.hawk #


p, covm = curve_fit(f, x, y, p0, err) # Do the fit
a, b, c = p

chisq = sum(((f(x, a, b, c) - y)/err)**2) # Compute chi-squared
ndf = n -len(p) # no. of degrees of freedom
print "NDF is " , ndf
Q = 1. - chi2.cdf(chisq, ndf) # Quality of fit parameter : Q , More the better ! #
chisq = chisq / ndf # Compute chi-squared per DOF

aerr, berr, cerr = sqrt(diag(covm)/chisq) # Correct the error bars


# ======================================
# Print the results #

print "a = %10.4f +/- %7.4f" % (a, aerr)
print "b = %10.4f +/- %7.4f" % (b, berr)
print "c = %10.4f +/- %7.4f" % (c, cerr)
print "chi squared / NDF = %7.4lf" % chisq
print "Q = %10.4f" % Q

# Chi-square calculator : https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html # 



