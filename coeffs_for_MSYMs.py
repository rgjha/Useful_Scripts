# As a check, this code should print 7.41 for energy coefficient p=0 case of QM #
# arXiv : 0706.0188

# arXiv : 0706.3517 for phase structure of QM 
# For high-T case in 1-d MSYM, refer to arXiv : 0710.2188
#http://arxiv.org/abs/0710.2188

# For high-T case in 2-d MSYM, refer to PRD 60, 106010 (below eq. 2.12)
# http://journals.aps.org/prd/pdf/10.1103/PhysRevD.60.106010

# Maldacena loop was used for predicting the Schwarzchild radius.
# This code prints 1.89 as the coefficient when p=0.
# https://arxiv.org/abs/0811.2081

# Tested the coefficients for p=2 as well from PRD 58, 046004# 
from math import *
import numpy as np
from scipy import *
from scipy.special import gamma

p = float(input('p : '))

dp = pow(2,(7-2*p))*pow(pi, (9-3*p)/(2))*gamma((7-p)/2)
dummy = (16*pi*pi*dp)/((7-p)*(7-p)) 
G = gamma((9-p)/2)
PI = pow(pi,(13-3*p)/2)
Pi_1=pow(pi,(18))
Pi_3=pow(pi,(22 - (4*p)))
Pi_2=pow(pi,(9-p))
OMEGA=2*pow(pi,(9-p)/2)/(G)
TWO = pow(2,(11-2*p))
NEW1=pow(2,(27-p))*Pi_1*(1/(7-p)**6)*(1/(OMEGA**2))
dummy1 = (9-p)/(G*PI*TWO)
dummy2 = (5-p)/(G*PI*TWO)
U = pow(dummy,((7-p)/(5-p)))
coeff= dummy1* U
EbyS=(9-p)/(7-p)*0.5 # Energy/Entropy with a factor of T ofcourse..
Free = -((coeff/EbyS) - coeff)
Free1= -dummy2*U  # Alternative
pr=((5-p)/(9-p))*coeff
S=(1/EbyS)*coeff 

ML = pow((7-p)/(4*pi*sqrt(dp)),(-1 + (3-p)/(5-p)))
ML = ML/(2*pi)


if p < 3:
    FACTOR=(1/sqrt(p+1))*pow((p+2), (4-p)/(6-2*p))
    print ("rt-lat/rt-cont = %8.4f" %(FACTOR))

if p == 3:
   FACTOR=sqrt(p+2)
   print ("lambda-lat/lambda-cont = %8.4f " %(FACTOR))

print ("Leading coeffecient of the energy is %8.4f" % (coeff))
print ("Leading coeffecient of the entropy is %8.4f" % (S))
print ("Leading coeffecient of the Free energy is %8.4f" % (Free))
print ("Leading coeffecient of the pressure (p) is %8.4f" % (pr))
# 2^(5/3)  Sqrt[Pi]  Gamma[0.8]^(5/3)/(3^(1/3)  Gamma[0.3]^(5/3))
# 4  Pi^2  Sqrt[2] / (Gamma[0.25]^4)
#print "Free energy directly evaluated is %8.4f" % (Free1)
#print "ALTERNATE : (E - ",p,"p) is %8.4f" % (coeff - p*pr)
print ("(E - ",p,"p) is %8.4f" % (coeff + p*Free))
print ("Coeffecient of the Maldacena-Wilson Loop is %8.4f" % (ML))
