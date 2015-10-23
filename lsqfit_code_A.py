''' This fits the defined function to the data. One needs X (independent axis), Y (data) and
    the covariance matrix between measurements. An intial guess (known as 'prior') is necess-
    ary. The chi2/dof is printed in the output which establishes the trustability of the code
    along with the Q value. You need to install lsqfit from GitHub, and also gvar, cython etc. 
'''


import numpy as np
import scipy as sp
from numpy import ndarray
import os
import fileinput
import glob
import random
import time
import math
import collections
from scipy import sparse
from scipy.sparse import coo_matrix
from scipy import linalg
from numpy import matrix, rank
import scipy.sparse.linalg as spsl
import random
import gvar as gv
import lsqfit

''' Define variables/constants '''


list = []
origin = []
r2 = []
r3 = []
r4 = []
r5 = []
r6 = []
r7 = []
r8 = []
r9 = []
r10 = []
r11 = []
r12 = []
r13 = []
r14 = []
r15 = []
r16 = []

    
        
'''with open("Error_A1.dat") as f1:
    origin = f1.readlines()
    origin = [float(x.strip('\n')) for x in origin]
with open("Error_A2.dat") as f2:
    r2 = f2.readlines()
    r2 = [float(x.strip('\n')) for x in r2]'''



with open("Error_A3.dat") as f3:
    r3 = f3.readlines()
    r3 = [float(x.strip('\n')) for x in r3]
with open("Error_A4.dat") as f4:
    r4 = f4.readlines()
    r4 = [float(x.strip('\n')) for x in r4]
with open("Error_A5.dat") as f5:
    r5 = f5.readlines()
    r5 = [float(x.strip('\n')) for x in r5]
with open("Error_A6.dat") as f6:
    r6 = f6.readlines()
    r6 = [float(x.strip('\n')) for x in r6]
with open("Error_A7.dat") as f7:
    r7 = f7.readlines()
    r7 = [float(x.strip('\n')) for x in r7]
with open("Error_A8.dat") as f8:
    r8 = f8.readlines()
    r8 = [float(x.strip('\n')) for x in r8]
with open("Error_A9.dat") as f9:
    r9 = f9.readlines()
    r9 = [float(x.strip('\n')) for x in r9]
with open("Error_A10.dat") as f10:
    r10 = f10.readlines()
    r10 = [float(x.strip('\n')) for x in r10]
with open("Error_A11.dat") as f11:
    r11 = f11.readlines()
    r11 = [float(x.strip('\n')) for x in r11]
with open("Error_A12.dat") as f12:
    r12 = f12.readlines()
    r12 = [float(x.strip('\n')) for x in r12]
with open("Error_A13.dat") as f13:
    r13 = f13.readlines()
    r13 = [float(x.strip('\n')) for x in r13]
with open("Error_A14.dat") as f14:
    r14 = f14.readlines()
    r14 = [float(x.strip('\n')) for x in r14]
with open("Error_A15.dat") as f15:
    r15 = f15.readlines()
    r15 = [float(x.strip('\n')) for x in r15]
with open("Error_A16.dat") as f16:
    r16 = f16.readlines()
    r16 = [float(x.strip('\n')) for x in r16]






#print row
#print col
#print val

#AA = zip(row,col,val)
#B = np.row_stack((origin,r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12 ))
#B = np.row_stack((r3, r4, r5, r6, r7, r8, r9, r10, r11, r12 ))
#B = np.row_stack((r3, r4, r5, r6, r7, r8))
B = np.row_stack((r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16))
#B = np.row_stack((origin,r2, r3, r4, r5, r6))
#print B
print B.shape


#g = gv.correlate(g,B)
#print (gv.evalcorr(g))


covariance = np.cov(B)      # Covariance matrix of B #
print covariance


correlation = np.corrcoef(B)
#print correlation          # Correlation matrix #
print correlation.shape

#Ybar = [0 for x in range(1,4)]
#Ybar = [6.97, 6.38, 4.68]
with open("meanA.dat") as f:
    mean = f.readlines()
    mean = [float(x.strip('\n')) for x in mean]

#print mean


Y = gv.gvar(mean, covariance)
#print Y

#X = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.])
#X = np.array([3., 4., 5., 6., 7., 8., 9., 10., 11., 12.])
X = np.array([3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.,16.])
#X = np.array([3., 4., 5., 6., 7., 8.])
#X = np.array([1., 2., 3., 4., 5., 6.])


''' Define the function required for fit '''

def f(x, p):
    a = p['a']   # ****** #
    b = p['b']
    E = p['E']
    
    return a * x**(-b) * np.exp(-E * x)

def make_prior():
    prior = {}
    prior['a'] = gv.gvar(1.0, 1.0)
    prior['E'] = gv.gvar(4.5, 100.0)
    prior['b'] = gv.gvar(-10.0, 50.0)
    
    return prior


if __name__ == '__main__':
    gv.ranseed([2009, 2010, 2011, 2012])
    p0 = None
    prior = make_prior()
    #p0 = {'a':0.1,'E':0.5}   # ***** #
    #print f(X,p0)
    fit = lsqfit.nonlinear_fit(data = (X, Y), fcn=f, prior=prior, p0=p0)
    print fit
    print(fit.format(maxline=True))











            
            
