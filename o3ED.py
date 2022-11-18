# Work in progress. 
# November 2022

import numpy as np
import math
import sys
from numpy import linalg as LA
from numpy.linalg import eig
from numpy.linalg import svd
from sympy import Matrix, init_printing
init_printing()

N=4 
SQRT2 = 1.4142135623730951
dimH=6
padup=int(dimH-2)
size=2**N

if len(sys.argv) < 3:
  print("Usage: python", str(sys.argv[0]), "l-max selected " "beta=1/g^2")
  sys.exit(1)

lmax = float(sys.argv[1])
beta = float(sys.argv[2])
sigmax = np.array([[0, 1],[1, 0]])
isigmay =  np.array([[0, 1],[-1, 0]])
sigmaz =  np.array([[1, 0],[0, -1]])

nplus = -(sigmax + isigmay)/(3*SQRT2)
nminus = -(sigmax - isigmay)/(3*SQRT2)

if lmax == 0.5:

   I2 = np.identity(2)
   I4 = np.identity(4)
   Ikin = np.identity(16)
   # beta=10 gives -2.18472

if lmax == 1.5:

   Ikin = np.identity(dimH**4)
   sigmaz = np.pad(sigmaz, (0,padup), 'constant', constant_values=(0)) 
   isigmay = np.pad(isigmay, (0,padup), 'constant', constant_values=(0)) 
   sigmax = np.pad(sigmax, (0,padup), 'constant', constant_values=(0)) 
   I2 = np.identity(dimH)
   I4 = np.identity(dimH**2)
   #nplus = np.kron(np.kron(nplus, nplus), nplus) 
   #nminus = np.kron(np.kron(nminus, nminus), nminus) 
   nplus = np.pad(nplus, (0,padup), 'constant', constant_values=(0)) 
   nminus = np.pad(nminus, (0,padup), 'constant', constant_values=(0)) 


HH = (N*(3./8.)*(1/beta)*Ikin) + \
beta*(np.kron(np.kron(nplus, nminus), I4) + np.kron(np.kron(nminus, nplus), I4) + (1./9.0)*np.kron(np.kron(sigmaz, sigmaz), I4)) + \
beta*(np.kron(np.kron(np.kron(I2, nplus), nminus),I2) + np.kron(np.kron(np.kron(I2, nminus), nplus),I2) + (1./9.0)*np.kron(np.kron(np.kron(I2, sigmaz), sigmaz), I2)) + \
beta*(np.kron(np.kron(I4, nplus),nminus) + np.kron(np.kron(I4, nminus),nplus) + (1./9.0)*np.kron(np.kron(I4, sigmaz), sigmaz)) + \
beta*(np.kron(np.kron(nminus, I4),nplus) + np.kron(np.kron(nplus, I4),nminus) + (1./9.0)*np.kron(np.kron(sigmaz, I4), sigmaz)) 


# Try padding H rather than what makes H 
#HH1= np.pad(HH, (0,4080), 'constant', constant_values=(0))

#with np.printoptions(precision=10, suppress=True, formatter={'float': '{:0.2f}'.format}, linewidth=100):
#    print(HH)

print ("Size of the Hamiltonian", np.shape(HH))
w, v = eig(HH)
print ("Ground state energy from ED of Pauli Hamiltonian", np.min(w.real)/N)


# Now think about H in terms of creation/annihilation operators(qumodes-H)
# Ideally, the ground state energy should match with both qubit-H
# and qumodes-H. We will try this! 

# First try the kinetic term. 
number_operator = np.zeros((size, size), int)

for i in range (0, size):
   number_operator[i][i] = (i+1)*2.

HH = (1/2)*(1/beta)*(0.25*number_operator*(0.25*number_operator+1))
w, v = eig(HH)
print ("Ground state energy using qumodes H (kinetic piece)", np.min(w))
