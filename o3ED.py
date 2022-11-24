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

from sympy import Symbol, Integer
from sympy.physics.quantum import (Dagger, qapply,represent,InnerProduct,Commutator)
from sympy.physics.quantum.sho1d import (RaisingOp, LoweringOp,NumberOp,Hamiltonian,SHOKet,SHOBra)
ad = RaisingOp('a')
a = LoweringOp('a')
NOP = NumberOp('N')
H = Hamiltonian('H')

N = 4 
SQRT2 = 1.4142135623730951
SQRT3 = 1.7320508075688772

if len(sys.argv) < 3:
  print("Usage: python", str(sys.argv[0]), "l-max selected " "beta=1/g^2")
  sys.exit(1)

lmax = float(sys.argv[1])
beta = float(sys.argv[2])
sigmax = np.array([[0, 1],[1, 0]])
isigmay =  np.array([[0, 1],[-1, 0]])
sigmaz =  np.array([[1, 0],[0, -1]])

def dagger(a):
    return np.transpose(a).conj()

nplus = -(sigmax + isigmay)/(3*SQRT2)
nminus = -(sigmax - isigmay)/(3*SQRT2)
nz = sigmaz 

if lmax == 0.5:

   dimH = 2
   I2 = np.identity(2)
   I4 = np.identity(4)
   Ikin = np.identity(16)
   # beta = 10 gives -2.18472 for N=4

if lmax == 1.5:

   dimH = 6
   Ikin = np.identity(dimH**4)
   I2 = np.identity(dimH)
   I4 = np.identity(dimH**2)
   # How to pad: np.pad(nplus, (0,padup), 'constant', constant_values=(0)) 

   nplus =  np.array([[0, 0, 1./SQRT3, 0, 0, 0],[-SQRT2/3., 0, 0, 1/3., 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, -SQRT2/(np.sqrt(3)*5.), 0, 0, 0], [-1/3., 0, 0, -SQRT2*2/15., 0, 0],[0, -1./SQRT3, 0, 0, -SQRT2/(np.sqrt(3)*5.), 0]])
   nminus = np.array([[0, -SQRT2/3., 0, 0, -1/3., 0],[0, 0, 0, 0, 0, -1/SQRT3],[1./SQRT3, 0, 0, -SQRT2/(np.sqrt(3)*5.), 0, 0],[0, 1/3., 0, 0, -SQRT2*2/15., 0], [0, 0, 0, 0, 0, -SQRT2/(np.sqrt(3)*5.)],[0, 0, 0, 0, 0, 0]])
   nz = np.array([[1, 0, 0, SQRT2, 0, 0],[0, -1, 0, 0, SQRT2, 0],[0, 0, 0.6, 0, 0, 0],[SQRT2, 0, 0, 0.2, 0, 0], [0, SQRT2, 0, 0, -0.2, 0],[0, 0, 0, 0, 0, -0.6]])
   # beta = 10 gives -5.668232717889893 for N=4

HH = (N*(3./8.)*(1/beta)*Ikin) + \
beta*(np.kron(np.kron(nplus, nminus), I4) + np.kron(np.kron(nminus, nplus), I4) + (1./9.0)*np.kron(np.kron(nz, nz), I4)) + \
beta*(np.kron(np.kron(np.kron(I2, nplus), nminus),I2) + np.kron(np.kron(np.kron(I2, nminus), nplus),I2) + (1./9.0)*np.kron(np.kron(np.kron(I2, nz), nz), I2)) + \
beta*(np.kron(np.kron(I4, nplus),nminus) + np.kron(np.kron(I4, nminus),nplus) + (1./9.0)*np.kron(np.kron(I4, nz), nz)) + \
beta*(np.kron(np.kron(nminus, I4),nplus) + np.kron(np.kron(nplus, I4),nminus) + (1./9.0)*np.kron(np.kron(nz, I4), nz))

# Check hermitian nature of H
print ("Is the Hamiltonian Hermitian?", np.allclose(HH, dagger(HH)))

#with np.printoptions(precision=10, suppress=True, formatter={'float': '{:0.2f}'.format}, linewidth=100):
#    print(HH)

print ("Size of the Hamiltonian", np.shape(HH))
w, v = eig(HH)
print ("Ground state energy from ED of Pauli Hamiltonian", np.min(w.real)/N)

sys.exit(1) 
# TODO from here!

# Now think about H in terms of creation/annihilation operators(qumodes-H)
# Ideally, the ground state energy should match with both qubit-H
# and qumodes-H. We will try this! 

# First try the kinetic term. 
number_operator = np.zeros((size, size), int)

for i in range (0, size):
   number_operator[i][i] = (i+1)*2.

HH = (1/2)*(1/beta)*(0.25*number_operator*(0.25*number_operator+np.eye(16)))
w, v = eig(HH)
print ("Ground state energy using qumodes H (kinetic piece)", np.min(w))


adag= np.zeros((size, size), float)
a = np.zeros((size, size), float)

bdag= np.zeros((size, size), float)
b = np.zeros((size, size), float)

for i in range (1, size):
   adag[i][i-1] = math.sqrt(i)
   a[i-1][i] = math.sqrt(i)
   bdag[i][i-1] = math.sqrt(i)
   b[i-1][i] = math.sqrt(i)

NN = adag @ a + bdag @ b

HH = (3/2)*(1/beta)*(1.0*(NN + 0.5*np.eye(16)@(1.0*(NN + 0.5*np.eye(16)))))
with np.printoptions(precision=10, suppress=True, formatter={'float': '{:0.2f}'.format}, linewidth=100):
   print (NN)
w, v = eig(HH)
print ("Ground state energy using qumodes H (kinetic piece)", np.min(w))
