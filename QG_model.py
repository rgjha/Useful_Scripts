import numpy as np
import math
import sys

N=10 
dim = int(pow(2, N/2)) 

print ("Dimension of H is", dim)
def areSame(A,B):

   for i in range(n):
      for j in range(n):
         if (A[i][j] != B[i][j]):
            return 0
   return 1

I32 = -1*np.identity(32)

sigma1 = [[0,1],[1,0]]
sigma2 = [[0,0-1.j],[0+1.j,0]]   
sigma3 = [[1,0],[0,-1]]
sigma4 = [[1,0],[0,1]]

gamma = np.zeros((N, dim, dim)) 
#Based on value of N, it constructs N Majorana fermions of size 2^N/2 x 2^N/2 

for i in range (0,N):


   if i ==0:
   	temp  = np.kron(sigma1, sigma1) 

   elif i < int(N/2):
   	temp = np.kron(temp, sigma1)

   #print ("Shape is ",  np.shape(gamma[i]))

#print ("Now check Clifford Algebra for our Gamma matrices")
#if (areSame(np.matmul(gamma1, gamma1), 1*I32)==1):
#   print("PASS")
#else:
#   print("ERROR")
