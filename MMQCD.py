# First pass: August 21, 2020
# Construction of bosonic basis states for variational Rayleigh-Ritz variational method 

import numpy as np
import sys
import collections
import itertools 
from sympy import LeviCivita
from numpy.random import permutation

Nc = 3     # Color group 
size_col = ((Nc**2)-1)*3 
O1state = np.zeros((3,(Nc**2)-1,3,(Nc**2)-1))
O1norm = np.zeros((3,(Nc**2)-1))
matrix_bos = np.zeros((3,(Nc**2)-1)) 


def action_A_or_Adag(v, in_vec, coeff, flag):
# Defines the action of A_ij or A_ij^dagger 

  if flag == 1:  # Creation
    in_vec[v[0]][v[1]] += 1
    coeff *= np.sqrt(in_vec[v[0]][v[1]]) 
  
  if flag == 0:  # Annihilation 
    coeff *= np.sqrt(in_vec[v[0]][v[1]])
    in_vec[v[0]][v[1]] -= 1
  
  return in_vec, coeff


def OP1():  # Operator 1 as per Aug 17 notes
# O1 = A_ij^dagger.A_ij^dagger 

  for i, a in itertools.product(range(3), range((Nc**2)-1)): 

      vac = np.zeros((3,(Nc**2)-1)) # Vacuum state
      O1state[i][a], O1norm[i][a] = action_A_or_Adag((i,a), vac, 1.0, 1)
      O1state[i][a], O1norm[i][a] = action_A_or_Adag((i,a),O1state[i][a], O1norm[i][a], 1) 

  return O1state, O1norm

# Setting up SU(3) structure constants 
def fabc(w):

  out = 0.0
  if w in list(itertools.permutations((1,2,3))): 
    out = 1.0
  if w in list(itertools.permutations((1,4,7))) or list(itertools.permutations((2,4,6))): 
    return 0.5 
  if w in list(itertools.permutations((2,5,7))) or w in list(itertools.permutations((3,4,5))): 
    return 0.5 
  if w in list(itertools.permutations((1,5,6))) or list(itertools.permutations((3,6,7))): 
    return -0.50
  if w in list(itertools.permutations((4,5,8))) or list(itertools.permutations((6,7,8))): 
    return np.sqrt(3)/2.0

  return out 


if __name__ == "__main__":

  OP1()
  #print(LeviCivita(3, 2, 1)) # Epsilon_abc 
  #print (fabc((4,3,5))) # f_abc 

  for i, a in itertools.product(range(3), range((Nc**2)-1)): 
    print (O1state[i][a])
  print ("Coefficient matrix \n", O1norm)
