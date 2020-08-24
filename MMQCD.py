# First pass: August 21, 2020
# Construction of bosonic basis states for variational Rayleigh-Ritz variational method 

import numpy as np
import sys
import collections
import itertools 
from sympy import LeviCivita
from numpy.random import permutation

Nc = 3     # Color group 
dimG = (Nc**2) - 1
size_col = dimG * 3 
n_O1_2boson = 24
n_O2_3boson = 324
n_04 = 3420
O1state = np.zeros((n_O1_2boson,3, dimG))
O1norm = np.zeros((n_O1_2boson))
O2state = np.zeros((n_O2_3boson, 3, dimG))
O2norm = np.zeros((n_O2_3boson))
O32state = np.zeros((n_O2_3boson, 3, dimG))
O32norm = np.zeros(n_O2_3boson)
O4state = np.zeros((n_04, 3, dimG))
O4norm = np.zeros((n_04))
O5state = np.zeros((n_04, 3, dimG))
O5norm = np.zeros((n_04))
matrix_bos = np.zeros((3, dimG))


# Setting up SU(3) structure constants 
# Note that standard f_123 is denoted as f_012 here.
def fabc(w):

  if len(w) != 3:
    print ("Needs three indices for computing f_abc, exiting!")
    sys.exit(1)

  out = 0.0
  if w in list(itertools.permutations((0,1,2))): 
    out = 1.0
  if w in list(itertools.permutations((0,3,6))) or w in list(itertools.permutations((1,3,5))): 
    out = 0.5
  if w in list(itertools.permutations((1,4,6))) or w in list(itertools.permutations((2,3,4))): 
    out = 0.5 
  if w in list(itertools.permutations((0,4,5))) or w in list(itertools.permutations((2,5,6))): 
    out = -0.5
  if w in list(itertools.permutations((3,4,7))) or w in list(itertools.permutations((5,6,7))): 
    out = np.sqrt(3)/2.0

  return out 


def action_A_or_Adag(v, in_vec, in_coeff, choice):
# Defines the action of A_ij or A_ij^dagger 

  out_vec = in_vec

  if len(v) != 2: 
    print ("Operator should carry two indices")
    sys.exit(1)

  if choice == 1:  # Creation
    out_vec[v[0]][v[1]] = in_vec[v[0]][v[1]] + 1
    out_coeff = in_coeff * np.sqrt(out_vec[v[0]][v[1]]) 
  
  if choice == 0:  # Annihilation 
    tmp = np.sqrt((out_vec[v[0]][v[1]]))
    out_vec[v[0]][v[1]] = out_vec[v[0]][v[1]] - 1

    if tmp >= 0:
      out_coeff = in_coeff * tmp 
    else:
      print ("Imaginary coefficient! Exit")
      sys.exit(1)


  return out_vec, out_coeff


def OP1():  # Operator 1 as per Aug 17 notes
# O1 = A_ij^dagger.A_ij^dagger 
   
  st_num = 0 
  for i, a in itertools.product(range(3), range(dimG)): 

      vac = np.zeros((3, dimG)) # Vacuum state
      O1state[st_num], O1norm[st_num] = action_A_or_Adag((i,a), vac, 1.0, 1)
      O1state[st_num], O1norm[st_num] = action_A_or_Adag((i,a),O1state[st_num], O1norm[st_num], 1) 

      st_num += 1 

  return O1state, O1norm


def OP2():  # Operator 2 
# f_abc ε_ijk (A†_ia * A†_jb * A†_kc) 
  
  st_num = 0 
  for i, a in itertools.product(range(3), range(dimG)): 
    for j, b in itertools.product(range(3), range(dimG)):
      for k, c in itertools.product(range(3), range(dimG)):

        tmp = fabc((a,b,c)) * LeviCivita(i, j, k)
        if tmp != 0:

            # This happens for only 324 (=6 x 54) combinations of a,b,c,i,j,k
            # Out of possible 8^3 x 27 = 512 x 27 = 13824  combinations
          
          vac = np.zeros((3, dimG)) # Vacuum state
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((k,c), vac, 1.0, 1)
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((j,b),O2state[st_num], O2norm[st_num], 1)
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((i,a),O2state[st_num], O2norm[st_num], 1) 
          O2norm[st_num] *= tmp 
          st_num += 1

  return  O2state, O2norm


def OP32():  # Operator 3 on 2 on vac
# f_abc ε_ijk (A_ia * A†_jb * A†_kc)
# For ex: f_584 ε_321 (A_35 * A†_28 * A†_14)
# -sqrt(3)/2 (A_35 * A†_28 * A†_14) | O2state >
  num = 0

  for i, a in itertools.product(range(3), range(dimG)):
    for j, b in itertools.product(range(3), range(dimG)):
      for k, c in itertools.product(range(3), range(dimG)):

        tmp = fabc((a,b,c)) * LeviCivita(i, j, k)
        if tmp != 0:

          O32state[num] = O2state[num]
          O32norm[num] = O2norm[num]

          O32state[num], O32norm[num] = action_A_or_Adag((k,c),O32state[num], O32norm[num], 1)
          O32state[num], O32norm[num] = action_A_or_Adag((j,b),O32state[num], O32norm[num], 1)
          O32state[num], O32norm[num] = action_A_or_Adag((i,a),O32state[num], O32norm[num], 0)
          O32norm[num] *= tmp
          num += 1

  return  O32state, O32norm


def OP4():  # Operator 4 on vac
# f_abc f_ade (A†_ib * A†_jc * A†_id * A†_je) | vac. >

  st_num = 0

  for a in range(dimG):
    for b in range(dimG):
      for c in range(dimG):
        for d in range(dimG):
          for e in range(dimG):

            tmp = fabc((a,b,c)) * fabc((a,d,e))
            if tmp != 0:


              for i in range(3):
                for j in range(3):

                  vac = np.zeros((3, dimG)) # Vacuum state
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((j,e), vac, 1.0, 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((i,d),O4state[st_num], O4norm[st_num], 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((j,c),O4state[st_num], O4norm[st_num], 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((i,b),O4state[st_num], O4norm[st_num], 1)
                  O4norm[st_num] *= tmp
                  st_num += 1


  return  O4state, O4norm


def OP5():  # Operator 4 on vac
# f_abc f_ade (A†_ib * A†_jc * A†_id * A†_je) | vac. >

  st_num = 0

  for a in range(dimG):
    for b in range(dimG):
      for c in range(dimG):
        for d in range(dimG):
          for e in range(dimG):

            tmp = fabc((a,b,c)) * fabc((a,d,e))
            if tmp != 0:


              for i in range(3):
                for j in range(3):

                  vac = np.zeros((3, dimG)) # Vacuum state
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,e), vac, 1.0, 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,d),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,c),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,b),O5state[st_num], O5norm[st_num], 0)
                  O5norm[st_num] *= tmp
                  st_num += 1


  return  O5state, O5norm


if __name__ == "__main__":

  OP1()
  OP2()
  OP32()
  OP4() 
  OP5() 

  for i in range (10): # Print first 10 states of each type for now. 
    print ("O1:", O1state[i].reshape(size_col), "with coefficient", O1norm[i]) # 24 states #2-boson state 
    print ("O2:", O2state[i].reshape(size_col), "with coefficient", O2norm[i]) # 324 states  #3-boson state 
    print ("O32:",O32state[i].reshape(size_col), "with coefficient", O32norm[i]) # 324 states  #5-boson state 
    print ("O4:",O4state[i].reshape(size_col), "with coefficient", O4norm[i]) # 3420 states # 4-boson state 
    print ("O5:",O5state[i].reshape(size_col), "with coefficient", O5norm[i]) # 3420 states # 2-boson state 



