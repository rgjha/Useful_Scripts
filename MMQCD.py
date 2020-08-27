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
n_h = 3420
O1state = np.zeros((n_O1_2boson,3, dimG))
O1norm = np.zeros((n_O1_2boson))
O2state = np.zeros((n_O2_3boson, 3, dimG))
O2norm = np.zeros((n_O2_3boson))
O3state = np.zeros((n_O2_3boson, 3, dimG))
O3norm = np.zeros(n_O2_3boson)
O4state = np.zeros((n_h, 3, dimG))
O4norm = np.zeros((n_h))
O5state = np.zeros((n_h, 3, dimG))
O5norm = np.zeros((n_h))
O6state = np.zeros((n_h, 3, dimG))
O6norm = np.zeros((n_h))
O6totalstate = np.zeros((n_h, 3, dimG))
O6totalnorm = np.zeros((n_h))
vac = np.zeros((3, dimG)) # Vacuum state

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
      print ("Imaginary coefficient!", tmp, " Exit")
      #sys.exit(1)
      return 0 


  return out_vec, out_coeff


def OP1(in_state, in_state_norm):  
# Operator 1 as per Aug 17 notes
# O1 = A_ij^dagger.A_ij^dagger | in state, in_state_norm >
   
  st_num = 0 

  for i, a in itertools.product(range(3), range(dimG)): 

    O1state[st_num] = in_state
    O1norm[st_num] = in_state_norm

    O1state[st_num], O1norm[st_num] = action_A_or_Adag((i,a),O1state[st_num], O1norm[st_num], 1)
    O1state[st_num], O1norm[st_num] = action_A_or_Adag((i,a),O1state[st_num], O1norm[st_num], 1) 

    st_num += 1 

  return O1state, O1norm


def OP2(in_state, in_state_norm): 
# O2 | . > = f_abc ε_ijk (A†_ia * A†_jb * A†_kc) | in state, in_state_norm >
# Note that out of possible 8^3 x 27 = 512 x 27 = 13824  combinations
# of a,b,c,i,j,k -- only 324 (=6 x 54) are non-zero!
  
  st_num = 0 
  for i, a in itertools.product(range(3), range(dimG)): 
    for j, b in itertools.product(range(3), range(dimG)):
      for k, c in itertools.product(range(3), range(dimG)):

        tmp = fabc((a,b,c)) * LeviCivita(i, j, k)
        if tmp != 0:
          
          O2state[st_num] = in_state
          O2norm[st_num] = in_state_norm
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((k,c), O2state[st_num], O2norm[st_num], 1)
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((j,b),O2state[st_num], O2norm[st_num], 1)
          O2state[st_num], O2norm[st_num] = action_A_or_Adag((i,a),O2state[st_num], O2norm[st_num], 1) 
          O2norm[st_num] *= tmp 
          st_num += 1

  return  O2state, O2norm


def OP3(in_state, in_state_norm):  
# f_abc ε_ijk (A_ia * A†_jb * A†_kc) | in state, in_state_norm >
# For ex: Operator 3 on 2 on vac
# f_abc ε_ijk (A_ia * A†_jb * A†_kc)
# f_584 ε_321 (A_35 * A†_28 * A†_14)
# -sqrt(3)/2 (A_35 * A†_28 * A†_14) | O2state >
  num = 0

  for i, a in itertools.product(range(3), range(dimG)):
    for j, b in itertools.product(range(3), range(dimG)):
      for k, c in itertools.product(range(3), range(dimG)):

        tmp = fabc((a,b,c)) * LeviCivita(i, j, k)
        if tmp != 0:

          O3state[num] = in_state[num]
          O3norm[num] = in_state_norm[num]

          O3state[num], O3norm[num] = action_A_or_Adag((k,c),O3state[num], O3norm[num], 1)
          O3state[num], O3norm[num] = action_A_or_Adag((j,b),O3state[num], O3norm[num], 1)
          O3state[num], O3norm[num] = action_A_or_Adag((i,a),O3state[num], O3norm[num], 0)
          O3norm[num] *= tmp
          num += 1

  return  O3state, O3norm


def OP4(in_state, in_state_norm): 
# f_abc f_ade (A†_ib * A†_jc * A†_id * A†_je) | in state, in_state_norm >

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

                  O4state[st_num] = in_state
                  O4norm[st_num] = in_state_norm

                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((j,e), O4state[st_num], O4norm[st_num], 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((i,d),O4state[st_num], O4norm[st_num], 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((j,c),O4state[st_num], O4norm[st_num], 1)
                  O4state[st_num], O4norm[st_num] = action_A_or_Adag((i,b),O4state[st_num], O4norm[st_num], 1)
                  O4norm[st_num] *= tmp
                  st_num += 1


  return  O4state, O4norm


def OP5(in_state, in_state_norm):
# f_abc f_ade (A_ib * A†_jc * A†_id * A†_je) | in state, in_state_norm >

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

                  O5state[st_num] = in_state
                  O5norm[st_num] = in_state_norm

                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,e),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,d),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,c),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,b),O5state[st_num], O5norm[st_num], 0)
                  O5norm[st_num] *= tmp
                  st_num += 1


  return  O5state, O5norm


def OP6(in_state, in_state_norm):  
# f_abc f_ade [(A_ib * A_jc * A†_id * A†_je) + (A_ib * A_id * A†_jc * A†_je) +  \ 
# (A_ib * A_je * A†_id * A†_jc)] | in state, in_state_norm >

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


                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]

                  print ("with coeff . ", O6norm[st_num])
                  
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c), O6state[st_num], O6norm[st_num], 1)
                  print ("with coeff .. ", O6norm[st_num])
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 1)
                  print ("with coeff ... ", O6norm[st_num])
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e),O6state[st_num], O6norm[st_num], 0)
                  print ("with coeff .... ", O6norm[st_num])
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)
                  O6totalstate[st_num] = O6state[st_num]
                  O6totalnorm[st_num] = O6norm[st_num]

                  print (O6state[st_num].reshape(24), "with coeff", O6norm[st_num])


                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]

                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e), O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c),O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 0)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)

                  O6totalstate[st_num] += O6state[st_num]
                  O6totalnorm[st_num] += O6norm[st_num]

                  print (O6state[st_num].reshape(24), "with coeff", O6norm[st_num])

                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]

                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e), O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c),O6state[st_num], O6norm[st_num], 0)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)

                  print (O6state[st_num].reshape(24), "with coeff", O6norm[st_num])

                  O6totalstate[st_num] += O6state[st_num]
                  O6totalnorm[st_num] += O6norm[st_num]

                  O6totalnorm[st_num] *= tmp
                  st_num += 1

                  sys.exit(1)


  return  O6totalstate, O6totalnorm


if __name__ == "__main__":


  out_state1, norm_out_state1 = OP1(vac, 1.0) 
  out_state2, norm_out_state2 = OP2(vac, 1.0)
  out_state32, norm_out_state32 = OP3(O2state, O2norm)
  out_state4, norm_out_state4 = OP4(vac, 1.0) 
  out_state5, norm_out_state5 = OP5(vac, 1.0) 
  out_state64, norm_out_state64 = OP6(O4state, O4norm) 

  for i in range (10): # Print first 10 states of each type for now. 
    print ("O1 |vac> :", out_state1[i].reshape(size_col), "with coefficient", round(norm_out_state1[i],3)) # 24 states #2-boson state 
    print ("O2 |vac> :", out_state2[i].reshape(size_col), "with coefficient", round(norm_out_state2[i],3)) # 324 states  #3-boson state 
    print ("O3 O2 |vac> :",out_state32[i].reshape(size_col), "with coefficient", round(norm_out_state32[i],3)) # 324 states  #5-boson state 
    print ("O4 |vac> :",out_state4[i].reshape(size_col), "with coefficient", round(norm_out_state4[i],3)) # 3420 states # 4-boson state 
    print ("O5 |vac> :",out_state5[i].reshape(size_col), "with coefficient", round(norm_out_state5[i],3)) # 3420 states # 2-boson state 
    print ("O6 O4 |vac> :",out_state64[i].reshape(size_col), "with coefficient", round(norm_out_state64[i],3)) # 3420 states



