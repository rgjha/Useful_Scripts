# First pass: August 21, 2020
# Construction of bosonic basis states for variational Rayleigh-Ritz variational method 

import numpy as np
import sys
import collections
import itertools 
from numpy import linalg as LA
from sympy import LeviCivita
from numpy.random import permutation
#from opt_einsum import contract
from ncon import ncon 


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
Lambda = np.zeros((dimG+1,3,3),dtype=complex) # Gell-Mann matrices 
Sigma = np.zeros((4,2,2),dtype=complex) # Pauli matrices 
Smatrix = np.zeros((6,6),dtype=complex)


def list_even(a,b,c):
  return (list(itertools.permutations((a,b,c)))[0], list(itertools.permutations((a,b,c)))[3], list(itertools.permutations((a,b,c)))[4]) 
def list_odd(a,b,c):
  return (list(itertools.permutations((a,b,c)))[1], list(itertools.permutations((a,b,c)))[2], list(itertools.permutations((a,b,c)))[5]) 
# Setting up SU(3) structure constants 
# Note that standard f_123 is denoted as f_012 here.
# Similarly, for d_abc 


def fSU3(w):

  if len(w) != 3:
    print ("Needs three indices for computing SU(3) structure constant, exiting!")
    sys.exit(1)

  out = 0.0
  if w in list_even(0,1,2): 
    out = 1.0
  if w in list_odd(0,1,2):
    out = -1.0
  if w in list_even(0,3,6) or w in list_even(1,3,5): 
    out = 0.5
  if w in list_odd(0,3,6)  or w in list_odd(1,3,5): 
    out = -0.5
  if w in list_even(1,4,6) or w in list_even(2,3,4): 
    out = 0.5 
  if w in list_odd(1,4,6)  or w in list_odd(2,3,4): 
    out = -0.5
  if w in list_even(0,4,5) or w in list_even(2,5,6): 
    out = -0.5 
  if w in list_odd(0,4,5)  or w in list_odd(2,5,6): 
    out = 0.5
  if w in list_even(3,4,7) or w in list_even(5,6,7): 
    out = np.sqrt(3)/2.0 
  if w in list_odd(3,4,7)  or w in list_odd(5,6,7): 
    out = -np.sqrt(3)/2.0 

  return out 


def dSU3(w):

  if len(w) != 3:
    print ("Needs three indices for computing SU(3) structure constant, exiting!")
    sys.exit(1)

  out = 0.0
  if w in list(itertools.permutations((0,0,7))) or w in list(itertools.permutations((1,1,7))) or w in list(itertools.permutations((2,2,7))):
    out = 1.0/np.sqrt(3) 
  if w in list(itertools.permutations((7,7,7))):
    out = -1.0/np.sqrt(3) 
  if w in list(itertools.permutations((3,3,7))) or w in list(itertools.permutations((4,4,7))):
    out = -0.50/np.sqrt(3)
  if w in list(itertools.permutations((5,5,7))) or w in list(itertools.permutations((6,6,7))): 
    out = -0.50/np.sqrt(3)
  if w in list(itertools.permutations((0,3,5))) or w in list(itertools.permutations((0,4,6))) or w in list(itertools.permutations((1,4,5))): 
    out = 0.50
  if w in list(itertools.permutations((2,3,3))) or w in list(itertools.permutations((2,4,4))): 
    out = 0.50
  if w in list(itertools.permutations((1,3,6))) or w in list(itertools.permutations((2,5,5))) or w in list(itertools.permutations((2,6,6))): 
    out = -0.50

  return out


def Pauli():
  Sigma[0] = [[1.0,0],[0,1.0]]
  Sigma[1] = [[0,1.0],[1.0,0]] 
  Sigma[2] = [[0, complex(0,-1.0)],[complex(0,1.0),0]]
  Sigma[3] = [[1,0],[0,-1]] 
  return Sigma


def makeLambda():
  Lambda[0] = [[1,0,0],[0,1,0],[0,0,1]]
  Lambda[1] = [[0,1,0],[1,0,0],[0,0,0]]
  Lambda[2] = [[0,complex(0,-1),0],[complex(0,1),0,0],[0,0,0]]
  Lambda[3] = [[1,0,0],[0,-1,0],[0,0,0]]
  Lambda[4] = [[0,0,1],[0,0,0],[1,0,0]]
  Lambda[5] = [[0,0,complex(0,-1)],[0,0,0],[complex(0,1),0,0]]
  Lambda[6] = [[0,0,0],[0,0,1],[0,1,0]]
  Lambda[7] = [[0,0,0],[0,0,complex(0,-1)],[0,complex(0,1),0]]
  Lambda[8] = [[1.0/(np.sqrt(3)),0,0],[0,1.0/(np.sqrt(3)),0],[0,0,-2/np.sqrt(3)]]

  return Lambda


def action_A_or_Adag(v, in_vec, in_coeff, choice):
# Defines the action of A_ij or A_ij^dagger 

  out_vec = in_vec

  if len(v) != 2 or v[0] > 2 or v[1] > 7: 
    print ("Operator should carry two indices with proper range")
    sys.exit(1)

  if choice == 1:  # Creation
    out_vec[v[0]][v[1]] = in_vec[v[0]][v[1]] + 1
    out_coeff = in_coeff * np.sqrt(out_vec[v[0]][v[1]]) 
  
  elif choice == 0 and out_vec[v[0]][v[1]] > 0:  # Annihilation 
    tmp = np.sqrt((out_vec[v[0]][v[1]]))
    out_vec[v[0]][v[1]] = out_vec[v[0]][v[1]] - 1
    out_coeff = in_coeff * tmp

  elif choice == 0 and out_vec[v[0]][v[1]] == 0:
    print ("WARNING: Annihilation operator acting on |0>", in_vec, v)
    sys.exit(1) 
       


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

        tmp = fSU3((a,b,c)) * LeviCivita(i, j, k)
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

        tmp = fSU3((a,b,c)) * LeviCivita(i, j, k)
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
# How many states should be here? 
# 3420 = 380*3*3 # 380 = fSU3, 3 = i , 3=j 

  st_num = 0
  count = 0 

  for a in range(dimG):
    for b in range(dimG):
      for c in range(dimG):
        for d in range(dimG):
          for e in range(dimG):

            tmp = fSU3((a,b,c)) * fSU3((a,d,e)) 
            # This 'tmp' is non-zero 3420 times..

            if tmp != 0:

              if a == 0:  # 36+36+36+64+64+64+64+16= 380
                count += 1 

              for i in range(3):
                for j in range(3):

                    # TODO! Aug 28 
                    
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

            tmp = fSU3((a,b,c)) * fSU3((a,d,e))
            if tmp != 0:


              for i in range(3):
                for j in range(3):

                  O5state[st_num] = in_state
                  O5norm[st_num] = in_state_norm


                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,e),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,d),O5state[st_num], O5norm[st_num], 1)
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((j,c),O5state[st_num], O5norm[st_num], 1)

                  print ("X", O5state[st_num], "and", i, b) # .... 
                  # TODO # CHECKPOINT 
                  O5state[st_num], O5norm[st_num] = action_A_or_Adag((i,b),O5state[st_num], O5norm[st_num], 0)
                  O5norm[st_num] *= tmp
                  st_num += 1


  return  O5state, O5norm


def OP6(in_state, in_state_norm):  
# f_abc f_ade [(A_ib * A_jc * A†_id * A†_je) + (A_ib * A_id * A†_jc * A†_je) +  \ 
# (A_ib * A_je * A†_id * A†_jc)] | in state, in_state_norm >
# Ex: f_123 f_123 [(A_ib * A_jc * A†_id * A†_je) + (A_ib * A_id * A†_jc * A†_je) +  \ 
# (A_ib * A_je * A†_id * A†_jc)] | in state, in_state_norm >

  st_num = 0

  for a in range(dimG):
    for b in range(dimG):
      for c in range(dimG):
        for d in range(dimG):
          for e in range(dimG):

            tmp = fSU3((a,b,c)) * fSU3((a,d,e))
            if tmp != 0:


              for i in range(3):
                for j in range(3):

                  if st_num < 10:
                    print (a+1, b+1, c+1, a+1, d+1, e+1)


                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]
                  
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c), O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e),O6state[st_num], O6norm[st_num], 0)
                  #print ("with coeff .... ", O6norm[st_num])
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)
                  O6totalstate[st_num] = O6state[st_num]
                  O6totalnorm[st_num] = O6norm[st_num]

                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]

                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e), O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c),O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 0)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)

                  O6totalstate[st_num] += O6state[st_num]
                  O6totalnorm[st_num] += O6norm[st_num]

                  O6state[st_num] = in_state[st_num]
                  O6norm[st_num] = in_state_norm[st_num]

                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,e), O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,d),O6state[st_num], O6norm[st_num], 1)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((j,c),O6state[st_num], O6norm[st_num], 0)
                  O6state[st_num], O6norm[st_num] = action_A_or_Adag((i,b),O6state[st_num], O6norm[st_num], 0)

                  O6totalstate[st_num] += O6state[st_num]
                  O6totalnorm[st_num] += O6norm[st_num]

                  O6totalnorm[st_num] *= tmp
                  st_num += 1

                  #sys.exit(1)


  return  O6totalstate, O6totalnorm


def pretty_print_matrix(matrix):
  print(('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix])))


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")

def Op_2ferm(w):

  if len(w) != 3:
    print ("Needs three indices --> flavor spin color")
    sys.exit(1)

  OF = np.kron(np.kron(Lambda[w[0]], Sigma[w[1]]), Lambda[w[2]])
  return OF 

def Op_4ferm(w1, w2):

  if len(w1) != 3  or len(w2) != 3:
    print ("Needs three indices --> flavor spin color")
    sys.exit(1)

  OF = np.kron(np.kron(Lambda[w1[0]], Sigma[w1[1]]), Lambda[w1[2]])
  OF2 = np.kron(np.kron(Lambda[w2[0]], Sigma[w2[1]]), Lambda[w2[2]])
  res = np.dot(OF, OF2)
  return res

def vev(s1,s2):

  # https://numpy.org/doc/stable/reference/generated/numpy.outer.html 
  # Outer product 

  if len(s1) != len(s2):
    return 0
  else:
    return np.linalg.det(np.einsum('i,j->ij', s1, s2)) 


def di_meson_S_element(P1,P2,O1,O2):

  #out = np.trace(P1 @ O1) * np.trace(P2 @ O2)
  #out +=  np.trace(P1 @ O2) * np.trace(P2 @ O1)
  #out -= np.trace(P1 @ O1 @ P2 @ O2)
  #out -= np.trace(P1 @ O2 @ P2 @ O1)

  Id = np.eye(18)  
  t1 = ncon((Id, Id),([-1,-2],[-3,-4]))
  tmp1 = np.subtract(t1, t1.transpose(0,3,2,1))
  tmp2 = np.subtract(t1, t1.transpose(1,2,3,0)) 
  #out = np.einsum('pm, qn, ik, jl, pkql, imjn', a, b, c, d, tmp1, tmp2)
  out = ncon((P1,P2,O1,O2,tmp1,tmp2),([1,2],[3,4],[5,6],[7,8],[1,6,3,8],[5,2,7,4]))

  return out 


if __name__ == "__main__":


  makeLambda()
  Pauli()

  # Construct O_f like 
  O1 = Op_2ferm((1,1,0))  # T_1 ⊗ \sigma_1 ⊗ I_3
  O2 = Op_2ferm((0,0,0))
  O3 = Op_2ferm((1,0,0))
  O4 = Op_2ferm((0,1,0))


  O5 = np.zeros((18,18))
  O6 = np.zeros((18,18))


  for a in range(1, 9):

    O5 = np.add(O5,Op_2ferm((1,1,a)))
    O6 = np.add(O6,Op_2ferm((0,0,a)))

  
  O7 = np.zeros((18,18))
  O8 = np.zeros((18,18))

  for a in range(1, 9):

    O7 = np.add(O7,Op_2ferm((1,0,a)))
    O8 = np.add(O8,Op_2ferm((0,1,a)))


  O9 = np.zeros((18,18))
  O10 = np.zeros((18,18))

  for j in range (1,3):
    for k in range (1,3):
      if j != k and j == 1:
        O9 = np.add(O9, Op_2ferm((1,j+1,0)))
        O10 = np.add(O10, Op_2ferm((0,k+1,0)))
      if j != k and j == 2:
        O9 = np.add(O9, -1.0*Op_2ferm((1,j+1,0)))
        O10 = np.add(O10, -1.0*Op_2ferm((0,k+1,0)))


  O11 = np.zeros((18,18))
  O12 = np.zeros((18,18))

  for j in range (1,3):
  	for k in range (1,3):

  		if j != k:
  			for a in range(1, 9):
  				O11 = np.add(O11, Op_2ferm((1,j+1,a)))
  				O12 = np.add(O12, Op_2ferm((0,k+1,a)))


  for i in range (6):
    for j in range (6):

      if i == j == 0:
        Smatrix[i][j] = di_meson_S_element(O1, O2, O1, O2) 
      if i == j == 1:
        Smatrix[i][j] = di_meson_S_element(O5, O6, O5, O6)
      if i == j == 2:
        Smatrix[i][j] = di_meson_S_element(O3, O4, O3, O4)
      if i == j == 3:
        Smatrix[i][j] = di_meson_S_element(O7, O8, O7, O8)
      if i == j == 4:
        Smatrix[i][j] = di_meson_S_element(O9, O10, O9, O10)
      if i == j == 5:
        Smatrix[i][j] = di_meson_S_element(O11, O12, O11, O12)

      if i == 0 and j == 1 or i == 1 and j == 0:
        Smatrix[i][j] = di_meson_S_element(O5, O6, O1, O2)
      if i == 0 and j == 2 or i == 2 and j == 0:
        Smatrix[i][j] = di_meson_S_element(O3, O4, O1, O2)
      if i == 0 and j == 3 or i == 3 and j == 0:
        Smatrix[i][j] = di_meson_S_element(O7, O8, O1, O2)
      if i == 0 and j == 4 or i == 4 and j == 0:
        Smatrix[i][j] = di_meson_S_element(O9, O10, O1, O2)
      if i == 0 and j == 5 or i == 5 and j == 0:
        Smatrix[i][j] = di_meson_S_element(O11, O12, O1, O2)

      if i == 1 and j == 2 or i == 2 and j == 1:
        Smatrix[i][j] = di_meson_S_element(O3, O4, O5, O6)
      if i == 1 and j == 3 or i == 3 and j == 1:
        Smatrix[i][j] = di_meson_S_element(O7, O8, O5, O6)

      if i == 2 and j == 3 or i == 3 and j == 2:
        Smatrix[i][j] = di_meson_S_element(O3, O4, O5, O6)
      

  pretty_print_matrix(Smatrix.real)


  #print ("NNZ in TMP2", np.count_nonzero(tmp2))
  #print ("Mat 2f", matprint(tmp2))
  #tmp3 = Op_4ferm((1,0,0),(0,1,0))


  '''

  #out_state1, norm_out_state1 = OP1(vac, 1.0) 
  #out_state2, norm_out_state2 = OP2(vac, 1.0)
  #out_state32, norm_out_state32 = OP3(O2state, O2norm)
  #out_state4, norm_out_state4 = OP4(vac, 1.0) 
  #out_state5, norm_out_state5 = OP5(vac, 1.0) # ISSUE HERE!
  #out_state64, norm_out_state64 = OP6(O4state, O4norm) 
  
  for i in range (10): # Print first 10 states of each type for now. 
    print ("O1 |vac> :", out_state1[i].reshape(size_col), "with coefficient", round(norm_out_state1[i],3)) # 24 states #2-boson state 
    print ("O2 |vac> :", out_state2[i].reshape(size_col), "with coefficient", round(norm_out_state2[i],3)) # 324 states  #3-boson state 
    print ("O3 O2 |vac> :",out_state32[i].reshape(size_col), "with coefficient", round(norm_out_state32[i],3)) # 324 states  #5-boson state 
    print ("O4 |vac> :",out_state4[i].reshape(size_col), "with coefficient", round(norm_out_state4[i],3)) # 3420 states # 4-boson state 
    print ("O5 |vac> :",out_state5[i].reshape(size_col), "with coefficient", round(norm_out_state5[i],3)) # 3420 states # 2-boson state 
    print ("O6 O4 |vac> :",out_state64[i].reshape(size_col), "with coefficient", round(norm_out_state64[i],3)) # 3420 states
    # O6 O4 |vac> : [0. 4. 8. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.] with coefficient 45.798
    # for f_123*f_123 * [(A_12 * A_13 * A†_12 * A†_13) + (A_12 * A_12 * A†_13 * A†_13) +  (A_12 * A_13 * A†_12 * A†_13)] * O4 | vac. >
    # is for i=j=1 
    # for f_123*f_123 * [(A_12 * A_23 * A†_12 * A†_23) + (A_12 * A_12 * A†_23 * A†_23) +  (A_12 * A_23 * A†_12 * A†_23)] * O4 | vac. >
    # i=1, j=2 f_123 * f_123
  '''





