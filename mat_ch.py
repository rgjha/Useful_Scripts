#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
import math
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
import numpy as np
import random
from scipy.linalg import expm
from numpy import linalg as LA
import time 
import datetime 
import sys
startTime = time.time()
print(("STARTED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

if len(sys.argv) < 3:
  print("Usage:", str(sys.argv[0]), "READIN " " SAVE_or_NOT")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])


#************************************************

NCOL = 3   # Rank of gauge group 
g=0.1
c=0.1
kappa=0.0
GENS = NCOL**2 - 1
NMAT = 3 
eps=0.01 
nsteps = int(0.1/eps)
Niters_sim=100
#************************************************

X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM = []
expDS = [] 
EOS = [] 
ACT = [] 
MOM = []

#************************************************

print ("Matrix chain simulation with number of matrices = %4.0f" %(NMAT)) 
print ("NCOL=" "%3.0f " ","  " and 'g' = " " %4.2f" % (NCOL, g)) 
print ("---------------------------------------------------------------------------------")

#************************************************
def dagger(a):

    return np.transpose(a).conj()
#************************************************
def comm(A,B):
    
    return np.dot(A,B) - np.dot(B,A)
#************************************************
def unit_matrix():
    
    matrix = np.zeros((NCOL, NCOL), dtype=complex)
    for i in range (NCOL):
        matrix[i][i] = complex(1.0,0.0)
    return matrix
#************************************************
def copy_fields(b):
    
    for j in range(NMAT):
        X_bak[j] = b[j]

    return X_bak
#************************************************
def rejected_go_back_old_fields(a):
        
    for j in range(NMAT):
        X[j] = a[j]

    return X
#************************************************
def refresh_mom():
    
    for j in range (NMAT):
        mom_X[j] = random_hermitian()

    return mom_X
#************************************************

def pretty_print_matrix(matrix):

    print(('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix])))


#************************************************
# Following Page (33) of 1808.08490

def random_hermitian():

    tmp = np.zeros((NCOL, NCOL), dtype=complex)

    for i in range (NCOL):

        for j in range (i+1, NCOL):

            r1 = np.random.normal(0,1)
            r2 = np.random.normal(0,1)

            tmp[i][j] = complex(r1, r2)/math.sqrt(2)
            tmp[j][i] = complex(r1, -r2)/math.sqrt(2)

    for i in range (NCOL):

        r1 = np.random.normal(0,1)
        tmp[i][i] = complex(r1, 0.0)


    return tmp 

#************************************************
def kinetic_energy(mom_X):

    s = 0.0 
    for j in range (NMAT):
        s += 0.50 * np.trace(np.dot(dagger(mom_X[j]),mom_X[j]))
    return s.real  
#************************************************          
def bosonic_action(X):

    b_action = 0.0 

    #print ("X1 = ", pretty_print_matrix(X[0]))
    #print ("X2 = ", pretty_print_matrix(X[1]))
    #print ("X3 = ", pretty_print_matrix(X[2]))

    for i in range (NMAT):
        b_action -= np.dot(X[i],X[i]) 

    for i in range (NMAT):
        b_action -= (g/NCOL)* np.dot(np.dot(X[i],X[i]),np.dot(X[i],X[i]))

    for i in range (NMAT-1):
        b_action += (2.0*c)* np.dot(X[i],X[i+1]) 

    b_action += kappa*np.dot(X[NMAT-1],X[0])

    #print ("Act is", np.trace(b_action).real)

    return np.trace(b_action).real
#************************************************

def bosonic_force(X): 

    tmp_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
    
    for i in range (NMAT): 

        tmp_X[i] = np.zeros((NCOL, NCOL), dtype=complex)

        for j in range (NMAT):

            ..... TODO .... 

        f_X[i] = ...

# Make AH

        for k in range (NCOL):

            f_X[i][k][k] =  complex(0.0,0.0)

            for l in range (k+1, NCOL):

                tr = 0.5 * (f_X[i][k][l].real - f_X[i][l][k].real)
                f_X[i][k][l] = complex(tr, f_X[i][k][l].imag)
                f_X[i][l][k] = complex(-tr, f_X[i][l][k].imag)
                tr = 0.5 * (f_X[i][k][l].imag + f_X[i][l][k].imag)
                f_X[i][k][l] = complex(f_X[i][k][l].real,tr) 
                f_X[i][l][k] = complex(f_X[i][l][k].real,tr) 

        # Check if forces f_X are hermtian => should be!
        # Tracelessness need not be checked as imposed some lines above
        if LA.norm(dagger(f_X[i])-f_X[i]) > 1e-14: 
            print ("f_X is not H", LA.norm(dagger(f_X[i])+f_X[i]))


    return f_X 

#************************************************

def leapfrog(X,mom_X, eps):


    for j in range(NMAT):
        
        X[j] = X[j] + (mom_X[j] * eps/2.0)

    f_X = bosonic_force(X)


    # Middle steps 
    for step in range(nsteps):

        for j in range(NMAT):

            
            mom_X[j] = mom_X[j] + (f_X[j] * eps) 
            X[j] = X[j] + (mom_X[j] * eps)
            
        f_X = bosonic_force(X)

    # Last 'half' step

    for j in range(NMAT):
        
        mom_X[j] = mom_X[j] + (f_X[j] * eps)
        X[j] = X[j] + (mom_X[j] * eps/2.0)
        

    return X, mom_X, f_X


def update(X):

    mom_X = refresh_mom()
    KE = kinetic_energy(mom_X)
    ba = bosonic_action(X)
    start_act =  ba + KE
    sys.exit(1)
    #print("start action: commutator = " , ba , "and gauge+scalar momenta =" , KE)
    X_bak = copy_fields(X) 
    X, mom_X, f_X = leapfrog(X,mom_X,eps)
    KE = kinetic_energy(mom_X)
    ba = bosonic_action(X)
    end_act = ba + KE
    #print("end action: commutator = " , ba , "and gauge+scalar momenta =" , KE)
    change = end_act - start_act
    HAM.append(abs(change))
    expDS.append(np.exp(-1.0*change))   
    # <exp(Hold-Hnew)> = 1 # https://www.osti.gov/servlets/purl/6871614

    if np.exp(-change) < random.uniform(0,1):
        X = rejected_go_back_old_fields(X_bak)
        print(("REJECT: deltaS = " "%8.7f " " startS = " "%8.7f" " endS = " "%8.7f" % (change, start_act, end_act)))
    else:   
        print(("ACCEPT: deltaS = " "%8.7f " "startS = " "%8.7f" " endS = " "%8.7f" % (change, start_act, end_act)))


    ACT.append(ba)
    MOM.append(float(KE/GENS))

    if MDTU%1 == 0:

        tmp = float(ba/(NCOL*NCOL))
        f4.write("%4.8f\n" % tmp)

    return X


#***************The main routine****************************

if __name__ == '__main__':


    if READIN ==0:

        for i in range (NMAT):  
            X[i] = random_hermitian()
            # Checked that X = X^dagger
            #print (pretty_print_matrix(X[i])) 
            #sys.exit(1)
            
    if READIN ==1:
        print ("Reading old config.")

        with open("config.txt") as f2:
            A = np.loadtxt(f2).view(complex)
        f2.close()

        for i in range (NMAT):
            for a in range (NCOL):
                for b in range (NCOL):

                    X[i][a][b] = A[(NCOL*i)+a][b] 


    f4 = open("action.txt", "w")

    for MDTU in range (Niters_sim):
        X = update(X) 

    f4.close()


    if SAVE ==1:

        print ("Saving config.")
        f1 = open("config.txt", "w")
        for i in range (NMAT):
            np.savetxt(f1, X[i].view(float), delimiter= " ")  
        f1.close()

    ACT = [x/GENS for x in ACT]
      
    MDTU = np.linspace(0, Niters_sim, Niters_sim, endpoint=True)
    plt.ylabel(r'$\Delta$ Action')
    plt.xlabel('MDTU')
    plt.figure(1)
    plot(MDTU, ACT)
    plot(MDTU, MOM) 
    plt.show()

    print(("<Action>", np.mean(ACT), "+/-", (np.std(ACT)/np.sqrt(np.size(ACT) - 1.0))))
    print(("<exp(-deltaH)>", np.mean(expDS), "+/-", np.std(expDS)/np.sqrt(np.size(expDS) - 1.0))) # This should be 1 within errors!
    print ("Delta 'H' is", np.mean(HAM), "+/-", np.std(HAM)/np.sqrt(np.size(HAM) - 1.0))  # This should scale with dtau!
    print(("COMPLETED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
