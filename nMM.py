#!/usr/bin/python
# -*- coding: utf-8 -*-
# WORK IN PROGRESS! 
# Group integration over compact Lie groups for multi-matrix models 
# especially looking at chains with X_n+1 = X_1 

import random
import math
from matplotlib.pyplot import *
from matplotlib import pyplot as plt
import numpy as np
import random
from scipy.linalg import expm
from numpy import linalg as LA
import np.random
import scipy as sp
import scipy.linalg
import time 
import datetime 
import sys
startTime = time.time()
print(("STARTED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


#************************************************

NCOL = 3  # Size of matrix 
c= 0.1
g= 0.1
NMAT = 10
eps=0.001   
nsteps = int(0.01/eps)

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

print ("Matrix model simulation")
print ("N=" "%3.0f " ","  " and coefficient = " " %4.2f" % (NCOL, COUPLING)) 
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
        mom_X[j] = random_anti_hermitian()

    return mom_X
#************************************************

def pretty_print_matrix(matrix):

    print(('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix])))

#************************************************
# Construct random hermitian matrices 

def random_hermitian(): 

    tmp = np.zeros((NCOL, NCOL), dtype=complex)
 
    for i in range (NCOL):
        for j in range (i+1, NCOL): 

            tmp[i][j] = complex(np.random.normal(0,1), np.random.normal(0,1))/math.sqrt(2.0)
            tmp[j][i] = complex(tmp[i][j].real, -1.0*tmp[i][j].imag)

    for i in range (NCOL):
        for j in range (NCOL): 

            if i ==j:
                tmp[i][i] = complex(np.random.normal(0,1),0.0)

    return tmp

#****************************************


def random_unitary(N):
    """Return a Haar distributed random unitary from U(N)"""
    Z = np.random.randn(N,N) + 1.0j * np.random.randn(N,N)
    [Q,R] = sp.linalg.qr(Z)
    D = np.diag(np.diagonal(R) / np.abs(np.diagonal(R)))
    return np.dot(Q, D)



#****************************************

def kinetic_energy(mom_X):

    s = 0.0 
    for j in range (NMAT):
        s += 0.50 * np.trace(np.dot(dagger(mom_X[j]),mom_X[j]))  # mom * mom^dag also works!
    return s.real  
#***************************************           
def bosonic_action(X):

    b_action = 0.0 
    for i in range (NMAT):
        
        tmp1 = np.dot(X[i],X[i])
        tmp2 = np.dot(tmp1, tmp1)
        b_action += np.trace(tmp1) + ((G/NCOL)* np.trace(tmp2))
    
    for i in range (NMAT-1):
        tmp1 = 2*c*np.trace(np.dot(X[i], X[i+1]))
        b_action -= tmp1
    
    return b_action
#************************************************

def bosonic_force(X):

    # Force from commutator term -- sum{j} [X_j, [X_i, X_j]]^dagger * 2 * coupling  

    tmp_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
    
    for i in range (NMAT):


TODO: 
            else:
                temp = comm(X[i], X[j])
                tmp_X[i] += comm(X[j], temp)

        # Alt: 2.0 * np.dot(np.dot(X[j], X[i]),X[j]) -  np.dot(np.dot(X[j], X[j]), X[i]) - np.dot(np.dot(X[i], X[j]), X[j])
        f_X[i] = 2.0*COUPLING*dagger(tmp_X[i])

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



        # Check if forces f_X are traceless and AH 
        if LA.norm(dagger(f_X[i])+f_X[i]) > 1e-14: 
            print ("f_X is not AH", LA.norm(dagger(f_X[i])+f_X[i]))


    return f_X 

#************************************************

def leapfrog(X,mom_X, eps):

    f_X = bosonic_force(X)

    for j in range(NMAT):
        
        mom_X[j] = mom_X[j] + (f_X[j] * eps/2.0) 
        X[j] = X[j] + (mom_X[j] * eps)


    f_X = bosonic_force(X)


    # Middle steps 
    for step in range(nsteps-1):

        for j in range(NMAT):

            X[j] = X[j] + (mom_X[j] * eps)
            mom_X[j] = mom_X[j] + (f_X[j] * eps) 
            

        f_X = bosonic_force(X)

    # Last 'half' step

    for j in range(NMAT):
        
        mom_X[j] = mom_X[j] + (f_X[j] * eps/2.0)
        X[j] = X[j] + (mom_X[j] * eps)
        

    return X, mom_X, f_X


def update(X):

    mom_X = refresh_mom()
    KE = kinetic_energy(mom_X)
    ba = bosonic_action(X)
    start_act =  ba + KE
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
    MOM.append(KE)

    if MDTU%1 == 0:

        tmp = float(ba/(NCOL*NCOL))
        f4.write("%4.8f\n" % tmp)

    return X


#***************The main routine****************************

if __name__ == '__main__':

    A = random_hermitian()
    print(pretty_print_matrix(A))


    '''
    f4 = open("action.txt", "w")

    for MDTU in range (Niters_sim):
        X = update(X) 

    f4.close()

    MDTU = np.linspace(0, Niters_sim, Niters_sim, endpoint=True)
    plt.ylabel(r'$\Delta$ Action')
    plt.xlabel('MDTU')
    plt.figure(1)
    plot(MDTU, ACT)
    plot(MDTU, MOM) 
    plt.show()

    ACT = [x/GENS for x in ACT]


    print(("<Action>", np.mean(ACT), "+/-", (np.std(ACT)/np.sqrt(np.size(ACT) - 1.0))))
    print(("<exp(-deltaH)>", np.mean(expDS), "+/-", np.std(expDS)/np.sqrt(np.size(expDS) - 1.0))) # This should be 1 within errors!
    print ("Delta 'H' is", np.mean(HAM), "+/-", np.std(HAM)/np.sqrt(np.size(HAM) - 1.0))  # This should scale with dtau!
    print(("COMPLETED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    '''
