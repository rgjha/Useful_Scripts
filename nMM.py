#!/usr/bin/python
# -*- coding: utf-8 -*-
# Last edit: January 26, 2021  
# Group integration over compact Lie groups for multi-matrix models 
# especially looking at chains with X_n+1 = X_1 (i.e. kappa != 0) 
# S = sum_{i=1}^{i=3} (1/2)*Tr(X_{i}^2) + (g/4.0)*Tr(X_{i}^4) + c Tr(X_0 X_1 + X_1 X_2) + kappa Tr(X_2 X_0) 


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
print ("STARTED:" , datetime.datetime.now().strftime("%d %B %Y %H:%M:%S"))


if len(sys.argv) < 3:
  print("Usage:", str(sys.argv[0]), "READIN " " SAVE_or_NOT")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])


#************************************************

NCOL = 200  # Size of matrix 
Niters_sim = 100
g = 1 
c = 1.35
kappa = 1.35 
NMAT = 3
eps=3e-4   
nsteps = 30

#************************************************

X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
mom_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
X_bak = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
HAM = []
expDS = [] 
EOS = [] 
ACT = [] 
ACT1 = []
ACT2 = []
MOM = []

#************************************************

print ("Matrix chain simulation with number of matrices = %4.0f" %(NMAT))
print ("NCOL=" "%3.0f " ","  " and '(g, c, kappa)' = " " (%4.2f,%4.2f,%4.2f)" % (NCOL, g, c1, kappa))
print ("---------------------------------------------------------------------------------")

#************************************************
def dagger(a):

    return np.transpose(a).conj()
#************************************************
def box_muller():

    PI = 2.0*math.asin(1.0);
    r = random.uniform(0,1)
    s = random.uniform(0,1)
    p = np.sqrt(-2.0*np.log(r)) * math.sin(2.0*PI*s)
    q = np.sqrt(-2.0*np.log(r)) * math.cos(2.0*PI*s)

    return p,q
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

            r1, r2 = box_muller()

            tmp[i][j] = complex(r1, r2)/math.sqrt(2)
            tmp[j][i] = complex(r1, -r2)/math.sqrt(2)

    for i in range (NCOL):
        
        r1, r2 = box_muller()
        tmp[i][i] = complex(r1, 0.0)

    return tmp

#****************************************

def random_unitary(N):
    """Return a Haar distributed random unitary from U(N)"""
    Z = np.random.randn(N,N) + 1.0j * np.random.randn(N,N)
    [Q,R] = sp.linalg.qr(Z)
    D = np.diag(np.diagonal(R) / np.abs(np.diagonal(R)))
    return np.dot(Q, D)
#****************************************

def hamil(X,mom_X):

    ham = action(X)
    for j in range (NMAT):
        ham += 0.50 * np.trace(np.dot(mom_X[j],mom_X[j])).real

    return ham

#****************************************


def kinetic_energy(mom_X):

    s = 0.0 
    for j in range (NMAT):
        s += 0.50 * np.trace(np.dot(dagger(mom_X[j]),mom_X[j]))  # mom * mom^dag also works!
    return s.real  
#***************************************           
def action(X):

    b_action = 0.0

    for i in range (NMAT):
        b_action += 0.50 * np.trace(np.dot(X[i],X[i])).real
        b_action += (g/4.0)* np.trace((matrix_power(X[i], 4))).real

        if i == NMAT-1:
            pre = kappa
        else:
            pre = c

        b_action -= pre*np.trace(np.dot(X[i],X[(i+1)%NMAT])).real

    return b_action*NCOL
#************************************************

def bosonic_force(X):

    f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)

    for i in range (NMAT):

        if i == 0:
            f_X[i] = (X[i]+(g*(matrix_power(X[i], 3))) - (c1*X[(i+1)%NMAT]) - (kappa*X[(i-1+NMAT)%NMAT]))*NCOL

        if i == 1:
            f_X[i] = (X[i]+(g*(matrix_power(X[i], 3))) - (c1*X[(i+1)%NMAT]) - (c1*X[(i-1+NMAT)%NMAT]))*NCOL

        if i == 2:
            f_X[i] = (X[i]+(g*(matrix_power(X[i], 3))) - (kappa*X[(i+1)%NMAT]) - (c1*X[(i-1+NMAT)%NMAT]))*NCOL

    return f_X 

#************************************************

def leapfrog(X,eps):

    mom_X = refresh_mom()
    ham_init = hamil(X,mom_X)


    for j in range(NMAT):
        X[j] = X[j] + (mom_X[j] * eps * 0.50)

    for j in range(NMAT):

        for i in range(1, nsteps):
            f_X = bosonic_force(X)
            mom_X[j] = mom_X[j] - (f_X[j]*eps)
            X[j] = X[j] + (mom_X[j]*eps)


    f_X = bosonic_force(X)

    for j in range(NMAT):

        mom_X[j] = mom_X[j] - (f_X[j] * eps)
        X[j] = X[j] + (mom_X[j] * eps * 0.50)


    ham_final = hamil(X,mom_X)
    
    return X, ham_init, ham_final

#***************The main routine****************************

if __name__ == '__main__':
    
    
    if READIN ==0:
        for i in range (NMAT):
            for j in range (NCOL):
                for k in range (NCOL):
                    X[i][j][k] = complex(0.0,0.0)
                    
    if READIN ==1:
        print ("Reading old config.")
        with open("config.txt") as f2:
            A = np.loadtxt(f2).view(complex)
        f2.close()

        for i in range (NMAT):
            for j in range (NCOL):
                for k in range (NCOL):
                    X[i][j][k] = A[(NCOL*i)+j][k]
                    
                    
     for MDTU in range (Niters_sim):

        X_bak = copy_fields(X)
        X, start, end = leapfrog(X, eps)

        change = end - start
        # <exp(Hold-Hnew)> = 1 # https://www.osti.gov/servlets/purl/6871614

        if np.exp(-change) < random.uniform(0,1):
            X = rejected_go_back_old_fields(X_bak)
            print(("REJECT: deltaS = " "%8.7f " " startS = " "%8.7f" " endS = " "%8.7f" % (change, start, end)))
        else:
            print(("ACCEPT: deltaS = " "%8.7f " "startS = " "%8.7f" " endS = " "%8.7f" % (change, start, end)))
            #continue  

        tmp0 = np.trace(np.dot(X[0],X[0])).real
        ACT.append(tmp0/NCOL)    
        
        if NMAT == 2:
            tmp1 = np.trace(np.dot(X[1],X[1])).real
            ACT1.append(tmp1/NCOL)

        if NMAT == 3:
            tmp1 = np.trace(np.dot(X[1],X[1])).real
            ACT1.append(tmp1/NCOL)
            tmp2 = np.trace(np.dot(X[2],X[2])).real
            ACT2.append(tmp2/NCOL)
            
            
    if SAVE ==1:
        
        print ("Saving config.")
        f1 = open("config.txt", "w")
        for i in range (NMAT):
            np.savetxt(f1, X[i].view(float), delimiter= " ")
        f1.close()

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    MDTU = np.linspace(0, Niters_sim, Niters_sim, endpoint=True)
    plt.ylabel(r'Tr($X_{1,2,3}^{2}$)')
    plt.xlabel('No. of iteration')
    plt.title(r"$\rm{g} = 1, \rm{c} = 1.35, \rm{\kappa} = 1.35, $N$ = 200$", fontsize=16, color='black')
    plt.figure(1)


    if NMAT == 2:
        plot(MDTU, ACT)
        plot(MDTU, ACT1)


    if NMAT == 3:
        plot(MDTU, ACT)
        plot(MDTU, ACT1)
        plot(MDTU, ACT2)

    plt.show()
    #plt.savefig('plot_3MM_chain.pdf')
    print(("<Tr X1^2 / NCOL>", np.mean(ACT), "+/-", (np.std(ACT)/np.sqrt(np.size(ACT) - 1.0))))
    print(("<Tr X2^2 / NCOL>", np.mean(ACT1), "+/-", (np.std(ACT1)/np.sqrt(np.size(ACT1) - 1.0))))
    print(("<Tr X3^2 / NCOL>", np.mean(ACT2), "+/-", (np.std(ACT2)/np.sqrt(np.size(ACT2) - 1.0))))
    print(("COMPLETED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    
    
    


    MDTU = np.linspace(0, Niters_sim, Niters_sim, endpoint=True)
    plt.ylabel(r'Tr($X_{1,2,3}^{2}$)')
    plt.xlabel('No. of iteration')
    plt.title(r"$\rm{g} = 1, \rm{c} = 1.35, \rm{\kappa} = 1.35, $N$ = 200$", fontsize=16, color='black')
    plt.figure(1)    
