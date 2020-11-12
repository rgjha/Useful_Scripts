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
from numpy.linalg import matrix_power
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

NCOL = 100   # Rank of gauge group 
Niters_sim=500
g=1.1
c=0.1
kappa=0.0
NMAT = 1 
eps=0.011 
nsteps = 20

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

def box_muller():

    PI = 2.0*math.asin(1.0);    
    r = random.uniform(0,1)
    s = random.uniform(0,1)
    p = np.sqrt(-2.0*np.log(r)) * math.sin(2.0*PI*s)
    q = np.sqrt(-2.0*np.log(r)) * math.cos(2.0*PI*s)

    return p,q

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

            r1, r2 = box_muller()

            tmp[i][j] = complex(r1, r2)/math.sqrt(2)
            tmp[j][i] = complex(r1, -r2)/math.sqrt(2)

    for i in range (NCOL):

        r1, r2 = box_muller()
        tmp[i][i] = complex(r1, 0.0)


    return tmp 

#************************************************

def hamil(X,mom_X):

    ham = action(X) 
    for j in range (NMAT):
        ham += 0.50 * np.trace(np.dot(mom_X[j],mom_X[j])).real 

    return ham  

#************************************************ 

def action(X):

    b_action = 0.0 

    for i in range (NMAT):
        b_action = 0.50 * np.trace(np.dot(X[i],X[i])).real  
        b_action += (g/4.0)* np.trace((matrix_power(X[i], 4))).real

    return b_action*NCOL

#************************************************

def bosonic_force(X): 

    f_X = np.zeros((NMAT, NCOL, NCOL), dtype=complex)
    
    for i in range (NMAT): 

        f_X[i] = (X[i] + (g*(matrix_power(X[i], 3))))*NCOL

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


        if MDTU%2 == 0:

            tmp = np.trace(np.dot(X[i],X[i])).real
            ACT.append(tmp/NCOL)


    if SAVE ==1:

        print ("Saving config.")
        f1 = open("config.txt", "w")
        for i in range (NMAT):
            np.savetxt(f1, X[i].view(float), delimiter= " ")  
        f1.close()

      
    MDTU = np.linspace(0, Niters_sim, Niters_sim, endpoint=True)
    plt.ylabel(r'Tr(Xsq)')
    plt.xlabel('MDTU')
    plt.figure(1)
    plot(MDTU, ACT)
    plt.show()

    print(("<Tr X^2 / NCOL>", np.mean(ACT), "+/-", (np.std(ACT)/np.sqrt(np.size(ACT) - 1.0))))
    print(("COMPLETED: " , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))


