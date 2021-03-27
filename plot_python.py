import sys
import math
from math import sqrt
import numpy as np
import scipy as sp  
from numpy import ndarray
from matplotlib import pyplot as plt
import time
import datetime
import itertools

if len( sys.argv ) == 3 :
    filename = sys.argv[1]
    filename1 = sys.argv[2]
if len( sys.argv ) == 0 or len( sys.argv ) != 3:
    print("Requires one argument : FILE ")
    print("Run this as : python xxx.py output_num digitize ")
    sys.exit()

beta=[]
dfdT=[]
beta1=[] 
d2FdT2 = [] 

file = open(filename, "r")
for line in itertools.islice(file, 0, None):
    line = line.split()
    beta.append(float(line[0]))
    dfdT.append(float(line[1]))


file1 = open(filename1, "r")
for line in itertools.islice(file1, 0, None):
    line = line.split()
    beta1.append(float(line[0]))
    d2FdT2.append(float(line[1]))



plt.rc('text', usetex=True)
plt.rc('font', family='serif')

f = plt.figure()
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$\beta$',fontsize=13)
ax1.set_ylabel(r'$\langle S_{\rm{this ~work}} \rangle$', color=color,fontsize=13)
ax1.plot(beta, dfdT, marker="*", color=color, label=r'$\chi=61$')

plt.grid(which='major', axis='y', linestyle='--')
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx() 
color2 = 'tab:blue'
ax2.set_ylabel(r'$\langle S_{\rm{arXiv:1302.1908}} \rangle$', color=color2,fontsize=13) 
ax2.plot(beta1, d2FdT2, marker="o", color=color2, label=r'$\chi=22$')
ax2.tick_params(axis='y', labelcolor=color2)
plt.grid(which='minor', axis='x', linestyle='--')
plt.title(r"3d classical XY model using",fontsize=16, color='black')
plt.title(r"$\chi = 50, D_{n} = 10, \rm{Vol.} = 2^{15}$",fontsize=16, color='black')
fig.tight_layout()
plt.show()
plt.savefig('plot_name.pdf')
