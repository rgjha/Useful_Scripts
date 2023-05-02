import os
import sys
import numpy as np
# Last updated: March 11, 2023
# Computes h-index and g-index given a list with citation of papers 

cit = [43, 12, 14, 3, 3, 4, 13, 45, 11, 1, 10, 18, 15, 11, 2, 5, 1, 1, 9, 6, 0.]
n = np.shape(cit)[0]

print ("Publication record of Raghav Govind Jha as per iNSPIRE-HEP") 

def h_index(citations):
    
    citations = np.array(citations)
    n         = citations.shape[0]
    array     = np.arange(1, n+1)
    
    # reverse sorting
    citations = np.sort(citations)[::-1]
    
    # intersection of citations and k
    h_idx = np.max(np.minimum(citations, array))
    
    return h_idx

def g_index(citations):
    citations_sorted = sorted(citations, reverse=True)
    n = len(citations_sorted)
    sum = 0
    for i in range (n):

        sum += citations_sorted[i] 
        if sum <= (i+1)**2: 
            out = i
            break 

    return out


if __name__ == '__main__':

    print ("h-index", h_index(cit)) 
    print ("g-index", g_index(cit)) 
    print ("Total citations", round(np.sum(cit)))
    print ("Average citation per paper", round(np.sum(cit)/n,2))
