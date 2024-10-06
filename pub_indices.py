import os
import sys
import numpy as np
# Last updated: 06 October 2024 
# Computes h-index and g-index given a list with citations of papers 

#cit = [47, 14, 16, 3, 3, 4, 13, 52, 11, 1, 14, 20, 18, 14, 5, 6, 1, 5, 11, 8, 0, 2, 1, 1, 4, 1, 1.] as on 09/03/24 
#cit = [47, 14, 16, 3, 3, 4, 13, 54, 11, 2, 17, 21, 18, 16, 5, 6, 1, 7, 11, 8, 0, 3, 1, 1, 4, 1, 1, 0, 0.] # as on 29/04/24
#cit = [48, 14, 16, 3, 3, 4, 13, 58, 12, 2, 18, 21, 18, 18, 5, 7, 1, 8, 12, 8, 0, 3, 3, 3, 9, 1, 2, 2, 0, 3, 3, 1.] # as on 07/08/24
#cit = [48, 14, 16, 3, 3, 4, 13, 59, 12, 2, 18, 21, 18, 18, 5, 7, 1, 8, 12, 8, 0, 3, 3, 4, 9, 1, 2, 2, 0, 3, 3, 2.] # as on 09/20/24
cit = [48, 14, 16, 3, 3, 4, 14, 59, 12, 2, 19, 21, 18, 19, 5, 7, 1, 9, 12, 8, 0, 3, 4, 5, 10, 1, 2, 2, 0, 3, 4, 2.] 




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
