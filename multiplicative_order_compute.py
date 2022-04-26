import numpy as np
from math import * 

def multiplicative_group(n):

# Returns the multiplicative group (MG) modulo n.
# n: Modulus of the MG.
                
    assert n > 1
    group = [1]
    for x in range(2, n):
        if gcd(x, n) == 1:
            group.append(x)
    return group

n = 21
print(f"The multiplicative group modulo n = {n} is:")
print(multiplicative_group(n))
