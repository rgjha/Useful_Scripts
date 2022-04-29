# Find GCD using Euclid's algorithm
# Classical segment of Shor's algorithm 

# "[The Euclidean algorithm] is the granddaddy of all algorithms, because it is the oldest nontrivial 
# algorithm that has survived to the present day."
# - Donald Knuth, The Art of Computer Programming, Vol. 2: Seminumerical Algorithms, 
# 2nd edition (1981), p. 318. 

import sys 

if len(sys.argv) < 2:
  print("Usage:", str(sys.argv[0]), "a " " b"), sys.exit(1)
  
a, b, iter = int(sys.argv[1]), int(sys.argv[2]), 0

def GCD(x , y, iter):

    if y == 0:
        print ("GCD found in", iter, "iterations") 
        return x 
    iter += 1 
    return GCD(y , int(x % y), iter)


print ("GCD is", GCD(a, b, iter))
