import math
from collections import Counter

# Computes the entropy of a image link string. 
# We can do fancier things but this will work well.
# The likelihood of two links which are different 
# with same entropy is close to zero. 

def compute_entropy(string):
    # Count the frequency of each character in the string
    char_counts = Counter(string)
    
    # Total length of the string
    total_length = len(string)
    
    entropy = 0.0
    for count in char_counts.values():
        # Compute the probability of each character
        probability = count / total_length
        
        # Compute the contribution to entropy
        entropy += probability * math.log2(probability)
    
    # Negate the entropy to obtain the final result
    entropy = -entropy
    
    return entropy

# Example usage
input_string = input("Enter a string: ")
entropy = compute_entropy(input_string)
print("Entropy:", round(entropy,5))
