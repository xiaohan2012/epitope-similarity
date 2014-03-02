import math

def dotproduct(v1, v2):
    """
    dot product between two vector
    
    Paramter:
    v1: array of number
    v2: array of number

    Return: 
    number
    """
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    """
    the norm of vector v, aka, |v|
    """
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    """
    the angle between two vector.

    Paramter:
    v1: array of number
    v2: array of number

    Return: 
    float
    """
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
  
