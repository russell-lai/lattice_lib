# General utility functions

from sage.all import *

def next_power_two(x):
    # Input: integer x
    # Output: the next power-of-2 y of x 
    y = 2
    while y < x:
        y *= 2
    return y

def prime_power_divisors(f):
    # Input: integer f
    # Output: list containing all prime-power divisors of f
    factors = factor(f)
    return [factors[i][0]**factors[i][1] for i in range(len(factors))]

def max_prime_power_divisor(f):
    return max(prime_power_divisors(f))

def is_signed_perm(A,B):
    # Check if the input lists A and B are signed permutations of each other.
    assert len(A) == len(B)
    A_abs = [abs(a) for a in A]
    B_abs = [abs(b) for b in B]
    A_abs.sort()
    B_abs.sort()
    return A_abs == B_abs

def is_local_signed_perm(A,B,chunk_size):
    # Check if the input lists A and B satisfy the following condition: B is obtained by first permuting the chunks of A, and then performing a signed permutation within each chunk.
    assert len(A) == len(B)
    assert len(A) % chunk_size == 0    
    num_chunks = int(len(A)/chunk_size)
    A_abs_blocks = [[abs(a) for a in A[i*chunk_size : (i+1)*chunk_size]] for i in range(num_chunks)]
    B_abs_blocks = [[abs(b) for b in B[i*chunk_size : (i+1)*chunk_size]] for i in range(num_chunks)]
    for block in A_abs_blocks:
        block.sort()
    for block in B_abs_blocks:
        block.sort()
    A_abs_blocks.sort()
    B_abs_blocks.sort()
    return A_abs_blocks == B_abs_blocks