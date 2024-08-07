# Tools related to subtractive sets

from sage.all import *
from util import *
var('z')

def smallest_norm(f, print_ideals = False, verbose = False):
    # Input: conductor f
    # Output: the smallest ideal norm r in the f-th cyclotomic field. 
    # If print_ideals = True, output also the ideals of norm r. 
    # If verbose = True, more information is printed during the computation. 
    # Note: Slow for f > 200
    assert f > 1
    phi = euler_phi(f)
    K = CyclotomicField(f)
    p = K.next_split_prime()
    if verbose:
        print("Smallest fully splitting prime is {}.".format(p))
    zeta_coeffs = list(K.zeta_coefficients(p))
    if verbose:
        print("The first {} coefficients of the Dedekind zeta function are {}.".format(p, zeta_coeffs))
    r = 2
    while zeta_coeffs[r-1] == 0:
        r += 1
    # the smallest ideal is of norm r

    assert is_prime_power(r)
    if print_ideals:
        ideal_factorisation = K.ideal(prime_divisors(r)[0]).factor()
        return r, [[ideal, multiplicity] for ideal, multiplicity in ideal_factorisation]
    else:
        return r

def subtractive_set(f):
    # Input: Conductor f
    # Output: Subtractive set over the f-th cyclotomic ring
    if is_prime_power(f):
        return subtractive_set_prime_power(f)
    else:
        return subtractive_set_non_prime_power(f)

def subtractive_set_prime_power(f):
    ## Input: Prime-power conductor f
    ## Output: Set C containing elements of the form (z**i-1)/(z-1) with i ranging from 0 to phi(rad(f))-1
    f_max = max_prime_power_divisor(f)
    K = CyclotomicField(f, 'z')
    C = [(z**i-1)/(z-1) for i in range(euler_phi(radical(f)))]
    return C

def subtractive_set_non_prime_power(f):
    ## Input: Non-prime-power conductor f
    ## Output: Set C containing first f/f_max roots of unity, where f_max is the maximum prime-power divisor of f 
    f_max = max_prime_power_divisor(f)
    K = CyclotomicField(f, 'z')
    C = [z**i for i in range(f/f_max)]
    return C

def subtractive_set_expansion_factor(C, t, trials=50):
    ## Input: Subtractive set C, recovery threshold t
    ## Output: Estimation of expansion factor of C with recovery threshold t
    beta_max = 0
    for _ in range(trials):
        S = [C[i] for i in Combinations(ZZ(len(C)),t).random_element()]
        V = matrix([[S[i]**j for i in range(t)] for j in range(t)])
        w = ~V * vector([1]+[0 for _ in range(t-1)])
        L = []
        for i in range(t):
            L += w[i].complex_embeddings()
        beta = vector(L).norm()
        if beta > beta_max:
            beta_max = beta
    return beta_max

def is_integral_hack(s):
    return (str(s)).find("/") == -1

def is_unit_hack(s):
    return is_integral_hack(1/s)

def is_subtractive(C, hack=True):
    ## Input: Set C over some number field K
    ## Output: Decide whether C is a subtractive set over the ring of integers of K
    if hack:
        for si in C:
            for sj in C:
                if si == sj:
                    continue
                if not is_unit_hack(si - sj):
                    return False
        return True
    else:
        K = C[0].parent()
        R = K.ring_of_integers()
        for c in C:
            for c_ in C:
                if c == c_:
                    continue
                if not R(c-c_).is_unit():
                    return False
        return True



