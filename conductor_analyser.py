# Tool for analysing conductors

from sage.all import *
import csv

def analyse(f):
    # Input: conductor f
    # Output: information of interest about conductor f
    phi = euler_phi(f)
    pmax = max(prime_divisors(f))
    fmax = max(prime_power_divisors(f))
    twist_degree = None
    twist_packing_ratio = None
    attributes = []
    if is_prime(f):
        conductor_type = "prime"
        subtractive_set_size = pmax
    elif is_prime_power(f):
        conductor_type = "prime-power"
        subtractive_set_size = pmax
    else:    
        conductor_type = "non-prime-power"
        subtractive_set_size = ZZ(f/fmax)
    if is_power_of_two(f):
        attributes += ["power-2"]
        twist_degree = phi
        twist_packing_ratio = 1
    elif is_squarefree(odd_part(f)):
        attributes += ["square-free odd part"]
        twist_degree = euler_phi(f/odd_part(f)) * prod([euler_phi(ppf)/2 for ppf in prime_power_divisors(odd_part(f))])  ## Degree of the subring supporting twisted-trace-based inner product computation
        twist_packing_ratio = phi/twist_degree  ## Ratio between cyclotomic ring degree and subring degree, smaller is better 
    if is_odd(f) and not is_squarefree(f):
        attributes += ["odd non-square-free"]

    ratio_degree_next_power_two = n(phi/next_power_two(phi), digits=3)
    ratio_subtractive_set_degree = n(subtractive_set_size/phi, digits=3)
    return [
        f, 
        len(factor(f)),
        pmax, 
        fmax, 
        phi, 
        twist_degree,
        twist_packing_ratio,
        subtractive_set_size,
        next_power_two(phi),
        ratio_subtractive_set_degree,
        ratio_degree_next_power_two, 
        factor(f), 
        conductor_type,
        *attributes
        ]

def gen_conductor_table(bound = 2**15+1):
    # Generate a file conductors.csv recording information about conductors from 2 to bound (default = 2**15+1).
    results = [[
        "conductor", 
        "#prime factors",
        "smoothness", 
        "power smoothness", 
        "degree",
        "twist degree",
        "twist packing ratio",
        "subtractive set size", 
        "next power-2",
        "subtractive set size/degree",
        "degree/next power-2", 
        "factorisation",
        "type",
        "attributes"
        ]]
    for f in range(2,bound):
        result = analyse(f)
        results += [result]

    with open('conductors.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(results)