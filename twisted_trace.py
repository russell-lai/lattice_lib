# Tools related to twisted trace

from sage.all import *
from sage.structure.element import Vector

def power_basis(f):
    K = CyclotomicField(f)
    B = K.integral_basis()
    return B

def real_power_basis(f):
    phi = euler_phi(f)
    K = CyclotomicField(f)
    z = K.gen()
    B = [1] + [z**i + z**(-i) for i in range(1,phi/2)]
    return B

def real_prefix_basis(f):
    phi = euler_phi(f)
    K = CyclotomicField(f)
    z = K.gen()
    B = [sum([z**j + z**(-j) for j in range(phi/2-i,phi/2+1)]) for i in range(phi/2)]
    return B

def gen_gram_matrix(f, algorithm = "fast"):
    phi = euler_phi(f)
    K = CyclotomicField(f)
    z = K.gen()
    B = [z**i for i in range(f)]
    if algorithm == "naive":
        G = matrix(ZZ, phi, phi)
        for i in range(phi):
            G[i,i] = phi
            for j in range(i+1,phi):
                G[i,j] = (B[i]*B[j].conjugate()).trace()
                G[j,i] = G[i,j]
    if algorithm == "fast":
        v = vector(ZZ, f)
        v[0] = phi
        for j in range(1,f):
            v[j] = (B[0]*B[j].conjugate()).trace()
        GExt = matrix.circulant(v)
        G = GExt[:phi,:phi]
    return G

def gen_twisted_gram_matrix(B, twist):
    l = len(B)
    G = matrix([[((twist * B[i]) * (twist * B[j]).conjugate()).trace() for j in range(l)] for i in range(l)])
    return G  

def trace_map(u, v, twist=1, normalise = True):
    if isinstance(u, Vector):
        assert isinstance(v, Vector) 
        l = len(u)
        assert len(v) == l
        result = (u*twist).inner_product( (v*twist).conjugate() ).trace()
        if normalise:
            result = result / (2 * u.base_ring().conductor())
    else:
        assert not isinstance(v, Vector) 
        result = ( (u*twist) * (v*twist).conjugate() ).trace()
        if normalise:
            result = result / (2 * u.parent().conductor())
    return result

################
# Example Code #
################

# # For prime conductors, the real prefix basis is orthogonal under the twisted trace map with twist 1-z.

# f = 11
# assert is_prime(f)
# K.<z> = CyclotomicField(f)
# a = 1-z
# Bp = real_prefix_basis(f)
# G = twisted_trace_map(Bp, a)
# print(G)

# # For prime conductors, it seems that all elements of the real prefix basis are of canonical norm (over the entire cyclotomic field) at most f/2. 
# # Furthermore, the bound seems to get tighter as f increases. 

# results = []
# for f in primes(3,100):
#     print(f)
#     phi = euler_phi(f)
#     K.<z> = CyclotomicField(f)
#     a = 1-z
#     Bp = real_prefix_basis(f)
#     n = max([(b * b.conjugate()).trace().sqrt().n() for b in Bp])
#     results += [[f,n]]
# for f,n in results:
#     print(f - 2*n)
    
