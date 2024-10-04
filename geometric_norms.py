from sage.all import *
from sage.structure.element import Vector, Matrix

def z_coeff(v):
    if isinstance(v, Vector):
        return vector(flatten([x.list() for x in v]))
    else: 
        return vector(v)

def K_vec(f, z_coeff):
    phi = euler_phi(f)
    assert len(z_coeff) % phi == 0
    K = CyclotomicField(f)
    return vector([K(z_coeff[i*phi:(i+1)*phi]) for i in range(len(z_coeff)/phi)])
    
def coeff_norm(v, order=oo):
    if isinstance(v, Matrix):
        return coeff_norm(vector(v), order)
    elif isinstance(v, Vector):
        return vector([coeff_norm(x, order) for x in v]).norm(order)
    else: 
        return vector(v).norm(order).n()

def canon_norm(v, order=2):
    if isinstance(v, Matrix):
        return canon_norm(vector(v), order)
    elif isinstance(v, Vector):
        return vector([canon_norm(x, order) for x in v]).norm(order)
    else: 
        if order == 2:
            return (v * v.conjugate()).trace().sqrt().n()
        else:
            return vector(v.complex_embeddings()).norm(order).n()
