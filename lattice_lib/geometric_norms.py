from sage.all import *
from sage.structure.element import Vector, Matrix

def z_coeff(v):
    """Integer/rational coefficient (power-basis) embedding of ``v``.

    For a Sage vector over a cyclotomic field, concatenates the power-basis
    coefficient lists of the entries into one long flat vector. Otherwise
    returns ``vector(v)`` (treating ``v`` as a plain coefficient iterable).
    Inverse of :func:`K_vec`.
    """
    if isinstance(v, Vector):
        return vector(flatten([x.list() for x in v]))
    else:
        return vector(v)

def K_vec(f, z_coeff):
    """Rebuild a vector of cyclotomic elements from a flat coefficient vector.

    Args:
        f: conductor of the cyclotomic field ``K = Q(zeta_f)``.
        z_coeff: a flat coefficient vector whose length is a multiple of
            ``phi(f)``.

    Returns:
        The vector over ``K`` obtained by grouping ``z_coeff`` into consecutive
        blocks of length ``phi(f)`` and reading each block as the power-basis
        coefficients of one element of ``K``. Inverse of :func:`z_coeff`.
    """
    phi = euler_phi(f)
    assert len(z_coeff) % phi == 0
    K = CyclotomicField(f)
    return vector([K(z_coeff[i*phi:(i+1)*phi]) for i in range(len(z_coeff)/phi)])

def coeff_norm(v, order=oo):
    """Norm of the coefficient (power-basis) embedding of ``v``.

    Args:
        v: a cyclotomic-field element, or a vector or matrix over such a field.
        order: the order of the norm to apply, defaulting to ``oo`` (the
            infinity / max-coefficient norm).

    Returns:
        A numerical value. For a single element this is the ``order``-norm of
        its power-basis coefficient vector. A matrix is flattened to a vector,
        and a vector aggregates the per-entry coefficient norms with another
        ``order``-norm.
    """
    if isinstance(v, Matrix):
        return coeff_norm(vector(v), order)
    elif isinstance(v, Vector):
        return vector([coeff_norm(x, order) for x in v]).norm(order)
    else:
        return vector(v).norm(order).n()

def canon_norm(v, order=2):
    """Norm of the canonical (complex-embedding) image of ``v``.

    Args:
        v: a cyclotomic-field element, or a vector or matrix over such a field.
        order: the order of the norm to apply, defaulting to ``2``.

    Returns:
        A numerical value. For a single element with ``order == 2`` this is the
        canonical 2-norm ``sqrt(Tr(v * conj(v)))``, equal to the l2-norm of the
        complex embeddings; for other orders it is the ``order``-norm of the
        complex embeddings. A matrix is flattened to a vector, and a vector
        aggregates the per-entry canonical norms with another ``order``-norm.
    """
    if isinstance(v, Matrix):
        return canon_norm(vector(v), order)
    elif isinstance(v, Vector):
        return vector([canon_norm(x, order) for x in v]).norm(order)
    else:
        if order == 2:
            return (v * v.conjugate()).trace().sqrt().n()
        else:
            return vector(v.complex_embeddings()).norm(order).n()
