# Circulant representation of cyclotomic elements

from sage.all import *
from sage.structure.element import Vector

def circulant_rep(f, a, q = None, scale = False, output_z_coeff = False):
    """Lift a cyclotomic element to a circulant representation in ``Q[x]/(x^f - 1)``.

    Writes ``x^f - 1 = Phi_f * Psi_f`` with ``Phi_f`` the ``f``-th cyclotomic
    polynomial. Using the CRT idempotent ``Psi_f * (Psi_f^{-1} mod Phi_f)``
    (which is ``1 mod Phi_f`` and ``0 mod Psi_f``), it maps ``a`` to the
    representative ``ahat`` in ``Q[x]/(x^f - 1)`` that agrees with ``a`` modulo
    ``Phi_f`` and vanishes modulo ``Psi_f``. The length-``f`` coefficient vector
    of ``ahat`` is the circulant (cyclic-convolution) representation of ``a``.

    Args:
        f: conductor of the cyclotomic field ``K = Q(zeta_f)``.
        a: an element of ``K``, or a Sage vector of such elements (handled
            entry by entry).
        q: if given (must be prime), the lift is computed over ``GF(q)`` and
            then read back into ``Q[x]``.
        scale: if ``True``, multiply the result by ``f`` to clear the
            denominators introduced by the modular inverse.
        output_z_coeff: if ``True``, return the length-``f`` coefficient list
            (zero-padded) instead of the polynomial; for a vector input the
            per-entry coefficient lists are flattened.

    Returns:
        The polynomial ``ahat`` (default), or its zero-padded length-``f``
        coefficient list when ``output_z_coeff`` is set; a vector/flattened
        list when ``a`` is a vector.
    """
    if isinstance(a, Vector):
        if output_z_coeff:
            return flatten([circulant_rep(f, entry, q, scale, output_z_coeff) for entry in a])
        else:
            return vector([circulant_rep(f, entry, q, scale, output_z_coeff) for entry in a])
    else:
        K = a.parent()
        var('x')
        P = PolynomialRing(QQ, 'x')
        Phi_f = P(cyclotomic_polynomial(f))
        Psi_f = P((x**f - 1)/Phi_f)
        conversion_factor = Psi_f * Psi_f.inverse_mod(Phi_f)
        if q:
            assert is_prime(q)
            var('xq')
            Pq = PolynomialRing(GF(q), 'xq')
            Phi_f_q = Pq(cyclotomic_polynomial(f))
            Psi_f_q = Pq(xq**f - 1)/Phi_f_q
            conversion_factor_mod_q = Pq(Psi_f) * Pq(Psi_f).inverse_mod(Pq(Phi_f))
            ahat = P((Pq(a.polynomial()) * conversion_factor_mod_q).mod(Pq(xq**f - 1)))
        else:
            ahat = (P(a.polynomial()) * conversion_factor).mod(P(x**f - 1))
        if scale:
            ahat = f * ahat
        if output_z_coeff:
            return ahat.list() + [0 for _ in range(f-ahat.degree()-1)]
        else:
            return ahat
