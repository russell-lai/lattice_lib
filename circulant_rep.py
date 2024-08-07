# Circulant representation of cyclotomic elements

from sage.all import *
from sage.structure.element import Vector

def circulant_rep(f, a, q = None, scale = False, output_z_coeff = False):
    # If scale = True, the output is scaled by f, making it integral.
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