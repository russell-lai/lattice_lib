# Test suite for lattice_lib. Run with:  sage tests/test_lattice_lib.py
#
# Plain-assert smoke / invariant tests covering each module, plus regression
# tests for the packaging fixes (conductor_analyser imports, package import).

import os
import sys

# Locate the repo root (the dir containing the lattice_lib package). Sage's file
# loader sets __file__ to its own all.py, so search cwd and its ancestors.
def _find_root():
    cand = os.getcwd()
    for _ in range(6):
        if os.path.exists(os.path.join(cand, 'lattice_lib', 'subtractive_set.py')):
            return cand
        parent = os.path.dirname(cand)
        if parent == cand:
            break
        cand = parent
    raise RuntimeError("could not locate the lattice_lib repo root from cwd; "
                       "run as: sage tests/test_lattice_lib.py")

_root = _find_root()
sys.path.insert(0, _root)   # so `import lattice_lib` / `from lattice_lib import *` resolve

from sage.all import *
from lattice_lib import *


def test_util():
    assert next_power_two(5) == 8
    assert next_power_two(8) == 8
    assert next_power_two(1) == 2
    assert set(prime_power_divisors(12)) == {4, 3}
    assert max_prime_power_divisor(12) == 4
    # BalanceInterval is the centered representative range of Z_q.
    assert BalanceInterval(8) == (-3, 4)
    assert BalanceInterval(7) == (-3, 3)
    assert is_signed_perm([1, -2, 3], [2, 3, -1])
    assert not is_signed_perm([1, 2], [1, 3])
    assert is_local_signed_perm([1, -2, 3, 4], [2, 1, 4, -3], 2)
    assert not is_local_signed_perm([1, 2, 3, 4], [1, 3, 2, 4], 2)


def test_geometric_norms():
    K = CyclotomicField(5)
    z = K.gen()
    # canon_norm(1) = sqrt(Tr(1)) = sqrt(phi).
    assert abs(canon_norm(K(1)) - sqrt(euler_phi(5)).n()) < 1e-9
    # A root of unity has canonical 2-norm sqrt(phi) and coeff inf-norm 1.
    assert abs(canon_norm(z) - sqrt(euler_phi(5)).n()) < 1e-9
    assert abs(coeff_norm(z, oo) - 1) < 1e-9
    # Vector norm aggregates entrywise.
    v = vector(K, [K(1), z])
    assert abs(coeff_norm(v, oo) - 1) < 1e-9


def test_circulant_rep():
    # ahat must reduce to a modulo Phi_f: the circulant rep is a lift to
    # Z[x]/(x^f - 1) of the element of Z[x]/Phi_f.
    for f in [5, 8, 12, 15]:
        K = CyclotomicField(f)
        z = K.gen()
        a = 1 + 2 * z + z**2
        ahat = circulant_rep(f, a)
        P = PolynomialRing(QQ, 'x')
        Phi = P(cyclotomic_polynomial(f))
        assert (P(ahat) - P(a.polynomial())) % Phi == 0


def test_twisted_trace_orthogonality():
    # Documented claim: for a prime conductor, the real prefix basis is
    # orthogonal under the twisted trace map with twist 1 - z.
    for f in [5, 7, 11]:
        K = CyclotomicField(f)
        z = K.gen()
        Bp = real_prefix_basis(f)
        G = gen_twisted_gram_matrix(Bp, 1 - z)
        n = G.nrows()
        for i in range(n):
            for j in range(n):
                if i != j:
                    assert G[i, j] == 0, (f, i, j, G[i, j])


def test_subtractive_sets_are_subtractive():
    for f in [7, 11, 8, 9, 16, 12, 15]:
        C = subtractive_set(f)
        assert len(C) >= 1
        assert is_subtractive(C), f
    # non_zero option drops the zero element (i = 0 term) for prime powers.
    C_all = subtractive_set(9)
    C_nz = subtractive_set(9, non_zero=True)
    assert len(C_nz) == len(C_all) - 1


def test_smallest_norm_methods_agree():
    for f in [5, 7, 8, 11, 12, 15]:
        a = smallest_norm(f, method='multiplicative_order')
        b = smallest_norm(f, method='zeta')
        assert a == b, (f, a, b)


def test_expansion_factor_runs():
    C = subtractive_set(11)
    beta = subtractive_set_expansion_factor(C, t=3, trials=5)
    assert beta > 0


def test_conductor_analyser_regression():
    # Regression: analyse() previously raised NameError on prime_power_divisors.
    row = analyse(12)
    assert row[0] == 12
    assert row[4] == euler_phi(12)
    # A range of conductor types should all return without error.
    for f in [7, 8, 12, 15, 16, 45]:
        analyse(f)


def test_package_import():
    import lattice_lib
    assert hasattr(lattice_lib, 'subtractive_set')
    assert hasattr(lattice_lib, 'analyse')
    assert hasattr(lattice_lib, 'next_power_two')


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    failures = 0
    for t in tests:
        try:
            t()
            print("PASS  {}".format(t.__name__))
        except Exception as e:
            failures += 1
            print("FAIL  {}: {!r}".format(t.__name__, e))
    print("\n{}/{} passed".format(len(tests) - failures, len(tests)))
    sys.exit(1 if failures else 0)


# Run unconditionally: Sage's file loader does not set __name__ to '__main__'.
main()
