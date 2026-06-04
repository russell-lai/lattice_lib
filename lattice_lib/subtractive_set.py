# Tools related to subtractive sets

from sage.all import *
from .util import *
var('z')

def smallest_norm(f, print_ideals = False, verbose = False, method = 'multiplicative_order'):
    """Smallest ideal norm in the ``f``-th cyclotomic field.

    Args:
        f: conductor (must be greater than 1).
        print_ideals: if ``True``, also return the ideals attaining the
            smallest norm. Only honoured by the ``"zeta"`` method.
        verbose: if ``True``, print intermediate information.
        method: selects the algorithm. ``"multiplicative_order"`` (default)
            finds the smallest prime power ``p^k`` where ``k`` is the
            multiplicative order of ``p`` modulo the part of ``f`` coprime to
            ``p``. ``"zeta"`` reads off the first nonzero higher coefficient of
            the Dedekind zeta function and is slow for ``f > 200``.

    Returns:
        The smallest ideal norm ``r``, or (with ``print_ideals`` under the
        ``"zeta"`` method) the pair ``(r, factors)`` where ``factors`` lists
        ``[ideal, multiplicity]`` entries.
    """

    assert f > 1

    if method == 'zeta':
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

    if method == 'multiplicative_order':
        p = 1
        current_low = Infinity
        while p <= current_low:
            p = next_prime(p)
            candidate = p**multiplicative_order(mod(p,ZZ(f).prime_to_m_part(p)))
            if candidate < current_low:
                current_low = candidate
        return current_low

def subtractive_set(f, non_zero = False):
    """Construct a subtractive set over the ``f``-th cyclotomic ring.

    Dispatches to :func:`subtractive_set_prime_power` for prime-power
    conductors and to :func:`subtractive_set_non_prime_power` otherwise. A
    subtractive set is one whose pairwise differences of distinct elements are
    all units. Use :func:`is_subtractive` to check this for the returned set.

    Args:
        f: conductor.
        non_zero: forwarded to the prime-power construction to drop the zero
            element (ignored for non-prime-power conductors).
    """
    if is_prime_power(f):
        return subtractive_set_prime_power(f, non_zero)
    else:
        return subtractive_set_non_prime_power(f)

def subtractive_set_prime_power(f, non_zero = False):
    """Subtractive set for a prime-power conductor ``f``.

    Returns the geometric sums ``(z^i - 1)/(z - 1)`` of powers of
    ``z = zeta_f`` for ``i`` in ``range(euler_phi(radical(f)))``. With
    ``non_zero=True`` the ``i = 0`` term (which equals 0) is omitted, so ``i``
    ranges over ``1 .. euler_phi(radical(f)) - 1``.
    """
    f_max = max_prime_power_divisor(f)
    K = CyclotomicField(f)
    z = K.gen()
    if non_zero:
        C = [(z**i-1)/(z-1) for i in range(1,euler_phi(radical(f)))]
    else:
        C = [(z**i-1)/(z-1) for i in range(euler_phi(radical(f)))]
    return C

def subtractive_set_non_prime_power(f):
    """Subtractive set for a non-prime-power conductor ``f``.

    Returns the first ``f / f_max`` powers ``z^i`` of ``z = zeta_f``, where
    ``f_max`` is the largest prime-power divisor of ``f``.
    """
    f_max = max_prime_power_divisor(f)
    K = CyclotomicField(f)
    z = K.gen()
    C = [z**i for i in range(f/f_max)]
    return C

def subtractive_set_expansion_factor(C, t, trials=50):
    """Estimate the expansion factor of a subtractive set by random sampling.

    For ``trials`` random size-``t`` subsets ``S = {s_0, ..., s_{t-1}}`` of
    ``C``, solves the Vandermonde-type system ``V w = (1, 0, ..., 0)`` where
    ``V[j, i] = s_i^j``, and records the largest l2-norm of ``w`` taken over
    the complex embeddings of its entries. The returned value is a sampled
    lower bound on the worst-case expansion factor, not an exact value.

    Args:
        C: a subtractive set.
        t: recovery threshold (the subset size).
        trials: number of random subsets to try (default 50).
    """
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
    """Heuristic integrality test: ``True`` if ``str(s)`` contains no ``"/"``.

    A fast textual check used in place of a ring-membership test. It can
    misjudge elements whose printed form does not match this pattern, so it is
    a heuristic, as the name indicates.
    """
    return (str(s)).find("/") == -1

def is_unit_hack(s):
    """Heuristic unit test: ``True`` if ``1/s`` passes :func:`is_integral_hack`.

    Treats ``s`` as a unit when its inverse looks integral. Shares the
    heuristic caveats of :func:`is_integral_hack`.
    """
    return is_integral_hack(1/s)

def is_subtractive(C, hack=True):
    """Decide whether ``C`` is a subtractive set over the ring of integers.

    Checks that every pairwise difference of distinct elements of ``C`` is a
    unit. With ``hack=True`` (default) it uses the fast textual
    :func:`is_unit_hack`. With ``hack=False`` it uses the exact
    ring-of-integers unit test, which is slower but reliable.

    Args:
        C: a list of elements of a number field.
        hack: select the fast heuristic (``True``) or the exact test
            (``False``).
    """
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
