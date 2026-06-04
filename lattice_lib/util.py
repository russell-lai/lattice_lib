# General utility functions

from sage.all import *
_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)

def next_power_two(x):
    """Smallest power of two that is at least ``x``.

    Args:
        x: an integer.

    Returns:
        The least ``y = 2**k`` with ``y >= x``. Always at least 2, so
        ``next_power_two(1) == 2`` and ``next_power_two(2) == 2``.
    """
    y = 2
    while y < x:
        y *= 2
    return y

def prime_power_divisors(f):
    """Maximal prime-power divisors of ``f``.

    Args:
        f: a positive integer.

    Returns:
        One entry ``p**e`` per prime ``p`` dividing ``f``, where ``p**e`` is
        the exact power of ``p`` in the factorisation of ``f``. For example
        ``prime_power_divisors(12)`` returns ``[4, 3]``.
    """
    factors = factor(f)
    return [factors[i][0]**factors[i][1] for i in range(len(factors))]

def max_prime_power_divisor(f):
    """Largest maximal prime-power divisor of ``f`` (its power-smoothness)."""
    return max(prime_power_divisors(f))

def is_signed_perm(A,B):
    """Decide whether ``A`` and ``B`` are signed permutations of each other.

    Returns ``True`` iff one list is obtained from the other by permuting the
    entries and flipping signs, i.e. the sorted multisets of absolute values
    agree. Requires ``len(A) == len(B)``.
    """
    assert len(A) == len(B)
    A_abs = [abs(a) for a in A]
    B_abs = [abs(b) for b in B]
    A_abs.sort()
    B_abs.sort()
    return A_abs == B_abs

def is_local_signed_perm(A,B,chunk_size):
    """Decide whether ``B`` is a block-wise signed permutation of ``A``.

    Splits each list into consecutive chunks of length ``chunk_size`` and
    returns ``True`` iff ``B`` is obtained from ``A`` by permuting the chunks
    and then applying a signed permutation within each chunk. Requires
    ``len(A) == len(B)`` and ``chunk_size`` to divide ``len(A)``.
    """
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

def BalanceInterval(q):
    """Centred representative interval ``(lo, hi)`` for ``Z_q``.

    Returns the inclusive bounds of the balanced residue range modulo ``q``:
    ``(-q/2 + 1, q/2)`` when ``q`` is even, and ``(ceil(-q/2), floor(q/2))``
    when ``q`` is odd. For example ``BalanceInterval(8) == (-3, 4)`` and
    ``BalanceInterval(7) == (-3, 3)``.
    """
    if is_even(q):
        return -q/_sage_const_2 +_sage_const_1 , q/_sage_const_2
    else:
        return ceil(-q/_sage_const_2 ), floor(q/_sage_const_2 )
