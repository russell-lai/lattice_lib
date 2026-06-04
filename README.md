# lattice_lib

SageMath tools for lattice-based cryptography over cyclotomic rings.

## Requirements

SageMath (developed and tested against Sage 10.9). All modules assume `from sage.all import *`.

## How to use

The importable package lives in the `lattice_lib/` subdirectory of this repository. There are two ways to use it.

**Without installing** (vendoring as a git submodule or clone): add the repository directory, not the package directory, to `sys.path` and import.

```python
import sys
sys.path.append('path/to/lattice-lib')   # the repo root, the parent of the package
from lattice_lib import *
```

**Installing** into your Sage environment:

```
sage -pip install -e path/to/lattice-lib
```

after which `from lattice_lib import *` works from anywhere.

Either way, `from lattice_lib import *` re-exports the public names of every module into one namespace:

- `util`: general utility functions
- `conductor_analyser`: tools for analysing conductors
- `geometric_norms`: tools related to geometric norms
- `subtractive_set`: tools related to subtractive sets
- `twisted_trace`: tools related to the twisted trace
- `circulant_rep`: circulant representations of cyclotomic elements

## API overview

Every public function carries a docstring, so `help(name)` works in Sage. The main entry points by module:

**util**
- `next_power_two(x)`: smallest power of two that is at least `x`.
- `prime_power_divisors(f)`, `max_prime_power_divisor(f)`: maximal prime-power divisors and their maximum.
- `is_signed_perm(A, B)`, `is_local_signed_perm(A, B, chunk_size)`: signed-permutation tests.
- `BalanceInterval(q)`: centred residue interval for `Z_q`.

**geometric_norms**
- `z_coeff(v)`, `K_vec(f, z_coeff)`: convert between flat coefficient vectors and vectors of cyclotomic elements.
- `coeff_norm(v, order=oo)`, `canon_norm(v, order=2)`: coefficient-embedding and canonical-embedding norms.

**circulant_rep**
- `circulant_rep(f, a, ...)`: circulant representation of cyclotomic elements in `Q[x]/(x^f - 1)`.

**twisted_trace**
- `power_basis(f)`, `real_power_basis(f)`, `real_prefix_basis(f)`: bases of the cyclotomic field and its maximal real subfield.
- `gen_gram_matrix(f, ...)`, `gen_twisted_gram_matrix(B, twist)`, `trace_map(u, v, ...)`: trace and twisted-trace Gram matrices and inner product.

**subtractive_set**
- `subtractive_set(f, ...)`: construct a subtractive set over the `f`-th cyclotomic ring.
- `smallest_norm(f, ...)`: smallest ideal norm of the cyclotomic field.
- `subtractive_set_expansion_factor(C, t, ...)`: sampled expansion-factor estimate.
- `is_subtractive(C, ...)`: check the subtractive property.

**conductor_analyser**
- `analyse(f)`: summary row for a single conductor.
- `gen_conductor_table(bound=...)`: write `conductors.csv` over a range of conductors.

## Tests

Run the test suite from the repository root:

```
sage tests/test_lattice_lib.py
```

It prints a PASS/FAIL line per test and exits non-zero on any failure.

## License

This library is licensed under the [LGPLv3+](https://www.gnu.org/licenses/lgpl-3.0.en.html). See [LICENSE](LICENSE) for the notice, and `COPYING` / `COPYING.LESSER` for the full GPLv3 and LGPLv3 texts.
