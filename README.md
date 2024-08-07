## How to use

- Clone `lattice_lib` to `<project>/lattice_lib`.
- In `<project>`, create a file called lattice_lib.py with the following content (for example):
```python
import sys
sys.path.append('lattice_lib')
from util import *                  # General utility functions
from conductor_analyser import *    # Tool for analysing conductors
from geometric_norms import *       # Tools related to geometric norms
from subtractive_set import *       # Tools related to subtractive sets
from twisted_trace import *         # Tools related to twisted trace
from circulant_rep import *         # Tools related to circulant representations of cyclotomic elements
```
- Use `from lattice_lib import *` to import the library.
