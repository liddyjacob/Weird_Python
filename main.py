from __future__ import print_function

import sympy
from tree import *
from sympy import sieve

tree = PAOtree()

test2 = set(sympy.ntheory.generate.primerange(5,10))

print(test2)

tree.display()
tree.add_power()
print(tree.divisors)
tree.display()
tree.add_power()
print(tree.divisors)
n = tree.display()
print(n)
tree.add_power()
print(tree.divisors)
tree.display()

