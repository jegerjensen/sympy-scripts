from sympy import Symbol
from sympy.physics.secondquant import substitute_dummies, _substitute
q = Symbol('q', dummy=True)
i = Symbol('i', below_fermi=True, dummy=True)
a = Symbol('a', above_fermi=True, dummy=True)

reverse = lambda x: reversed(x)
_substitute(a, [q], reverse)   # will succeed
_substitute(a, [q], reverse, above_fermi=True)   # will succeed
_substitute(i, [q], reverse, above_fermi=True)   # will not succeed
_substitute(i, [q], reverse, above_fermi=False)  # will succeed

