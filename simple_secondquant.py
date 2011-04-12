from sympy import symbols, Basic, Symbol, Function, Rational

# Define classes for Second quantization operators
#
# To avoid messing with the internals of SymPy, we use the @property decorator
# to overload Basic.is_commutative.  (SymPy prevents you to change .is_commutative
# in the __init__ method.)

class Creator(Basic):

    def __init__(self, label, **kw_args):
        assert isinstance(label, Symbol), "label must be a Symbol instance"
        self.label = label
        #self.is_commutative = False  <--- this is not allowed...

    @property
    def is_commutative(self):#       <---- ...so here is how we do it.
        return False

    def _sympystr(self, p):
        return "C_%s" % self.label

class Annihilator(Basic):

    def __init__(self, label, **kw_args):
        assert isinstance(label, Symbol), "label must be a Symbol instance"
        self.label = label

    @property
    def is_commutative(self):
        return False

    def _sympystr(self, p):
        return "A_%s" % self.label

# Shorter names will be convienient in formulas
Ann = Annihilator
Cre = Creator

# In Sympy, every symbol must be declared
p,q,r,s = symbols('p q r s')

# You can attach assumptions to symbols
a,b = symbols('a b', above_fermi=True)
i,j = symbols('i j', below_fermi=True)

# the assumptions are stored in a dictionary called assumptions0
assert a.assumptions0.get('above_fermi') == True
assert a.assumptions0.get('below_fermi') == None
assert i.assumptions0.get('below_fermi') == True
assert i.assumptions0.get('above_fermi') == None

# Note that symbols that differ only in assumptions are still considered
# equal in the mathematical sense
assert Symbol('i', is_positive=True) == Symbol('i', is_positive=False)
assert Symbol('i', is_commutative=True) == Symbol('i', is_commutative=False)
assert Symbol('i', above_fermi=True) == Symbol('i', below_fermi=True)

# ... but they are different as python objects
assert hash(Symbol('i')) != hash(Symbol('i', above_fermi=True))
# for comparison:
assert hash(Symbol('i')) == hash(Symbol('i'))

# Let's setup a hamiltonian
#
# Terms in the hamiltonian can be represented by general function symbols
t = Function('t')
f = Function('f')
v = Function('v')

# 1/4 must be written as Rational(1,4) to avoid integer trunkation
H = f(p,q)*Cre(p)*Ann(q) + Rational(1, 4)*v(p, q, r, s)*Cre(p)*Cre(q)*Ann(s)*Ann(r)

print "Hamiltonian is: %s" % H

# setup a T1 Coupled Cluster operator
T1 = t(a, i)*Cre(a)*Ann(i)

# We need a simple implementation of a commutator
def commutator(A, B):
    expr = A*B - B*A
    # call expand() to return an expression without parenthesis
    return expr.expand()

expr = commutator(H, T1)
print "Commutator expression is: %s" %expr
print

# the terms in the expression are available in the tuple expr.args
for i,term in enumerate(expr.args):
    # the factors in each term are available in the tuple term.args
    print "Term %i has factors %s" % (i, ", ".join(map(str, term.args)))
    for factor in term.args:
        print "   Factor %s is of type %s" % (str(factor), type(factor))
    print

# Now we just need an implementation of wicks theorem, e.g.
def wicks(expr):

    # A clever implementation goes here

    return expr

evaluated = wicks(expr)
