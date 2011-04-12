from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex
)
import cPickle as pickle

def get_CC_operators():
	"""
	Returns a tuple (T1,T2) of unique operators.
	"""
	i = symbols('i',below_fermi=True,dummy=True)
	a = symbols('a',above_fermi=True,dummy=True)
	t_ai = SymmetricTensor('t',a,i)
	ai = NO(Fd(a)*F(i))
	i,j = symbols('ij',below_fermi=True,dummy=True)
	a,b = symbols('ab',above_fermi=True,dummy=True)
	t_abij = AntiSymmetricTensor('t',(a,b),(i,j))
	abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

	T1 = t_ai*ai
	T2 = Number((1,4))*t_abij*abji
	return (T1,T2)

print "Setting up hamiltonian"
p,q,r,s = symbols('pqrs',dummy=True)
f = SymmetricTensor('f',p,q)
pr = NO((Fd(p)*F(q)))
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))

H=f*pr +Number(1,4)*v*pqsr

print "Calculating nested commutators"
C = Commutator

# T1,T2 = get_CC_operators()
# T = T1+ T2
# print "comm1..."
# comm1 = Wicks(C(H,T),simplifyDummies=True, simplifyKroneckerDeltas=True)
# 
# T1,T2 = get_CC_operators()
# T = T1+ T2
# print "comm2..."
# comm2 = Wicks(C(comm1,T),simplifyDummies=True, simplifyKroneckerDeltas=True)
# 
# T1,T2 = get_CC_operators()
# T = T1+ T2
# print "comm3..."
# comm3 = Wicks(C(comm2,T),simplifyDummies=True, simplifyKroneckerDeltas=True)
# 
# T1,T2 = get_CC_operators()
# T = T1+ T2
# print "comm4..."
# comm4 = Wicks(C(comm3,T),simplifyDummies=True, simplifyKroneckerDeltas=True)
# 
# print "construct Hausdoff expansion..."
# eq = H + comm1+comm2/2+comm3/6+comm4/24

filename='private/testing.pkl'
# f1 = open(filename,'w')
# pickle.dump(eq,f1)
# f1.close()
f2 = open(filename,'r')
eq2 = pickle.load(f2)
f2.close()
eq = eq2

eq = eq.expand()
eq = evaluate_deltas(eq)
eq = substitute_dummies(eq, newIndices=True, reverse_order=False)
print "*********************"
print

print "extracting CC equations from full Hbar"
i,j,k,l = symbols('ijkl',below_fermi=True)
a,b,c,d = symbols('abcd',above_fermi=True)
print
print "CC Energy:"
print latex(Wicks(eq, simplifyDummies=True,
    keepOnlyFullyContracted=True))
print
print "CC T1:"
eqT1 = Wicks(NO(Fd(i)*F(a))*eq, simplifyKroneckerDeltas=True, keepOnlyFullyContracted=True)
eqT1 = substitute_dummies(eqT1,reverse_order=False)
print latex(eqT1)
print
print "CC T2:"
eqT2 = Wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*eq,simplifyDummies=True, keepOnlyFullyContracted=True, simplifyKroneckerDeltas=True)
P = PermutationOperator
eqT2 = simplifyIndexPermutations(eqT2,[P(a,b),P(i,j)])
print latex(eqT2)

