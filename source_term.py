from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol
)

P = PermutationOperator
C = Commutator

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
H= Number(1,4)*v*pqsr
x = Symbol('x')
H = C(F(x),H)

print "Calculating nested commutators"

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm1..."
comm1 = Wicks(C(H,T),simplifyDummies=True, simplifyKroneckerDeltas=True)

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm2..."
comm2 = Wicks(C(comm1,T),simplifyDummies=True, simplifyKroneckerDeltas=True)

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm3..."
comm3 = Wicks(C(comm2,T),simplifyDummies=True, simplifyKroneckerDeltas=True)

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm4..."
comm4 = Wicks(C(comm3,T),simplifyDummies=True, simplifyKroneckerDeltas=True)

print "construct Hausdoff expansion..."
eq = H + comm1+comm2/2+comm3/6+comm4/24
eq = eq.expand()
eq = evaluate_deltas(eq)
eq = substitute_dummies(eq, newIndices=True, reverse_order=False)
source = eq
print
print "*********************"
print "Evaluating source:"
print latex(source)



print "*********************"
print
print "Setting up Left/Right solutions"





# new symbols to construct matrix elements
i1= Symbol('i',below_fermi=True, dummy=True)
i2= Symbol('j',below_fermi=True, dummy=True)
i3= Symbol('k',below_fermi=True, dummy=True)
i4= Symbol('l',below_fermi=True, dummy=True)
a1= Symbol('a',above_fermi=True, dummy=True)
a2= Symbol('b',above_fermi=True, dummy=True)
a3= Symbol('c',above_fermi=True, dummy=True)
a4= Symbol('d',above_fermi=True, dummy=True)


print

r_i = AntiSymmetricTensor('r',tuple([]),(i1,))
r_aij = AntiSymmetricTensor('r',(a1,),(i1,i2))
R_Am1 = r_i*F(i1) + r_aij*NO(Fd(a1)*F(i2)*F(i1))/2
l_i = AntiSymmetricTensor('l',(i1,),tuple([]))
l_aij = AntiSymmetricTensor('l',(i1,i2),(a1,))
L_Am1 = l_i*Fd(i1) + l_aij*NO(Fd(i1)*Fd(i2)*F(a1))/2

print
print "*********************"
print "A-1 system:"
print latex(R_Am1)

print "*********************"

print
# new symbols again
i1= Symbol('i',below_fermi=True, dummy=True)
i2= Symbol('j',below_fermi=True, dummy=True)
i3= Symbol('k',below_fermi=True, dummy=True)
i4= Symbol('l',below_fermi=True, dummy=True)
a1= Symbol('a',above_fermi=True, dummy=True)
a2= Symbol('b',above_fermi=True, dummy=True)
a3= Symbol('c',above_fermi=True, dummy=True)
a4= Symbol('d',above_fermi=True, dummy=True)

l0,r0 = symbols('lr')

L_A=l0 + AntiSymmetricTensor('l',(i1,),(a1,))*NO(Fd(i1)*F(a1)) + AntiSymmetricTensor('l',(i1,i2),(a1,a2))*NO(Fd(i1)*Fd(i2)*F(a2)*F(a1))/4
R_A=r0 + AntiSymmetricTensor('r',(a3,),(i3,))*NO(Fd(a3)*F(i3)) + AntiSymmetricTensor('r',(a3,a4),(i3,i4))*NO(Fd(a3)*Fd(a4)*F(i4)*F(i3))/4

# eq = L_A*crebar*R_Am1
eq = L_Am1*source*R_A



# preview(eq,output='png')
eq = eq.expand()
eq = substitute_dummies(eq, newIndices=True)
overlaps = Wicks(eq, keepOnlyFullyContracted=True)
overlaps=evaluate_deltas(overlaps)
overlaps=substitute_dummies(overlaps, newIndices=True)
print latex(overlaps)

preview(overlaps,output='png')
