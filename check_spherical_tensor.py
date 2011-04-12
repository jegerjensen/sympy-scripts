from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol
)

def get_CC_operators():
	"""
	Returns a tuple (T1,T2) of unique operators.
	"""
	i = symbols('i',below_fermi=True,dummy=True)
	a = symbols('a',above_fermi=True,dummy=True)
	t_ai = AntiSymmetricTensor('t',(a,),(i,))
	ai = NO(Fd(a)*F(i))
	i,j = symbols('ij',below_fermi=True,dummy=True)
	a,b = symbols('ab',above_fermi=True,dummy=True)
	t_abij = AntiSymmetricTensor('t',(a,b),(i,j))
	abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

	T1 = t_ai*ai
	T2 = Number((1,4))*t_abij*abji
	return (T1,T2)

print "Setting up J2 operator"
k,l = symbols('kl',dummy=True,below_fermi=True)
c,d = symbols('cd',dummy=True,above_fermi=True)
p,q = symbols('pq',dummy=True)


Jabove = AntiSymmetricTensor('J',(d,),(c,))
Jbelow = AntiSymmetricTensor('J',(l,),(k,))
Jop = Jabove*NO(Fd(d)*F(c)) + Jbelow*NO(Fd(l)*F(k))


C = Commutator
P = PermutationOperator


# [J,T]
T1,T2 = get_CC_operators()
T = T1+ T2


comm = wicks(C(Jop,T))
comm = evaluate_deltas(comm)
comm = substitute_dummies(comm,reverse_order=False, new_indices=True)



print
print "*********************"
print "Evaluating comutator:"
print latex(C(Jop,T))
print
print "Result"
print latex(comm)
comm.subs([(l,k+1),(c,c+1)])
print latex(comm)
preview(comm,output='png')

kjkjn

# new symbols to construct matrix elements
i1= Symbol('i',below_fermi=True)
i2= Symbol('j',below_fermi=True)
a1= Symbol('a',above_fermi=True)
a2= Symbol('b',above_fermi=True)


# Now, when all J-operators are setup, we define new dummy symbols to use with R
k,l,m = symbols('klm',dummy=True,below_fermi=True)
c,d = symbols('cd',dummy=True,above_fermi=True)



print
print "*********************"
print "Evaluating Jop on R_{A-1}"

r_i = AntiSymmetricTensor('r',tuple([]),(k,))
r_aij = AntiSymmetricTensor('r',(c,),(k,l))
R = r_i*F(k) + r_aij*NO(Fd(c)*F(l)*F(k))/2

print
print "r_i"
val_i = Wicks(Fd(i1)*JopBar*R,keepOnlyFullyContracted=True)
val_i = evaluate_deltas(val_i)
val_i = substitute_dummies(val_i, newIndices=True)
print latex(val_i)
preview(val_i,output='png')

print
print "r_aij"
val_aij = Wicks(NO(Fd(i1)*Fd(i2)*F(a1))*JopBar*R,keepOnlyFullyContracted=True, simplifyKroneckerDeltas=True)
val_aij = substitute_dummies(val_aij,newIndices=True, reverse_order=True)
val_aij = simplifyIndexPermutations(val_aij,[P(i1,i2)])
print latex(val_aij)
preview(val_aij,output='png')




print
print "*********************"
print "Evaluating Jop on L_{A-1}"

l_i = AntiSymmetricTensor('l',(k,),tuple([]))
l_aij = AntiSymmetricTensor('l',(k,l),(c,))
L = l_i*Fd(k) + l_aij*NO(Fd(k)*Fd(l)*F(c))/2


print
print "l_i"
val_i = Wicks(L*JopBar*F(i1),keepOnlyFullyContracted=True)
val_i = evaluate_deltas(val_i)
val_i = substitute_dummies(val_i, newIndices=True)
print latex(val_i)
preview(val_i,output='png')

print
print "l_aij"
val_aij = Wicks(L*JopBar*NO(Fd(a1)*F(i2)*F(i1)),keepOnlyFullyContracted=True, simplifyKroneckerDeltas=True)
val_aij = substitute_dummies(val_aij,newIndices=True, reverse_order=True)
val_aij = simplifyIndexPermutations(val_aij,[P(a1,a2),P(i1,i2)])
print latex(val_aij)
preview(val_aij,output='png')
