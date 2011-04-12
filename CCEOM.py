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
	t_ai = SymmetricTensor('t',a,i)
	ai = NO(Fd(a)*F(i))
	i,j = symbols('ij',below_fermi=True,dummy=True)
	a,b = symbols('ab',above_fermi=True,dummy=True)
	t_abij = AntiSymmetricTensor('t',(a,b),(i,j))
	abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

	T1 = t_ai*ai
	T2 = Number((1,4))*t_abij*abji
	return (T1,T2)

p,q,r,s,t,u = symbols('pqrstu',dummy=True)

Hbar = ( AntiSymmetricTensor('H',(p,),(q,))*NO(Fd(p)*F(q)) +
		AntiSymmetricTensor('H',(p,q),(r,s))*
		NO(Fd(p)*Fd(q)*F(s)*F(r))/4 +
		AntiSymmetricTensor('H',(p,q,r),(s,t,u))*
		NO(Fd(p)*Fd(q)*Fd(r)*F(u)*F(t)*F(s))/36
		)


C = Commutator
P = PermutationOperator


print
print "*********************"
print "Evaluating operator:"
print latex(Hbar)
print


# New symbols to construct matrix elements
i1= Symbol('i',below_fermi=True)
i2= Symbol('j',below_fermi=True)
a1= Symbol('a',above_fermi=True)
a2= Symbol('b',above_fermi=True)
a3= Symbol('c',above_fermi=True)


# New dummy symbols to use with R
k,l,m = symbols('klm',dummy=True,below_fermi=True)
c,d,e = symbols('def',dummy=True,above_fermi=True)



print
print "*********************"
print "Evaluating Hbar on R_{A-1}"

r_i = AntiSymmetricTensor('r',tuple([]),(k,))
r_aij = AntiSymmetricTensor('r',(c,),(k,l))
R = r_i*F(k) + r_aij*NO(Fd(c)*F(l)*F(k))/2


print
print "r_i"
val_i = wicks(Fd(i1)*Hbar*R,keep_only_fully_contracted=True)
val_i = evaluate_deltas(val_i)
val_i = substitute_dummies(val_i, new_indices=True)
print latex(val_i)
# preview(val_i,output="dvi")

print
print "r_aij"
val_aij = wicks(NO(Fd(i1)*Fd(i2)*F(a1))*Hbar*R,keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_aij = substitute_dummies(val_aij,new_indices=True, reverse_order=True)
val_aij = simplify_index_permutations(val_aij,[P(i1,i2)])
print latex(val_aij)
# preview(val_aij,output="dvi")




print
print "*********************"
print "Evaluating Hbar on L_{A-1}"

l_i = AntiSymmetricTensor('l',(k,),tuple([]))
l_aij = AntiSymmetricTensor('l',(k,l),(c,))
L = l_i*Fd(k) + l_aij*NO(Fd(k)*Fd(l)*F(c))/2


print
print "l_i"
val_i = wicks(L*Hbar*F(i1),keep_only_fully_contracted=True)
val_i = evaluate_deltas(val_i)
val_i = substitute_dummies(val_i, new_indices=True)
print latex(val_i)
# preview(val_i,output="dvi")

print
print "l_aij"
val_aij = wicks(L*Hbar*NO(Fd(a1)*F(i2)*F(i1)),keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_aij = substitute_dummies(val_aij,new_indices=True, reverse_order=True)
val_aij = simplify_index_permutations(val_aij,[P(a1,a2),P(i1,i2)])
print latex(val_aij)
# preview(val_aij,output="dvi")




print
print "*********************"
print "Evaluating Hbar on R_{A}"

r_ai = AntiSymmetricTensor('r',(c,),(k,))
r_abij = AntiSymmetricTensor('r',(c,d),(k,l))
R = r_ai*Fd(c)*F(k) + r_abij*NO(Fd(c)*Fd(d)*F(l)*F(k))/4

print
print "r_ai"
val_ai = wicks(Fd(i1)*F(a1)*Hbar*R,keep_only_fully_contracted=True)
val_ai = evaluate_deltas(val_ai)
val_ai = substitute_dummies(val_ai, new_indices=True)
print latex(val_ai)
# preview(val_ai,output="dvi")

print
print "r_abij"
val_abij = wicks(NO(Fd(i1)*Fd(i2)*F(a2)*F(a1))*Hbar*R,keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_abij = substitute_dummies(val_abij,new_indices=True, reverse_order=True)
val_abij = simplify_index_permutations(val_abij,[P(i1,i2), P(a1,a2)])
print latex(val_abij)
# preview(val_abij,output="dvi")




print
print "*********************"
print "Evaluating Hbar on L_{A}"

l_ai = AntiSymmetricTensor('l',(k,),(c,))
l_abij = AntiSymmetricTensor('l',(k,l),(c,d))
L = l_ai*NO(Fd(k)*F(c)) + l_abij*NO(Fd(k)*Fd(l)*F(d)*F(c))/4


print
print "l_ai"
val_ai = wicks(L*Hbar*Fd(a1)*F(i1),keep_only_fully_contracted=True)
val_ai = evaluate_deltas(val_ai)
val_ai = substitute_dummies(val_ai, new_indices=True)
print latex(val_ai)
# preview(val_ai,output="dvi")

print
print "l_abij"
val_abij = wicks(L*Hbar*NO(Fd(a1)*Fd(a2)*F(i2)*F(i1)),keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_abij = substitute_dummies(val_abij,new_indices=True, reverse_order=True)
val_abij = simplify_index_permutations(val_abij,[P(a1,a2),P(i1,i2)])
print latex(val_abij)
# preview(val_abij,output="dvi")



print
print "*********************"
print "Evaluating Hbar on R_{A+2}"

r_ab = AntiSymmetricTensor('r',(c,),(d,))
r_abci = AntiSymmetricTensor('r',(c,d,e),(k,))
R = r_ab*Fd(c)*Fd(d)/2 + r_abci*NO(Fd(c)*Fd(d)*Fd(e)*F(k))/6

print
print "r_ab"
val_ab = wicks(F(a2)*F(a1)*Hbar*R,keep_only_fully_contracted=True)
val_ab = evaluate_deltas(val_ab)
val_ab = substitute_dummies(val_ab, new_indices=True)
print latex(val_ab)
# preview(val_ab,output="dvi")

print
print "r_abci"
val_abci = wicks(NO(Fd(i1)*F(a3)*F(a2)*F(a1))*Hbar*R,keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_abci = substitute_dummies(val_abci,new_indices=True, reverse_order=True)
val_abci = simplify_index_permutations(val_abci,[P(i1,i2), P(a1,a2)])
print latex(val_abci)
# preview(val_abci,output="dvi")
