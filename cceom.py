from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol, Eq
)

from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import ExprSplitter

from utilities.indexify import indexify
from utilities.undummy import undummy
from utilities.ccm import get_CC_operators

p,q,r,s,t,u = symbols('pqrstu',dummy=True)

Hbar = ( AntiSymmetricTensor('h',(p,),(q,))*NO(Fd(p)*F(q)) +
		AntiSymmetricTensor('h',(p,q),(r,s))*
		NO(Fd(p)*Fd(q)*F(s)*F(r))/4 +
		AntiSymmetricTensor('h',(p,q,r),(s,t,u))*
		NO(Fd(p)*Fd(q)*Fd(r)*F(u)*F(t)*F(s))/36
		)


C = Commutator
P = PermutationOperator


print
print "Evaluating operator:"
print latex(Hbar, mode='inline')
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

"""

print
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
val_aij = substitute_dummies(val_aij,new_indices=True)
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
val_aij = substitute_dummies(val_aij,new_indices=True)
val_aij = simplify_index_permutations(val_aij,[P(a1,a2),P(i1,i2)])
print latex(val_aij)
# preview(val_aij,output="dvi")


"""


print
print "Evaluating Hbar on R_{A}"
print

r_ai = AntiSymmetricTensor('r',(c,),(k,))
r_abij = AntiSymmetricTensor('r',(c,d),(k,l))
R = r_ai*Fd(c)*F(k) + r_abij*NO(Fd(c)*Fd(d)*F(l)*F(k))/4

print
print "r0"
print
val_0 = wicks(C(Hbar, R),keep_only_fully_contracted=True)
val_0 = evaluate_deltas(val_0)
val_0 = substitute_dummies(val_0, new_indices=True)
print latex(val_0, mode='inline')

print
print "r_ai"
print
val_ai = wicks(Fd(i1)*F(a1)*C(Hbar, R),keep_only_fully_contracted=True)
val_ai = evaluate_deltas(val_ai)
val_ai = substitute_dummies(val_ai, new_indices=True)
print latex(val_ai, mode='inline')
lhs1 = AntiSymmetricTensor('lhs', (a1,), (i1,))




print
print "r_abij"
print
val_abij = wicks(NO(Fd(i1)*Fd(i2)*F(a2)*F(a1))*C(Hbar, R),keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
val_abij = substitute_dummies(val_abij,new_indices=True)
val_abij = simplify_index_permutations(val_abij,[P(i1,i2), P(a1,a2)])
print latex(val_abij, mode='inline')
lhs2 = AntiSymmetricTensor('lhs', (a1,a2), (i1,i2))

stop

expr_ai = indexify(undummy(Eq(lhs1, val_ai)))
print expr_ai
expr_abij = indexify(undummy(Eq(-lhs2, -val_abij)))
print expr_abij

spl = ExprSplitter('diagram_RA')
routines = spl.spawn_routines(expr_ai)
routines.extend(spl.spawn_routines(expr_abij))
codegen(routines, 'f95', 'code/eomcc/sympy_eomcc', to_files=True)

stop


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
val_abij = substitute_dummies(val_abij,new_indices=True)
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
val_abci = substitute_dummies(val_abci,new_indices=True)
val_abci = simplify_index_permutations(val_abci,[P(i1,i2), P(a1,a2)])
print latex(val_abci)
# preview(val_abci,output="dvi")
