from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol
)

from utilities.ccm import get_CC_operators

print "Setting up J2 operator"
k,l = symbols('kl',dummy=True,below_fermi=True)
c,d = symbols('cd',dummy=True,above_fermi=True)


Jabove = AntiSymmetricTensor('J',(d,),(c,))
Jbelow = AntiSymmetricTensor('J',(l,),(k,))
Jop = Jabove*NO(Fd(d)*F(c)) + Jbelow*NO(Fd(l)*F(k))


C = Commutator
P = PermutationOperator


# [J,T]
T1,T2 = get_CC_operators()
T = T1+ T2
comm1 = Wicks(C(Jop,T))
comm1 = evaluate_deltas(comm1)

# [[J,T],T]
T1,T2 = get_CC_operators()
T = T1+ T2
comm2 = Wicks(C(comm1,T))
comm2 = evaluate_deltas(comm2)

# \bar J
JopBar = Jop + comm1 + comm2/2
JopBar = JopBar.expand()
JopBar = substitute_dummies(JopBar,reverse_order=False, newIndices=True)



print
print "*********************"
print "Evaluating operator:"
print latex(Jop)
print
print latex(JopBar)
# preview(JopBar,output='png')


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
