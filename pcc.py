from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Add
)
from utilities.ccm import get_CC_operators

class pypar_dummy:
    def size(self):
        return 1
    def rank(self):
        return 0

pypar = pypar_dummy()

numprocs = pypar.size()
iam = pypar.rank()
master = 0

def parallel_Wicks(expr):
    result = []
    for ind in xrange(len(expr.args)):
        if ind%numprocs == iam:
            wicked = Wicks(expr.args[ind])
            if not wicked:
                continue
            # print "raw length:",len(wicked.args)
            wicked = evaluate_deltas(wicked)
            # print "after eval_delta:",len(wicked.args)
            wicked = substitute_dummies(wicked)
            # print "after subst_dummies:",len(wicked.args)
            result.append(wicked)
    result = Add(*result)

    my_len = len(result.args)

    # pypar.reduce(my_len, sum, )
    return result




print "Setting up hamiltonian"
p,q,r,s = symbols('pqrs',dummy=True)
f = SymmetricTensor('f',p,q)
pr = NO((Fd(p)*F(q)))
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))

H=f*pr +Number(1,4)*v*pqsr

# simplify H immediately
H=H.doit()
H=H.expand()
print len(H.args)
H=evaluate_deltas(H)
print len(H.args)
H=substitute_dummies(H)
print len(H.args)


print "Calculating nested commutators"
C = Commutator

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm1..."
comm1 = parallel_Wicks(C(H,T))

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm2..."
comm2 = parallel_Wicks(C(comm1,T))

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm3..."
comm3 = parallel_Wicks(C(comm2,T))

T1,T2 = get_CC_operators()
T = T1+ T2
print "comm4..."
comm4 = parallel_Wicks(C(comm3,T))

print "construct Hausdoff expansion..."
eq = H + comm1+comm2/2+comm3/6+comm4/24
eq = eq.expand()


# import pickle
# filename='testing.pkl'
# f1 = open(filename,'w')
# pickle.dump(eq,f1)
# f1.close()
# f2 = open(filename,'r')
# eq2 = pickle.load(f2)
# f2.close()
# eq = eq2


print "simplifying"
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

