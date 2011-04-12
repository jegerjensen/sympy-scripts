from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex
)
import cPickle as pickle

from utilities.ccm import get_CC_operators

class io_obj:
    filename='private/pickled_cc4body.pkl'
    def open(self, mode='r'):
        try:
            self.file = open(self.filename,mode)
        except IOError:
            return False
        return True

    def close(self):
        self.file.close()

    def save(self,descript,obj):
        self.open('a')
        pickle.dump((descript,obj),self.file)
        self.close()
    def load(self):
        descript, obj = pickle.load(self.file)
        return descript, obj
io = io_obj()


print "Setting up hamiltonian"
p,q,r,s,t,u,v,w = symbols('pqrstuvw',dummy=True)
f = SymmetricTensor('f',p,q)
pr = NO((Fd(p)*F(q)))
v2 = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
v3= AntiSymmetricTensor('v',(p,q,r),(s,t,u))
pqruts = NO(Fd(p)*Fd(q)*Fd(r)*F(u)*F(t)*F(s))
v4= AntiSymmetricTensor('v',(p,q,r,s),(t,u,v,w))
pqrswvut = NO(Fd(p)*Fd(q)*Fd(r)*Fd(s)*F(w)*F(v)*F(u)*F(t))

# simplify H immediately
H=f*pr +v2*pqsr/4 + v3*pqruts/36 + v4*pqrswvut/576
print "Expanding normally ordered strings in H..."
H=H.doit()
H=H.expand()
print len(H.args)," terms in H"
H=evaluate_deltas(H)
print len(H.args)," terms in H"
H=substitute_dummies(H)
print len(H.args)," terms in H"

print "Calculating nested commutators"
C = Commutator


# load pickled commutators,calculate the others
commlist = ('comm1', 'comm2', 'comm3', 'comm4', 'comm5', 'comm6', 'comm7', 'comm8')
exec("=".join(commlist)+"=0")
if io.open():
    still_more=True
else:
    still_more=False
for i in range(len(commlist)):
    descript = commlist[i]
    if still_more:
        try:
            descr, obj = io.load()
            print "loading ", descr
            # replace dummy indices to get dummy-count correct
            if descr == descript:
                obj = substitute_dummies(obj,newIndices=True)
                exec(descr+'=obj')
            continue

        except EOFError:
            print " got exception at ",descript
            still_more=False
            io.close()

    if not i:
        T1,T2 = get_CC_operators()
        T = T1+ T2
        print descript, "With H"
        comm1 = Wicks(C(H,T),simplifyDummies=True, simplifyKroneckerDeltas=True)
        io.save(descript,comm1)
    else:
        descr_last = commlist[i-1]
        T1,T2 = get_CC_operators()
        T = T1+ T2
        print "doing ",descript
        execstring=(descript+"=Wicks(C("+descr_last+",T),simplifyDummies=True, simplifyKroneckerDeltas=True)")
        exec(descript+"=Wicks(C("+descr_last+",T),simplifyDummies=True, simplifyKroneckerDeltas=True)")
        io.save(descript, eval(descript))






print "construct Hausdoff expansion..."
eq = H + comm1+comm2/2+comm3/6+comm4/24+comm5/120+comm6/720+comm7/5040+comm8/40320
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

