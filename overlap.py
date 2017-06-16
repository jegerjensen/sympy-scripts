from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol
)
from sympy.utilities.codegen import codegen
from sympy.utilities.codegen import ExprSplitter

from utilities.indexify import indexify
from utilities.undummy import undummy
from sympy import Eq, IndexedBase, Idx

from utilities.ccm import get_CC_operators
P = PermutationOperator

def _report(expr):
    # return str(expr)
    return latex(expr, mode='inline')

print "Setting up creation/annihilation operators"
# p = Symbol('p', below_fermi=True)
# p = Symbol('p', above_fermi=True)
p = Symbol('p')
q = Symbol('q')

C = Commutator
P = PermutationOperator
T1,T2 = get_CC_operators()
T = T1+ T2

cre = Fd(p)
ann = F(p)
crebar = cre + wicks(C(cre,T))
annbar = ann + wicks(C(ann,T))

print
print "crebar=",_report(crebar)
print
print "annbar=",_report(annbar)
print

# preview(crebar,output='png')
# preview(annbar,output='png')



#  symbols to construct matrix elements
i1= Symbol('i',below_fermi=True, dummy=True)
i2= Symbol('j',below_fermi=True, dummy=True)
i3= Symbol('k',below_fermi=True, dummy=True)
i4= Symbol('l',below_fermi=True, dummy=True)
a1= Symbol('a',above_fermi=True, dummy=True)
a2= Symbol('b',above_fermi=True, dummy=True)
a3= Symbol('c',above_fermi=True, dummy=True)
a4= Symbol('d',above_fermi=True, dummy=True)


r_a = AntiSymmetricTensor('r', (a1,), tuple())
r_abi = AntiSymmetricTensor('r', (a1, a2), (i1,))
R_Ap1 = r_a*Fd(a1) + r_abi*NO(Fd(a1)*Fd(a2)*F(i1))/2
l_a = AntiSymmetricTensor('l', tuple(), (a1,))
l_abi = AntiSymmetricTensor('l', (i1,), (a1, a2))
L_Ap1 = l_a*F(a1) + l_abi*NO(Fd(i1)*F(a2)*F(a1))/2

r_i = AntiSymmetricTensor('r',tuple([]),(i1,))
r_aij = AntiSymmetricTensor('r',(a1,),(i1,i2))
R_Am1 = r_i*F(i1) + r_aij*NO(Fd(a1)*F(i2)*F(i1))/2
l_i = AntiSymmetricTensor('l',(i1,),tuple([]))
l_aij = AntiSymmetricTensor('l',(i1,i2),(a1,))
L_Am1 = l_i*Fd(i1) + l_aij*NO(Fd(i1)*Fd(i2)*F(a1))/2


l0,r0 = symbols('l0 r0')
L_A=l0 + AntiSymmetricTensor('l',(i3,),(a3,))*NO(Fd(i3)*F(a3)) + AntiSymmetricTensor('l',(i3,i4),(a3,a4))*NO(Fd(i3)*Fd(i4)*F(a4)*F(a3))/4
R_A=r0 + AntiSymmetricTensor('r',(a3,),(i3,))*NO(Fd(a3)*F(i3)) + AntiSymmetricTensor('r',(a3,a4),(i3,i4))*NO(Fd(a3)*Fd(a4)*F(i4)*F(i3))/4

print
print "*********************"
print "A system:"
print
print "RA = ",_report(R_A)
print
print "LA = ",_report(L_A)

print
print "A-1 system:"
print
print "RAm1 = ",_report(R_Am1)
print
print "LAm1 = ",_report(L_Am1)
print
print "A+1 system:"
print
print "RAp1 = ",_report(R_Ap1)
print
print "LAp1 = ",_report(L_Ap1)
print
print "*********************"

def generate_expressions(expr, dummies):

    eq = expr.expand()
    overlaps = wicks(eq, keep_only_fully_contracted=True)
    overlaps = evaluate_deltas(overlaps)
    overlaps = substitute_dummies(overlaps, new_indices=True, pretty_indices=dummies)

    return overlaps

def get_routines(expr, description=""):
    project_descr = 'overlaps'

    below_orbs = Symbol('below_orbs', integer=True)
    total_orbs = Symbol('total_orbs', integer=True)

    expr = undummy(expr)


    expr_i = expr.subs(Symbol('p'), Symbol('i', below_fermi=True))
    expr_i = evaluate_deltas(expr_i)
    lhs_i = IndexedBase('lhs')[Idx('i', below_orbs)]
    expr_i = indexify(Eq(lhs_i, expr_i))

    print "expr_i is", expr_i


    expr_a = expr.subs(Symbol('p'), Symbol('a', above_fermi=True))
    expr_a = evaluate_deltas(expr_a)
    lhs_a = IndexedBase('lhs')[Idx('a', (below_orbs, total_orbs-1))]
    expr_a = indexify(Eq(lhs_a, expr_a))
    print "expr_a is", expr_a

    spl = ExprSplitter('overlap_%s_i' %description)
    routines = spl.spawn_routines(expr_i)

    spl = ExprSplitter('overlap_%s_a' %description)
    routines.extend(spl.spawn_routines(expr_a))

    return routines



my_dummies={}
my_dummies['below'] = "ijklm"
my_dummies['above'] = "abcde"
routines = []

# < A | c' | A-1 >
overlaps = generate_expressions(L_A*crebar*R_Am1, my_dummies)
print
print "A a' A-1"
print
print _report(overlaps)

# < A-1 | c  | A >
overlaps = generate_expressions(L_Am1*annbar*R_A, my_dummies)
print
print "A-1 a A"
print
print
print _report(overlaps)

# < A | c  | A+1 >
overlaps = generate_expressions(L_A*annbar*R_Ap1, my_dummies)
print
print "A a A+1"
print
print _report(overlaps)
routines.extend(get_routines(overlaps, "_AAp1"))

# < A+1 | c' | A >
overlaps = generate_expressions(L_Ap1*crebar*R_A, my_dummies)
print
print "A+1 a' A"
print
print _report(overlaps)
routines.extend(get_routines(overlaps, "_Ap1A"))


codegen(routines, 'f95', 'code/overlaps/sympy_overlaps_Ap1', to_files=True, project='overlaps')

