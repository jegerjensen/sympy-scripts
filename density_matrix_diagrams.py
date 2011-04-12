from sympy.physics.secondquant import *
from sympy import (
    symbols, expand, pprint, Number, latex, Function, preview, Symbol, Tuple
)
from sympy.utilities.codegen import codegen

from utilities.indexify import indexify
from utilities.undummy import undummy
from sympy import Eq, IndexedBase, Idx

P = PermutationOperator

def _report(expr):
    # return str(expr)
    return latex(expr, mode='inline')

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

my_dummies={}
my_dummies['below'] = "ijklm"
my_dummies['above'] = "abcde"

print "Setting up creation/annihilation operators"
# p = Symbol('p', below_fermi=True)
# p = Symbol('p', above_fermi=True)
p = Symbol('p')
q = Symbol('q')

C = Commutator
P = PermutationOperator
T1,T2 = get_CC_operators()
T = T1+ T2

T1,T2 = get_CC_operators()
TT = T1+ T2
pqbar = NO(Fd(p)*F(q))
pqbar = pqbar + Commutator(pqbar, T) + Commutator(Commutator(pqbar, T), TT)/2
pqbar = wicks(pqbar)
pqbar = evaluate_deltas(pqbar)
pqbar = substitute_dummies(pqbar, new_indices=True, pretty_indices=my_dummies)
print latex(pqbar, mode='dmath')


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
l_a = AntiSymmetricTensor('l', tuple(), (a3,))
l_abi = AntiSymmetricTensor('l', (i3,), (a3, a4))
L_Ap1 = l_a*F(a3) + l_abi*NO(Fd(i3)*F(a4)*F(a3))/2

r_i = AntiSymmetricTensor('r',tuple([]),(i1,))
r_aij = AntiSymmetricTensor('r',(a1,),(i1,i2))
R_Am1 = r_i*F(i1) + r_aij*NO(Fd(a1)*F(i2)*F(i1))/2
l_i = AntiSymmetricTensor('l',(i3,),tuple([]))
l_aij = AntiSymmetricTensor('l',(i3,i4),(a3,))
L_Am1 = l_i*Fd(i3) + l_aij*NO(Fd(i3)*Fd(i4)*F(a3))/2


l0,r0 = symbols('l0 r0')
R_A=r0 + AntiSymmetricTensor('r',(a1,),(i1,))*NO(Fd(a1)*F(i1)) + AntiSymmetricTensor('r',(a1,a2),(i1,i2))*NO(Fd(a1)*Fd(a2)*F(i2)*F(i1))/4
L_A=l0 + AntiSymmetricTensor('l',(i3,),(a3,))*NO(Fd(i3)*F(a3)) + AntiSymmetricTensor('l',(i3,i4),(a3,a4))*NO(Fd(i3)*Fd(i4)*F(a4)*F(a3))/4

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

def get_routines(expr, description="", substitutions={}):
    project_descr = 'overlaps'

    below_orbs = Symbol('below_orbs', integer=True)
    total_orbs = Symbol('total_orbs', integer=True)

    expr = undummy(expr)

    expr = expr.subs(substitutions)
    expr = evaluate_deltas(expr)

    lhs_indices = []
    for i in Tuple(p, q).subs(substitutions):
        if i.assumptions0.get('above_fermi'):
            ranges = (below_orbs, total_orbs-1)
        elif i.assumptions0.get('below_fermi'):
            ranges = (0, below_orbs-1)
        label = Symbol(i.name, integer=True)
        lhs_indices.append(Idx(label, ranges))

    lhs = IndexedBase('lhs')[lhs_indices]
    expr = indexify(Eq(lhs, expr))

    print "indexified expr is", expr

    return ('mtransition_'+ description, expr)



# routines = []
# routines.extend(get_routines(overlaps, "_Ap1A"))
# codegen(routines, 'f95', 'code/sympy_overlaps', to_files=True, project='overlaps')

print
print "Transition matrix elements"
print


TMatrix = generate_expressions(L_Am1*pqbar*R_Am1, my_dummies)

print latex(TMatrix, mode='dmath')

pi = Symbol('p', below_fermi=True, integer=True)
qj = Symbol('q', below_fermi=True, integer=True)
pa = Symbol('p', above_fermi=True, integer=True)
qb = Symbol('q', above_fermi=True, integer=True)

routines = []
routines.append(get_routines(TMatrix, substitutions={p:pi, q:qj}, description="ij"))
routines.append(get_routines(TMatrix, substitutions={p:pa, q:qb}, description="ab"))
routines.append(get_routines(TMatrix, substitutions={p:pi, q:qb}, description="ia"))
routines.append(get_routines(TMatrix, substitutions={p:pa, q:qj}, description="ai"))
codegen(routines, 'f95', 'code/sympy_transition_density', to_files=True, project='transition_density')
