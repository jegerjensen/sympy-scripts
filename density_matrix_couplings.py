from sympy.physics.racahalgebra import (
        refine_phases, refine_tjs2sjs, convert_cgc2tjs, is_equivalent,
        SixJSymbol, ASigma, combine_ASigmas, evaluate_sums, apply_deltas,
        apply_orthogonality, ClebschGordanCoefficient, extract_symbol2dummy_dict,
        extract_dummy2symbol_dict, SphericalTensor, refine_tjs2njs
        )
from sympy.physics.braket import (
        MatrixElement, ReducedMatrixElement, apply_wigner_eckardt,
        rewrite_as_direct_product, DirectMatrixElement,
        SphericalTensorOperator, Dagger, SphFermKet, SphFermBra,
        FermKet, FermBra, ThreeTensorMatrixElement,
        inject_every_symbol_globally, rewrite_coupling,
        braket_assumptions, DualSphericalTensorOperator
        )
from sympy import (
        Symbol, symbols, global_assumptions, Assume, ask, Q, Mul, S, sqrt,
        pprint, Eq, pprint_use_unicode, latex, preview, fcode, Function, Dij
        )
import sys

pprint_use_unicode(False)
def _report(expr):

    if isinstance(expr, basestring):
        print r"\begin{verbatim}"
        print expr
        print r"\end{verbatim}"
        return

    pprint(str(expr))
    # print fcode(expr)
    # print(expr)
    # pprint(expr)

    # subsdict = extract_dummy2symbol_dict(expr)
    # expr = expr.subs(subsdict)
    # print latex(expr, mode='dmath', mul_symbol='dot')
    print

i, j, k, l = map(SphFermKet, 'ijkl')
a, b, c, d = map(SphFermKet, 'abcd')

# EOM operators

J_RAm1 = Symbol('J_RAm1', nonnegative=True)
M_RAm1 = Symbol('M_RAm1')
J_LAm1 = Symbol('J_LAm1', nonnegative=True)
M_LAm1 = Symbol('M_LAm1')
braket_assumptions.add(Assume(J_RAm1, 'half_integer'))
braket_assumptions.add(Assume(M_RAm1, 'half_integer'))
braket_assumptions.add(Assume(J_LAm1, 'half_integer'))
braket_assumptions.add(Assume(M_LAm1, 'half_integer'))
LAm1 = DualSphericalTensorOperator('L', J_LAm1, M_LAm1)
RAm1 = SphericalTensorOperator('R', J_RAm1, M_RAm1)
T = SphericalTensorOperator('T', S.Zero, S.Zero)

# pq operators

J_pq, j_p, j_q = symbols('J_pq j_p j_q', nonnegative=True)
M_pq, m_p, m_q = symbols('M_pq m_p m_q')
braket_assumptions.add(Assume(j_p, Q.half_integer))
braket_assumptions.add(Assume(m_p, Q.half_integer))
braket_assumptions.add(Assume(j_q, Q.half_integer))
braket_assumptions.add(Assume(m_q, Q.half_integer))
braket_assumptions.add(Assume(J_pq, Q.integer))
braket_assumptions.add(Assume(M_pq, Q.integer))
pd = SphericalTensorOperator('pd', j_p, m_p)
q = DualSphericalTensorOperator('q', j_q, m_q)
PdQ = SphericalTensor('PdQ', J_pq, M_pq, pd, q)

print
print "Defined tensor operators:"
_report(LAm1)
_report(RAm1)
_report(T)
_report(pd)
_report(q)
_report(PdQ)

# Matrix elements

print
print "*************** $<i| L | >$ *****************"
print

l_i = MatrixElement(Dagger(i), LAm1, "")
l_i_sph = ReducedMatrixElement(Dagger(i), LAm1, "")
_report( Eq(l_i, l_i_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(l_i_sph.get_direct_product_ito_self(tjs=1)))

print
print "*************** $<ij| L | a>$ *****************"
print

l_aij = MatrixElement((Dagger(i), Dagger(j)), LAm1, a)
l_aij_sph = ReducedMatrixElement(SphFermBra('ij', Dagger(i), Dagger(j)), LAm1, a)
_report(Eq(l_aij, l_aij_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(l_aij_sph.get_direct_product_ito_self(tjs=1)))

print
print "*************** $<| R | i>$ *****************"
print

r_i = MatrixElement(0, RAm1, i)
r_i_sph = ReducedMatrixElement(0, RAm1, i)
_report( Eq(r_i, r_i_sph.get_direct_product_ito_self(tjs=0)))
_report( Eq(r_i, r_i_sph.get_direct_product_ito_self(tjs=1)))

print
_report(fcode(r_i_sph.get_direct_product_ito_self(tjs=1)))
print
print "*************** $<a| R | ij>$ *****************"
print

r_aij = MatrixElement(Dagger(a), RAm1, (i, j))
r_aij_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ij', i, j, reverse=0))
_report( Eq(r_aij, r_aij_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(r_aij_sph.get_direct_product_ito_self(tjs=1)))

print
print "*************** $<a| T | i>$ *****************"
print

t_ai = MatrixElement(Dagger(a), T, i)
t_ai_sph = ReducedMatrixElement(Dagger(a), T, i)
_report(Eq(t_ai, t_ai_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(t_ai_sph.get_direct_product_ito_self(tjs=1)))

print
print "*************** $<ab| T | ij>$ *****************"
print

t_abij = MatrixElement((Dagger(a), Dagger(b)), T, (i, j))
t_abij_sph = ReducedMatrixElement(SphFermBra('ab', Dagger(a), Dagger(b), reverse=0), T, SphFermKet('ij', i, j, reverse=0))
_report(Eq(t_abij, t_abij_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(t_abij_sph.get_direct_product_ito_self(tjs=1)))


#
# Symbols
#

J_ij, J_Am1, j_a, j_b, j_i, j_j, j_k = symbols('J_ij J_Am1 j_a j_b j_i j_j j_k', nonnegative=True)
M_ij, M_Am1, m_a, m_b, m_i, m_j, m_k = symbols('M_ij M_Am1 m_a m_b m_i m_j m_k')

coupling_pq, directpq = PdQ.as_coeff_direct_product(use_dummies=0)
reduction = combine_ASigmas(coupling_pq*ASigma(M_RAm1, M_pq)
        *ClebschGordanCoefficient(J_RAm1, M_RAm1, J_pq, M_pq, J_LAm1, M_LAm1))

print
print "We are going to derive coupled expressions for the quantity:"
_report(reduction*Symbol('\\mathrm{mscheme\_diagram}'))

print
print "pq-coupling can be written:"
coupling_pq = convert_cgc2tjs(coupling_pq)
coupling_pq = refine_phases(coupling_pq)
_report(coupling_pq)
_report(fcode(coupling_pq))

print
print "geometrical reduction can be written:"
red_expr = ClebschGordanCoefficient(J_RAm1, M_RAm1, J_pq, M_pq, J_LAm1, M_LAm1)
red_expr = convert_cgc2tjs(red_expr)
red_expr = refine_phases(red_expr)
_report(red_expr)
_report(fcode(red_expr))

"""
print
print "Diagram $-l^{ij}_{a} r^{}_{i} t^{a}_{k} t^{b}_{j} \delta_{pk} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_b, m_q:m_b
        }

t_ak = t_ai.subs(i, k)
t_ak_sph = t_ai_sph.subs(i, k)
t_bj = t_ai.subs({Dagger(a):Dagger(b), i:j})
t_bj_sph = t_ai_sph.subs({Dagger(a):Dagger(b), i:j})

expr_msch = -l_aij*r_i*t_ak*t_bj*reduction.subs(pq_bindings)*ASigma(m_a, m_i, m_j)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_i_sph, l_aij_sph, t_ak_sph, t_bj_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [m_i, m_j])
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, verbose=True)
_report(expr_sph)


print
print "Diagram $-\\frac{1}{2} l^{ij}_{a} r^{b}_{ij} t^{a}_{k} \delta_{pk} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_b, m_q:m_b
        }

r_bij = r_aij.subs(Dagger(a), Dagger(b))
r_bij_sph = r_aij_sph.subs(Dagger(a), Dagger(b))
t_ak = t_ai.subs(i, k)
t_ak_sph = t_ai_sph.subs(i, k)

expr_msch = -l_aij*r_bij*t_ak*reduction.subs(pq_bindings)*ASigma(m_a, m_i, m_j)/2
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_bij_sph, l_aij_sph, t_ak_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [m_i, m_j])
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, verbose=True)
_report(expr_sph)

print
print "Diagram $ \\frac{1}{2} l^{ij}_{a} r^{}_{k} t^{ab}_{ij} \delta_{pk} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_b, m_q:m_b
        }

r_k = r_i.subs(i, k)
r_k_sph = r_i_sph.subs(i, k)

expr_msch = l_aij*r_k*t_abij*reduction.subs(pq_bindings)*ASigma(m_a, m_i, m_j)/2
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_k_sph, l_aij_sph, t_abij_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [m_i, m_j])
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [M_pq, M_RAm1, m_k])
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
d = extract_symbol2dummy_dict(expr_sph)
M_ab = d[Symbol('M_ab')]
expr_sph = apply_orthogonality(expr_sph, [m_a, M_ab])
_report(expr_sph)
expr_sph = refine_phases(expr_sph, verbose=True)
_report(expr_sph)


print
print "Diagram $ l^{ij}_{a} r^{}_{i} t^{b}_{j} \delta_{pa} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_a, m_p:m_a,
        j_q:j_b, m_q:m_b
        }

t_bj = t_ai.subs({Dagger(a):Dagger(b), i:j})
t_bj_sph = t_ai_sph.subs({Dagger(a):Dagger(b), i:j})

expr_msch = l_aij*r_i*t_bj*reduction.subs(pq_bindings)*ASigma(m_i, m_j)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_i_sph, l_aij_sph, t_bj_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, verbose=True)
_report(expr_sph)


print
print "Diagram $-l^{i}_{} r^{}_{j} t^{a}_{i} \delta_{pj} \delta_{qa}$"
print

pq_bindings = {
        j_p:j_j, m_p:m_j,
        j_q:j_a, m_q:m_a
        }

r_j = r_i.subs(i, j)
r_j_sph = r_i_sph.subs(i, j)

expr_msch = -l_i*r_j*t_ai*reduction.subs(pq_bindings)*ASigma(m_i)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_j_sph, l_i_sph, t_ai_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [M_pq, m_j])
_report(expr_sph)




print
print "diagram $-l^{ij}_{a} r^{}_{i} t^{a}_{k} \delta_{pk} \delta_{qj}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_j, m_q:m_j
        }

t_ak = t_ai.subs(i, k)
t_ak_sph = t_ai_sph.subs(i, k)

expr_msch = -l_aij*r_i*t_ak*reduction.subs(pq_bindings)*ASigma(m_a, m_i)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_i_sph, l_aij_sph, t_ak_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, verbose=True)
_report(expr_sph)

print
print "diagram $ l^{ij}_{a} r^{}_{i} \delta_{pa} \delta_{qj}$"
print

pq_bindings = {
        j_p:j_a, m_p:m_a,
        j_q:j_j, m_q:m_j
        }

expr_msch = l_aij*r_i*reduction.subs(pq_bindings)*ASigma(m_i)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_i_sph, l_aij_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)


print
print "diagram $l^{i}_{} r^{a}_{ij} \delta_{pj} \delta_{qa}$"
print

pq_bindings = {
        j_p:j_j, m_p:m_j,
        j_q:j_a, m_q:m_a
        }

expr_msch = l_i*r_aij*reduction.subs(pq_bindings)*ASigma(m_i)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_aij_sph, l_i_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)

_report(fcode(expr_sph))

print
print "diagram $- l^{i}_{} r^{}_{j} \delta_{pj} \delta_{qi} $"
print

pq_bindings = {
        j_p:j_j, m_p:m_j,
        j_q:j_i, m_q:m_i
        }

r_j = r_i.subs(i, j)
r_j_sph = r_i_sph.subs(i, j)

expr_msch = -r_j*l_i*reduction.subs(pq_bindings)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_j_sph, l_i_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [M_pq, m_j])
_report(expr_sph)

print
print "diagram $\\frac{1}{2} l^{ij}_{a} r^{b}_{ij} \delta_{pa} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_a, m_p:m_a,
        j_q:j_b, m_q:m_b
        }

r_bij = r_aij.subs(Dagger(a), Dagger(b))
r_bij_sph = r_aij_sph.subs(Dagger(a), Dagger(b))

expr_msch = r_bij*l_aij/2*ASigma(m_i, m_j)*reduction.subs(pq_bindings)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_bij_sph, l_aij_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = apply_orthogonality(expr_sph, [m_i, m_j])
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, verbose=True)
_report(expr_sph)
_report(fcode(expr_sph))
"""

print
print "diagram $+l^{i}_{} r^{}_{i} t^{a}_{j} \delta_{pj} \delta_{qa}$"
print "and"
print "diagram $\\frac{1}{2} l^{ij}_{a} r^{a}_{ij} t^{b}_{k} \delta_{pk} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_i, m_p:m_i,
        j_q:j_a, m_q:m_a
        }

expr_msch = t_ai*reduction.subs(pq_bindings)
_report(expr_msch)
expr_msch = expr_msch*Dij(M_RAm1, M_LAm1)*Dij(J_RAm1, J_LAm1)*Dij(M_pq, 0)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [t_ai_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_phases(expr_sph)
_report(expr_sph)

print "rewrite to simpler form hinders use of orthogonality relation"

_report(fcode(expr_sph))
stop

print
print "diagram $-l^{ij}_{a} r^{a}_{ik} \delta_{pk} \delta_{qj}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_j, m_q:m_j
        }
r_aik = r_aij.subs(j, k)
r_aik_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ik', i, k))

expr_msch = -l_aij*r_aik*reduction.subs(pq_bindings)*ASigma(m_a, m_i)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_aik_sph, l_aij_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_phases(expr_sph)
_report(expr_sph)



print "Trick: insert $1 =\\sum(M_ij, M_ik)(J_{ij} M_{ij} J_{ik} M_{ik}|J_{pq} M_{pq})^2$"

dummies = extract_symbol2dummy_dict(expr_sph)
jij = dummies[Symbol('J_ij', nonnegative=True)]
mij = dummies[Symbol('M_ij')]
jik = dummies[Symbol('J_ik', nonnegative=True)]
mik = dummies[Symbol('M_ik')]

expr_sph = expr_sph*ASigma(mij, mik)\
        *ClebschGordanCoefficient(jij,-mij,jik,mik,J_pq,M_pq)**2
_report(expr_sph)
expr_sph = refine_phases(convert_cgc2tjs(expr_sph))
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)


print
print "diagram $-l^{ij}_{a} r^{}_{i} t^{ab}_{kj}\delta_{pk}\delta_{qb}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_b, m_q:m_b
        }

t_abkj = t_abij.subs(i, k)
t_abkj_sph = ReducedMatrixElement(SphFermBra('ab', Dagger(a), Dagger(b)), T, SphFermKet('kj', k, j))

expr_msch = -l_aij*r_i*t_abkj*reduction.subs(pq_bindings)*ASigma(m_a, m_i, m_j)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_i_sph, l_aij_sph, t_abkj_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_phases(expr_sph, verbose=True)
_report(expr_sph)

print "Trick: insert $1 =\\sum(m_a,m_j)(j_am_aj_j{-m_j}|J_{pq}M{pq})^2$"
expr_sph = expr_sph*ASigma(m_j, m_a)\
        *ClebschGordanCoefficient(j_a,m_a,j_j,-m_j,J_pq,M_pq)**2
_report(expr_sph)
expr_sph = refine_phases(convert_cgc2tjs(expr_sph))
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, let_pass=0)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)



print
print "Diagram $-l^{ij}_{a} r^{a}_{ik} t^{b}_{j} \delta_{pk} \delta_{qb}$"
print

pq_bindings = {
        j_p:j_k, m_p:m_k,
        j_q:j_b, m_q:m_b
        }

t_bj = t_ai.subs({Dagger(a):Dagger(b), i:j})
t_bj_sph = t_ai_sph.subs({Dagger(a):Dagger(b), i:j})
r_aik = r_aij.subs(j, k)
r_aik_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ik', i, k))

expr_msch = -l_aij*r_aik*t_bj*reduction.subs(pq_bindings)*ASigma(m_a, m_i, m_j)
_report(expr_msch)
expr_sph = rewrite_coupling(expr_msch, [r_aik_sph, l_aij_sph, t_bj_sph])
expr_sph = combine_ASigmas(expr_sph)
_report(expr_sph)
expr_sph = convert_cgc2tjs(expr_sph)
_report(expr_sph)
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(expr_sph)
expr_sph = refine_phases(expr_sph, verbose=True)
_report(expr_sph)


print "Trick: insert $1 =\\sum(M_ij, M_ik)(J_{ij} M_{ij} J_{ik} M_{ik}|J_{pq} M_{pq})^2$"

dummies = extract_symbol2dummy_dict(expr_sph)
jij = dummies[Symbol('J_ij', nonnegative=True)]
mij = dummies[Symbol('M_ij')]
jik = dummies[Symbol('J_ik', nonnegative=True)]
mik = dummies[Symbol('M_ik')]

expr_sph = expr_sph*ASigma(mij, mik)\
        *ClebschGordanCoefficient(jij,-mij,jik,mik,J_pq,M_pq)**2
_report(expr_sph)
expr_sph = refine_phases(convert_cgc2tjs(expr_sph))
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph, let_pass=0)
_report(expr_sph)
expr_sph = refine_tjs2sjs(expr_sph)
_report(expr_sph)

