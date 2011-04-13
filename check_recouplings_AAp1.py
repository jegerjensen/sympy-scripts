from sympy.physics.racahalgebra import (
        refine_phases, refine_tjs2sjs, convert_cgc2tjs, is_equivalent,
        SixJSymbol, ASigma, combine_ASigmas, evaluate_sums, apply_deltas,
        apply_orthogonality, ClebschGordanCoefficient, extract_symbol2dummy_dict,
        extract_dummy2symbol_dict, apply_identity_tjs
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
        pprint, Eq, pprint_use_unicode, latex, preview, fcode, Function
        )
import sys

from utilities.redmats import ReducedMatrixElement_left, ReducedMatrixElement_right
from utilities.ccm_sph import *

pprint_use_unicode(False)
def _report(expr):

    if isinstance(expr, basestring):
        print r"\begin{verbatim}"
        print expr
        print r"\end{verbatim}"
        return

    # pprint(str(expr))
    # return

    # print fcode(expr)

    # print(expr)

    # pprint(expr)
    subsdict = extract_dummy2symbol_dict(expr)
    expr = expr.subs(subsdict)
    print latex(expr, mode='dmath', mul_symbol='dot')
    # print latex(expr, mode='inline')

    # preview(expr)

    print

i, j, k, l = map(SphFermKet, 'ijkl')
a, b, c, d = map(SphFermKet, 'abcd')


print
print "Defined tensor operators:"
_report( LA)
_report( RA)
_report( LAp1)
_report( RAp1)
_report( V)
_report( T)

# print
# print "*************** <i| L | a> *****************"
# print

l_ai = MatrixElement(Dagger(i), LA, a)
l_ai_sph = ReducedMatrixElement(Dagger(i), LA, a, definition='left')
# _report( Eq(l_ai, l_ai_sph.get_direct_product_ito_self(tjs=0)))

# print
# _report(fcode(l_ai_sph.get_direct_product_ito_self(tjs=1)))

# print
# print "*************** <ij| L | ab> *****************"
# print

l_abij = MatrixElement((Dagger(i), Dagger(j)), LA, (a, b))
l_abij_sph = ReducedMatrixElement(SphFermBra('ij', Dagger(i), Dagger(j)), LA, SphFermKet('ab', a, b), definition='left')
# _report( Eq(l_abij, l_abij_sph.get_direct_product_ito_self(tjs=0)))

# print
# _report(fcode(l_abij_sph.get_direct_product_ito_self(tjs=1)))

# print
# print "*************** <a| T | i> *****************"
# print

t_ai = MatrixElement(Dagger(a), T, i)
t_ai_sph = ReducedMatrixElement(Dagger(a), T, i)
# _report(Eq(t_ai, t_ai_sph.get_direct_product_ito_self(tjs=0)))

# print
# _report(fcode(t_ai_sph.get_direct_product_ito_self(tjs=1)))

# print
# print "*************** <ab| T | ij> *****************"
# print

t_abij = MatrixElement((Dagger(a), Dagger(b)), T, (i, j))
t_abij_sph = ReducedMatrixElement(SphFermBra('ab', Dagger(a), Dagger(b)), T, SphFermKet('ij', i, j))
# _report(Eq(t_abij, t_abij_sph.get_direct_product_ito_self(tjs=0)))

# print
# _report(fcode(t_abij_sph.get_direct_product_ito_self(tjs=1)))



# print
# print "*************** <| R | a> *****************"
# print

r_a = MatrixElement(Dagger(a), RAp1, 0)
r_a_sph = ReducedMatrixElement(Dagger(a), RAp1, 0, definition='right')
# _report( Eq(r_a, r_a_sph.get_direct_product_ito_self(tjs=0)))
# _report( Eq(r_a, r_a_sph.get_direct_product_ito_self(tjs=1)))

print
_report(fcode(r_a_sph.get_direct_product_ito_self(tjs=1)))
# print
# print "*************** <a| R | ij> *****************"
# print

r_abi = MatrixElement((Dagger(a), Dagger(b)), RAp1, i)
r_abi_sph = ReducedMatrixElement(SphFermBra('ab', Dagger(a), Dagger(b)), RAp1, i, definition='right')
# _report( Eq(r_abi, r_abi_sph.get_direct_product_ito_self(tjs=0)))

print
_report(fcode(r_abi_sph.get_direct_product_ito_self(tjs=1), source_format='free'))




print("\subsection{coupled elements}")
_report(Eq(r_a_sph,rewrite_coupling(r_a_sph, r_a)))
_report(Eq(r_abi_sph,rewrite_coupling(r_abi_sph, r_abi)))
_report(Eq(l_ai_sph,rewrite_coupling(l_ai_sph, l_ai)))
_report(Eq(l_abij_sph,rewrite_coupling(l_abij_sph, l_abij)))
_report(Eq(t_ai_sph,rewrite_coupling(t_ai_sph, t_ai)))
_report(Eq(t_abij_sph,rewrite_coupling(t_abij_sph, t_abij)))


l0 = Symbol('l_0')
r0 = Symbol('r_0')
SF = Symbol('SF')

m_a, m_b, m_c, m_i = symbols('m_a m_b m_c m_i')



print
print "\subsection{Recoupling $<{A}|a|{A+1}>$}"
print

j_a = Symbol('j_a', nonnegative=True)
m_a = Symbol('m_a')
braket_assumptions.add(Assume(j_a, 'half_integer'))
sf_reduction = (-1)**(j_a - m_a)*ClebschGordanCoefficient(J_Ap1, M_Ap1, j_a, -m_a, J_A, M_A)*ASigma(M_Ap1, m_a)


print
print "using reduction"
print fcode(convert_cgc2tjs(sf_reduction))


stop


print
print("\subsubsection{recoupling diagram $\\frac12 l^{ij}_{cb}r^bt_{ij}^{ac} $}")
print

l_cbij = l_abij.subs(a,c)
l_cbij_sph = l_abij_sph.subs(a,c)
r_b = r_a.subs(Dagger(a), Dagger(b))
r_b_sph = r_a_sph.subs(Dagger(a), Dagger(b))
t_acij = t_abij.subs(Dagger(b), Dagger(c))
t_acij_sph = t_abij_sph.subs(Dagger(b), Dagger(c))

expr_msc = combine_ASigmas(S(1)/2*l_cbij*r_b*t_acij*ASigma('m_c','m_b','m_i','m_j') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_cbij_sph, r_b_sph, t_acij_sph])
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_i', 'm_j'])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1, let_pass=0)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))

print
print ("\subsubsection{recoupling diagram $\\frac12l^{ij}_{cb}r^{cb}_it^a_j$}")
print
l_cbij = l_abij.subs(a,c)
l_cbij_sph = l_abij_sph.subs(a,c)
r_cbi = r_abi.subs(Dagger(a), Dagger(c))
r_cbi_sph = r_abi_sph.subs(Dagger(a), Dagger(c))
t_aj = t_ai.subs(i, j)
t_aj_sph = t_ai_sph.subs(i, j)

expr_msc = combine_ASigmas(
        S(1)/2*l_cbij*r_cbi*t_aj*ASigma('m_j','m_i','m_c', 'm_b')
        * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_cbij_sph, r_cbi_sph, t_aj_sph])
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_c', 'm_b'])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))


print
print ("\subsubsection{recoupling diagram $-l^i_br^bt_i^a$}")
print 
l_bi = l_ai.subs(a, b)
l_bi_sph = l_ai_sph.subs(a, b)
r_b = r_a.subs(Dagger(a), Dagger(b))
r_b_sph = r_a_sph.subs(Dagger(a), Dagger(b))


expr_msc = combine_ASigmas(
        t_ai*r_b*l_bi*ASigma('m_b','m_i') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_bi_sph, r_b_sph, t_ai_sph])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, [m_a, M_Ap1])
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))

print
print "\subsubsection{recoupling diagram $l_0 r^a$}"
print
expr_msc = combine_ASigmas(l0*r_a * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [r_a_sph])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = apply_identity_tjs(expr_sph)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))

print
print ("\subsection{recoupling diagram $l^i_br^{ab}_i$}")
print
l_bi = l_ai.subs(a, b)
l_bi_sph = l_ai_sph.subs(a, b)
expr_msc = combine_ASigmas(l_bi*r_abi*ASigma('m_b', 'm_i') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_bi_sph, r_abi_sph])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))
print
print "\subsection{Recoupling <{A}|i|{A+1}>}"
print

j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
braket_assumptions.add(Assume(j_i, 'half_integer'))
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient(J_Ap1, M_Ap1, j_i, -m_i, J_A, M_A)*ASigma(M_Ap1, m_i)

# print "reduction factor"
# _report( sf_reduction)
# _report(fcode(refine_phases(convert_cgc2tjs(sf_reduction))))

print
print("\subsubsection{recoupling diagram $-l_a^i r^a$}")
expr_msc = combine_ASigmas(-ASigma(m_a)*l_ai*r_a*sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_ai_sph, r_a_sph])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, [M_Ap1, m_i])
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))

print
print("\subsubsection{recoupling diagram $-\\frac12 l^{ij}_{ab}r^{ab}_j$}")
print
r_abj = r_abi.subs(i,j)
r_abj_sph = r_abi_sph.subs(i,j)

expr_msc = combine_ASigmas(-r_abj*l_abij*ASigma('m_a','m_b','m_j') * sf_reduction / 2)
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [r_abj_sph, l_abij_sph])
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, [m_a, m_b])
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))


print
print("\subsubsection{simplified diagram $\\frac12 l^{ij}_{cb}r^bt_{ij}^{ac} $}")
print "JA=0"
print

l_cbij = l_abij.subs(a,c)
l_cbij = l_cbij.subs(J_A, 0)
l_cbij = l_cbij.subs(M_A, 0)
l_cbij_sph = l_abij_sph.subs(a,c)
l_cbij_sph = l_cbij_sph.subs(J_A, 0)
l_cbij_sph = l_cbij_sph.subs(M_A, 0)
r_b = r_a.subs(Dagger(a), Dagger(b))
r_b_sph = r_a_sph.subs(Dagger(a), Dagger(b))
t_acij = t_abij.subs(Dagger(b), Dagger(c))
t_acij_sph = t_abij_sph.subs(Dagger(b), Dagger(c))

expr_msc = combine_ASigmas(S(1)/2*l_cbij*r_b*t_acij*ASigma('m_c','m_b','m_i','m_j') * sf_reduction.subs([(J_A, 0), (M_A, 0)]))
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_cbij_sph, r_b_sph, t_acij_sph])
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_i', 'm_j'])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, [M_Ap1, m_c])
expr_sph = refine_phases(expr_sph)
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))

print
print ("\subsubsection{simplified diagram $\\frac12l^{ij}_{cb}r^{cb}_it^a_j$}")
print
l_cbij = l_abij.subs(a,c)
l_cbij = l_cbij.subs(J_A, 0)
l_cbij = l_cbij.subs(M_A, 0)
l_cbij_sph = l_abij_sph.subs(a,c)
l_cbij_sph = l_cbij_sph.subs(J_A, 0)
l_cbij_sph = l_cbij_sph.subs(M_A, 0)
r_cbi = r_abi.subs(Dagger(a), Dagger(c))
r_cbi_sph = r_abi_sph.subs(Dagger(a), Dagger(c))
t_aj = t_ai.subs(i, j)
t_aj_sph = t_ai_sph.subs(i, j)

expr_msc = combine_ASigmas(
        S(1)/2*l_cbij*r_cbi*t_aj*ASigma('m_j','m_i','m_c', 'm_b')
        * sf_reduction.subs([(J_A, 0), (M_A, 0)]))
_report(Eq(SF, expr_msc))
expr_sph = rewrite_coupling(expr_msc, [l_cbij_sph, r_cbi_sph, t_aj_sph])
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_c', 'm_b'])
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph, all_deltas=True)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, [M_Ap1, m_i])
_report(Eq(SF, expr_sph))
_report(fcode(expr_sph))
