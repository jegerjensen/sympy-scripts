from sympy.physics.racahalgebra import (
        refine_phases, refine_tjs2sjs, convert_cgc2tjs, is_equivalent,
        SixJSymbol, ASigma, combine_ASigmas, evaluate_sums, apply_deltas,
        apply_orthogonality, ClebschGordanCoefficient, extract_symbol2dummy_dict
        )
from sympy.physics.braket import (
        MatrixElement, ReducedMatrixElement, apply_wigner_eckardt,
        rewrite_as_direct_product, DirectMatrixElement,
        SphericalTensorOperator, Dagger, SphFermKet, SphFermBra,
        FermKet, FermBra, ThreeTensorMatrixElement,
        inject_every_symbol_globally, rewrite_coupling
        )
from sympy import (
        Symbol, symbols, global_assumptions, Assume, ask, Q, Mul, S, sqrt,
        pprint, Eq, pprint_use_unicode, latex, preview, fcode, Function
        )
import sys

pprint_use_unicode(False)
def _report(expr):

    pprint(str(expr))

    # print fcode(expr)

    # print(expr)
    # pprint(expr)
    # print latex(expr)
    # preview(expr)

    print

def check_agains_gaute(expr_sjs):
    inject_every_symbol_globally(expr_sjs, quiet=0, force=True)
    Jab, j_eom, jt = symbols('Jab j_eom jt')
    gautes = SixJSymbol(j_i,j_eom,Jab, j_b,j_j,jt)*(2*jt+1)*(-1)**(Jab+j_j+j_eom+jt-jt+j_b+j_eom+1)*ASigma(jt)
    gautes = gautes.subs({jt:J_ij, j_eom:J_Am1, Jab:J_bj})
    factor, junk = expr_sjs.as_coeff_terms(ReducedMatrixElement)
    print "checking  equivalence of"
    print gautes
    print "and\n", factor
    if is_equivalent(gautes, factor, True):
        print "Recoupling expressions are matching!"
        print gautes
        print factor
    else:
        print "Not matching"

def make_tensor_operators(str1, assumption_str=None):
    rank_proj = symbols(str1)
    if assumption_str:
        for s in rank_proj:
            global_assumptions.add(Assume(s,assumption_str))

    tensors = []
    for i in range(0,len(rank_proj),2):
        tensors.append(SphericalTensorOperator('t',rank_proj[i], rank_proj[i+1]))

    return tensors

i, j, k, l = map(SphFermKet, 'ijkl')
a, b, c, d = map(SphFermKet, 'abcd')

LA, RA = make_tensor_operators('J_A M_A J_A M_A',Q.integer)
LAm1,RAm1 = make_tensor_operators('J_Am1 M_Am1 J_Am1 M_Am1','half_integer')
LA = Dagger(LA.subs(Symbol('t'), Symbol('L')))
# LA = SphericalTensorOperator('L', S.Zero, S.Zero)
RA = RA.subs(Symbol('t'), Symbol('R'))
LAm1 = LAm1.subs(Symbol('t'), Symbol('L'))
RAm1 = RAm1.subs(Symbol('t'), Symbol('R'))
V = SphericalTensorOperator('V', S.Zero, S.Zero)
T = SphericalTensorOperator('T', S.Zero, S.Zero)

print
print "Defined tensor operators:"
print LA
print RA
print LAm1
print RAm1
print V
print T


print
print "*************** <i| L | a> *****************"
print

l_ai = MatrixElement(Dagger(i), LA, a)
l_ai_sph = ReducedMatrixElement(Dagger(i), LA, a)
print Eq(l_ai, l_ai_sph.get_direct_product_ito_self(tjs=0))

print
print fcode(l_ai_sph.get_direct_product_ito_self(tjs=1))

print
print "*************** <ij| L | ab> *****************"
print

l_abij = MatrixElement((Dagger(i), Dagger(j)), LA, (a, b))
l_abij_sph = ReducedMatrixElement(SphFermBra('ij', Dagger(i), Dagger(j), reverse=0), LA, SphFermKet('ab', a, b, reverse=0))
print Eq(l_abij, l_abij_sph.get_direct_product_ito_self(tjs=0))

print
print fcode(l_abij_sph.get_direct_product_ito_self(tjs=1))

print
print "*************** <a| T | i> *****************"
print

t_ai = MatrixElement(Dagger(a), T, i)
t_ai_sph = ReducedMatrixElement(Dagger(a), T, i)
print Eq(t_ai, t_ai_sph.get_direct_product_ito_self(tjs=0))

print
print fcode(t_ai_sph.get_direct_product_ito_self(tjs=1))

print
print "*************** <ab| T | ij> *****************"
print

t_abij = MatrixElement((Dagger(a), Dagger(b)), T, (i, j))
t_abij_sph = ReducedMatrixElement(SphFermBra('ab', Dagger(a), Dagger(b), reverse=0), T, SphFermKet('ij', i, j, reverse=0))
print Eq(t_abij, t_abij_sph.get_direct_product_ito_self(tjs=0))

print
print fcode(t_abij_sph.get_direct_product_ito_self(tjs=1))



print
print "*************** <| R | i> *****************"
print

r_i = MatrixElement(0, RAm1, i)
r_i_sph = ReducedMatrixElement(0, RAm1, i)
print Eq(r_i, r_i_sph.get_direct_product_ito_self(tjs=0))
print Eq(r_i, r_i_sph.get_direct_product_ito_self(tjs=1))

print
print fcode(r_i_sph.get_direct_product_ito_self(tjs=1))
print
print "*************** <a| R | ij> *****************"
print

# r_aij = MatrixElement((Dagger(i), Dagger(j)), RAm1, b)
# r_aij_sph = ReducedMatrixElement(SphFermBra('ij', Dagger(i), Dagger(j), reverse=1), RAm1, b)
r_aij = MatrixElement(Dagger(a), RAm1, (i, j))
r_aij_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ij', i, j, reverse=0))
print Eq(r_aij, r_aij_sph.get_direct_product_ito_self(tjs=0))

print
print fcode(r_aij_sph.get_direct_product_ito_self(tjs=1))

# As Gaute told me on chat:
#  <b|R|ij> = (ji,mi,jj,mj|Jij Mij) (jb,mb,J,M|Jij Mij) <b||R||ij>
J_ij, J_Am1, j_b, j_i, j_j = symbols('J_ij J_Am1 j_b j_i j_j', nonnegative=True)
M_ij, M_Am1, m_b, m_i, m_j = symbols('M_ij M_Am1 m_b m_i m_j')

J_Am1 = RAm1.rank


print
print
print
gautes_raij = (ASigma(m_b, M_Am1, m_i, m_j)
        *MatrixElement(Dagger(b), RAm1, (i, j))
        *ClebschGordanCoefficient(j_b, m_b, J_Am1, M_Am1, J_ij, M_ij)
        *ClebschGordanCoefficient(j_i, m_i, j_j, m_j, J_ij, M_ij)
        )
print "As Gaute told me on chat:"
print ' "<b|R|ij> = (ji,mi,jj,mj|Jij Mij) (jb,mb,J,M|Jij Mij) <b||R||ij>" '
print " so Gautes reduced matrix element can be written"
print gautes_raij
print "The other way around, we can express the direct matrix as"
gautes = Function('Gautes')
direct_as_gautes = ( ASigma(J_ij, M_ij)
        *ClebschGordanCoefficient(j_b, m_b, J_Am1, M_Am1, J_ij, M_ij)
        *ClebschGordanCoefficient(j_i, m_i, j_j, m_j, J_ij, M_ij)
        *gautes(Dagger(b), SphFermKet('ij',i,j))
        )
print direct_as_gautes
print """
Here, we have inserted the R operator directly, but in order to make Gautes
reduced matrix element independent of projection, it is necessary to treat the
R operator as a dual SphericalTensorOperator, transforming like (-1)**(j-m)T(j,-m).
In the wikipedian reduced matrix element, we expect R to transform as a
tensor T(j,m)i.  Since we prefer to insert the R operator directly in the
wikipedian element, the R operator in gautes must be substituted with
(-1)**(j+m)R(j,-m).  The expression is then """
direct_as_gautes = ( ASigma(J_ij, M_ij)*(-1)**(J_Am1 + M_Am1)
        *ClebschGordanCoefficient(j_b, m_b, J_Am1, -M_Am1, J_ij, M_ij)
        *ClebschGordanCoefficient(j_i, m_i, j_j, m_j, J_ij, M_ij)
        *gautes(Dagger(b), SphFermKet('ij',i,j))
        )
print direct_as_gautes
print "We should use dummies for summation indices"
direct_as_gautes = direct_as_gautes.subs([(J_ij, J_ij.as_dummy()), (M_ij, M_ij.as_dummy())])
print direct_as_gautes
print "The wikipedian reduced matrix element can also be written in terms of the"
print "direct product:"
mine_as_direct = r_aij_sph.as_direct_product()
print mine_as_direct
subsdict = extract_symbol2dummy_dict(mine_as_direct)
direct_as_gautes = direct_as_gautes.subs(subsdict)
direct_w_dummies = r_aij_sph.get_related_direct_matrix().subs(subsdict)
relation_mine_ito_gautes = mine_as_direct.subs(direct_w_dummies, direct_as_gautes)
print "Inserting the expression for direct in terms of Gautes gives the relation"
print relation_mine_ito_gautes
print "Which we can simplify to"
relation_mine_ito_gautes = combine_ASigmas(relation_mine_ito_gautes)
relation_mine_ito_gautes = apply_orthogonality(relation_mine_ito_gautes, [ subsdict[m_i], subsdict[m_j]])
print relation_mine_ito_gautes
print relation_mine_ito_gautes
relation_mine_ito_gautes = evaluate_sums(relation_mine_ito_gautes)
print relation_mine_ito_gautes
relation_mine_ito_gautes = convert_cgc2tjs(relation_mine_ito_gautes)
relation_mine_ito_gautes = refine_phases(relation_mine_ito_gautes)
print relation_mine_ito_gautes
subsdict = extract_symbol2dummy_dict(relation_mine_ito_gautes)
# print "To use the other orthogonality, we need a sum over 'm_b'"
# relation_mine_ito_gautes = relation_mine_ito_gautes*ASigma(m_b)/(2*j_b + 1)
# print relation_mine_ito_gautes
# relation_mine_ito_gautes = apply_orthogonality(relation_mine_ito_gautes, [ subsdict[M_ij], m_b])
relation_mine_ito_gautes = apply_orthogonality(relation_mine_ito_gautes, [ subsdict[M_ij], subsdict[M_Am1]])
print relation_mine_ito_gautes
# print "And finally, we can remove the summation over M_Am1 multipliyng with (2J_Am1 + 1)"
# relation_mine_ito_gautes = evaluate_sums(relation_mine_ito_gautes, independent_of={subsdict[M_ij]:2*subsdict[J_ij]+1, subsdict[M_Am1]: 2*J_Am1+1 })
# print relation_mine_ito_gautes
relation_mine_ito_gautes = refine_phases(relation_mine_ito_gautes)
print relation_mine_ito_gautes
print
print fcode(relation_mine_ito_gautes)





SF = Symbol('SF')

"""




LA, RA = make_tensor_operators('J_A M_A J_A M_A',Q.integer)
LAm1,RAm1 = make_tensor_operators('J_Am1 M_Am1 J_Am1 M_Am1','half_integer')
# LA = LA.subs(Symbol('t'), Symbol('L'))
LA = Dagger(LA.subs(Symbol('t'), Symbol('L')))
# LA = SphericalTensorOperator('L', S.Zero, S.Zero)
RA = RA.subs(Symbol('t'), Symbol('R'))
LAm1 = Dagger(LAm1.subs(Symbol('t'), Symbol('L')))
RAm1 = RAm1.subs(Symbol('t'), Symbol('R'))
V = SphericalTensorOperator('V', S.Zero, S.Zero)
T = SphericalTensorOperator('T', S.Zero, S.Zero)
cre_p = SphericalTensorOperator('a','j','m')
ann_p = Dagger(SphericalTensorOperator('a','j','m'))
"""
print
print "Defined tensor operators:"
print LA
print RA
print LAm1
print RAm1
print V
print T

print
print "==== Recoupling L_{A} a' R_{A-1} ===="
print

# """
j_a = Symbol('j_a', nonnegative=True)
m_a = Symbol('m_a')
sf_reduction = (-1)**(j_a - m_a)*ClebschGordanCoefficient('J_A', 'M_A', j_a, -m_a, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_a')
# sf_reduction = (-1)**(j_a - m_a)*ClebschGordanCoefficient(0, 0, j_a, -m_a, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_a')

print "reduction factor"
print sf_reduction
print fcode(refine_phases(convert_cgc2tjs(sf_reduction)))

l0 = Symbol('l_0')
r0 = Symbol('r_0')

# r_i = MatrixElement("",RAm1, i)
# r_aij = MatrixElement(Dagger(a), RAm1, (i, j))
# r_i_sph = ReducedMatrixElement("",RAm1, i)
# r_aij_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ij',i, j))

# l_ai = MatrixElement(Dagger(i),LA,a)
# l_ai_sph = ReducedMatrixElement(Dagger(i),LA,a)
# l_abij = MatrixElement((Dagger(i),Dagger(j)),LA,(a,b))
# l_abij_sph = ReducedMatrixElement(SphFermBra('ij',Dagger(i),Dagger(j), reverse=0), LA, SphFermKet('ab', a, b,reverse=0))

# t_ai = MatrixElement(Dagger(a),T ,i)
# t_ai_sph = ReducedMatrixElement(Dagger(a),T ,i)
# t_abij = MatrixElement((Dagger(a), Dagger(b)), T,(i, j))
# t_abij_sph = ReducedMatrixElement(SphFermBra('ab',Dagger(a), Dagger(b), reverse=0), T, SphFermKet('ij', i, j,reverse=0))

_report("*** coupled elements: ***")
_report(Eq(r_i_sph,rewrite_coupling(r_i_sph, r_i)))
_report(Eq(r_aij_sph,rewrite_coupling(r_aij_sph, r_aij)))
_report(Eq(l_ai_sph,rewrite_coupling(l_ai_sph, l_ai)))
_report(Eq(l_abij_sph,rewrite_coupling(l_abij_sph, l_abij)))
_report(Eq(t_ai_sph,rewrite_coupling(t_ai_sph, t_ai)))
_report(Eq(t_abij_sph,rewrite_coupling(t_abij_sph, t_abij)))

"""
coupled_subs = {
        l_ai: rewrite_coupling(l_ai, l_ai_sph),
        r_i: rewrite_coupling(r_i, r_i_sph)
        }
print
_report("recoupling diagram 1")
expr_msc = combine_ASigmas(l_ai*r_i*ASigma('m_i') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = apply_deltas(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)


print
_report("recoupling diagram 2")
r_bij = r_aij.subs(Dagger(a), Dagger(b))
r_bij_sph = r_aij_sph.subs(Dagger(a), Dagger(b))

coupled_subs = {
        l_abij: rewrite_coupling(l_abij, l_abij_sph),
        r_bij: rewrite_coupling(r_bij, r_bij_sph)
        }

expr_msc = combine_ASigmas(S(1)/2*l_abij*r_bij*ASigma('m_b','m_i','m_j') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_b', 'm_i', 'm_j'])
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases((expr_sph))
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1, let_pass=0)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)


print
print "==== Recoupling L_{A} i' R_{A-1} ===="
print


j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
sf_reduction = 1/((-1)**(j_i - m_i)*ClebschGordanCoefficient(0, 0, j_i, -m_i, 'J_Am1', 'M_Am1').get_as_ThreeJSymbol())


coupled_subs = {
        r_i: rewrite_coupling(r_i, r_i_sph)
        }
print
_report("recoupling diagram 3")
expr_msc = combine_ASigmas(l0*r_i*sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)



print
_report("recoupling diagram 4")
l_aj = l_ai.subs(Dagger(i),Dagger(j))
l_aj_sph = l_ai_sph.subs(Dagger(i),Dagger(j))
coupled_subs = {
        l_aj: rewrite_coupling(l_aj, l_aj_sph),
        r_aij: rewrite_coupling(r_aij, r_aij_sph)
        }
j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient('J_A', 'M_A', j_i, -m_i, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_i')
print
expr_msc = combine_ASigmas(l_aj*r_aij*ASigma('m_a','m_j') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)

print
_report("recoupling diagram 5")
l_aj = l_ai.subs(Dagger(i),Dagger(j))
l_aj_sph = l_ai_sph.subs(Dagger(i),Dagger(j))
r_j = r_i.subs(i,j)
r_j_sph = r_i_sph.subs(i,j)
coupled_subs = {
        t_ai: rewrite_coupling(t_ai, t_ai_sph),
        l_aj: rewrite_coupling(l_aj, l_aj_sph),
        r_j: rewrite_coupling(r_j, r_j_sph)
        }
j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient('J_A', 'M_A', j_i, -m_i, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_i')
print
expr_msc = combine_ASigmas(l_aj*r_j*t_ai*ASigma('m_a','m_j') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = apply_deltas(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph, forbidden=[S('M_A'), S('M_Am1'), S('m_a')])
_report(Eq(SF, expr_sph))
print fcode(expr_sph)

"""

print
_report("recoupling diagram 6")
l_abjk = l_abij.subs(Dagger(j),Dagger(k))
l_abjk = l_abjk.subs(Dagger(i),Dagger(j))
l_abjk_sph = l_abij_sph.subs(SphFermBra('ij',Dagger(i),Dagger(j)), SphFermBra('jk', Dagger(j), Dagger(k)))
r_bjk = r_aij.subs(j,k)
r_bjk = r_bjk.subs(i,j)
r_bjk = r_bjk.subs(Dagger(a), Dagger(b))
r_bjk_sph = r_aij_sph.subs({Dagger(a): Dagger(b), SphFermKet('ij',i,j): SphFermKet('jk',j,k)})
coupled_subs = {
        t_ai: rewrite_coupling(t_ai, t_ai_sph),
        l_abjk: rewrite_coupling(l_abjk, l_abjk_sph),
        r_bjk: rewrite_coupling(r_bjk, r_bjk_sph)
        }
j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient('J_A', 'M_A', j_i, -m_i, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_i')
print
expr_msc = combine_ASigmas(l_abjk*r_bjk*t_ai*ASigma('m_a','m_j','m_k','m_b') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_j', 'm_k'])
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = apply_deltas(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)




print
_report("recoupling diagram 7")
l_abjk = l_abij.subs(Dagger(j),Dagger(k))
l_abjk = l_abjk.subs(Dagger(i),Dagger(j))
l_abjk_sph = l_abij_sph.subs(SphFermBra('ij',Dagger(i),Dagger(j)), SphFermBra('jk', Dagger(j), Dagger(k)))
r_j = r_i.subs(i,j)
r_j_sph = r_i_sph.subs(i,j)
t_abik = t_abij.subs(j,k)
t_abik_sph = t_abij_sph.subs(SphFermKet('ij',i,j), SphFermKet('ik',i,k))
coupled_subs = {
        t_abik: rewrite_coupling(t_abik, t_abik_sph),
        l_abjk: rewrite_coupling(l_abjk, l_abjk_sph),
        r_j: rewrite_coupling(r_j, r_j_sph)
        }
j_i = Symbol('j_i', nonnegative=True)
m_i = Symbol('m_i')
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient('J_A', 'M_A', j_i, -m_i, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_i')
print
expr_msc = combine_ASigmas(l_abjk*r_j*t_abik*ASigma('m_a','m_j','m_k','m_b') * sf_reduction)
_report(Eq(SF, expr_msc))
expr_sph = expr_msc.subs(coupled_subs)
_report(Eq(SF, expr_sph))
expr_sph = apply_orthogonality(expr_sph, ['m_b', 'm_a'])
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = apply_deltas(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)

