from sympy.physics.racahalgebra import (
        refine_phases, refine_tjs2sjs, convert_cgc2tjs, is_equivalent,
        SixJSymbol, ASigma, combine_ASigmas, evaluate_sums, apply_deltas,
        apply_orthogonality, ClebschGordanCoefficient
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
        pprint, Eq, pprint_use_unicode, latex, preview, fcode
        )
import sys
pprint_use_unicode(False)

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

def make_spherical_sp_states(str1):
    states = []
    for i in symbols(str1):
        states.append(SphFermKet(i))
    return states

def redmat(M):
    assert isinstance(M, ThreeTensorMatrixElement)
    return ReducedMatrixElement(M.left, M.operator, M.right)

def _report(expr):

    pprint(str(expr))

    # print fcode(expr)

    # print(expr)
    # pprint(expr)
    # print latex(expr)
    # preview(expr)

    print


SF = Symbol('SF')


Ket = FermKet
Bra = FermBra

i, j, k, l = make_spherical_sp_states('i j k l')
a, b, c, d = make_spherical_sp_states('a b c d')

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

print
print "Defined tensor operators:"
print LA
print RA
print LAm1
print RAm1
print V
print T
print cre_p
print ann_p

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

r_i = MatrixElement("",RAm1, i)
r_aij = MatrixElement(Dagger(a), RAm1, (i, j))
r_i_sph = ReducedMatrixElement("",RAm1, i)
r_aij_sph = ReducedMatrixElement(Dagger(a), RAm1, SphFermKet('ij',i, j))

l_ai = MatrixElement(Dagger(i),LA,a)
l_ai_sph = ReducedMatrixElement(Dagger(i),LA,a)
l_abij = MatrixElement((Dagger(i),Dagger(j)),LA,(a,b))
l_abij_sph = ReducedMatrixElement(SphFermBra('ij',Dagger(i),Dagger(j), reverse=0), LA, SphFermKet('ab', a, b,reverse=0))

t_ai = MatrixElement(Dagger(a),T ,i)
t_ai_sph = ReducedMatrixElement(Dagger(a),T ,i)
t_abij = MatrixElement((Dagger(a), Dagger(b)), T,(i, j))
t_abij_sph = ReducedMatrixElement(SphFermBra('ab',Dagger(a), Dagger(b), reverse=0), T, SphFermKet('ij', i, j,reverse=0))


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
sf_reduction = (-1)**(j_i - m_i)*ClebschGordanCoefficient('J_A', 'M_A', j_i, -m_i, 'J_Am1', 'M_Am1')*ASigma('M_A', 'm_i')


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
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_phases(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)

"""

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

print
_report("recoupling diagram 7")
l_abjk = l_abij.subs(Dagger(j),Dagger(k))
l_abjk = l_abjk.subs(Dagger(i),Dagger(j))
l_abjk_sph = l_abij_sph.subs(Dagger(j),Dagger(k))
l_abjk_sph = l_abjk_sph.subs(Dagger(i),Dagger(j))
r_j = r_i.subs(i,j)
r_j_sph = r_i_sph.subs(i,j)
t_abik = t_abij.subs(j,k)
t_abik_sph = t_abij_sph.subs(j,k)
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
expr_sph = convert_cgc2tjs(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = evaluate_sums(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = apply_deltas(expr_sph)
_report(Eq(SF, expr_sph))
expr_sph = refine_tjs2sjs(expr_sph, verbose=1)
_report(Eq(SF, expr_sph))
print fcode(expr_sph)

