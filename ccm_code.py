#!/usr/bin/env python
"""
Calculates the Coupled-Cluster energy- and amplitude equations
See 'An Introduction to Coupled Cluster Theory' by
T. Daniel Crawford and Henry F. Schaefer III.
http://www.ccc.uga.edu/lec_top/cc/html/review.html
"""

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator)
from sympy import (
    symbols, expand, pprint, Number, latex, Symbol, Tuple, Eq
)
from sympy import IndexedBase, Idx, fcode, S
from sympy.utilities.codegen import codegen
from sympy.core.symbol import Dummy

from utilities.ccm import get_CC_operators

pretty_dummies_dict={
        'above':'cdefgh',
        'below':'klmno',
        'general':'pqrstu'
        }
project_descr = 'CCSD_2bodyV'


def main():
    print
    print "Calculates the Coupled-Cluster energy- and amplitude equations"
    print "See 'An Introduction to Coupled Cluster Theory' by"
    print "T. Daniel Crawford and Henry F. Schaefer III"
    print "http://www.ccc.uga.edu/lec_top/cc/html/review.html"
    print

    # setup hamiltonian
    p,q,r,s = symbols('pqrs',dummy=True)
    f = AntiSymmetricTensor('f',(p,),(q,))
    pr = NO((Fd(p)*F(q)))
    v = AntiSymmetricTensor('v',(p,q),(r,s))
    pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))

    H=f*pr

    # Uncomment the next line to use a 2-body hamiltonian:
    H=f*pr + Number(1,4)*v*pqsr

    print "Using the hamiltonian:", latex(H)

    print "Calculating nested commutators"
    C = Commutator

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm1..."
    comm1 = wicks(C(H,T))
    comm1 = evaluate_deltas(comm1)
    comm1 = substitute_dummies(comm1)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm2..."
    comm2 = wicks(C(comm1,T))
    comm2 = evaluate_deltas(comm2)
    comm2 = substitute_dummies(comm2)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm3..."
    comm3 = wicks(C(comm2,T))
    comm3 = evaluate_deltas(comm3)
    comm3 = substitute_dummies(comm3)

    T1,T2 = get_CC_operators()
    T = T1+ T2
    print "comm4..."
    comm4 = wicks(C(comm3,T))
    comm4 = evaluate_deltas(comm4)
    comm4 = substitute_dummies(comm4)

    print "construct Hausdoff expansion..."
    eq = H + comm1+comm2/2  +comm3/6+comm4/24
    eq = eq.expand()
    eq = evaluate_deltas(eq)
    eq = substitute_dummies(eq, new_indices=True,
            pretty_indices=pretty_dummies_dict)
    eq = eq.evalf()
    print "*********************"
    print

    print "extracting CC equations from full Hbar"
    i,j,k,l = symbols('ijkl',below_fermi=True)
    a,b,c,d = symbols('abcd',above_fermi=True)
    print
    print "CC Energy:"
    eqEnergy = wicks(eq, simplify_dummies=True, keep_only_fully_contracted=True)
    print latex(eqEnergy)

    eqEnergy = indexify(eqEnergy)
    print eqEnergy
    codegen(('CC_energy', eqEnergy), 'f95', 'code/sympy_ccm_Energy', to_files=True, project=project_descr)


    print
    print "CC T1:"
    eqT1 = wicks(NO(Fd(i)*F(a))*eq, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
    eqT1 = substitute_dummies(eqT1)
    print latex(eqT1)

    lhs = AntiSymmetricTensor('lhs', (a,),  (i,))
    eqT1_idx = indexify(Eq(lhs, eqT1))
    print eqT1_idx
    codegen(('t1_diagrams', eqT1_idx), 'f95', 'code/sympy_ccm_T1', to_files=True, project=project_descr)


    print
    print "CC T2:"
    eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*eq,simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
    print latex(eqT2)

    lhs = AntiSymmetricTensor('lhs', (a, b),  (i, j))
    eqT2_idx = indexify(Eq(lhs, eqT2))
    print eqT2_idx
    codegen(('t2_diagrams', eqT2_idx), 'f95', 'code/sympy_ccm_T2', to_files=True, project=project_descr)


def indexify(expr):
    """converts an expression with AntiSymmetricTensor to Indexed and Idx

    This is useful for generation of code for Coupled-Cluster
    """

    all_AT = expr.atoms(AntiSymmetricTensor)

    total_orbs, below_orbs = symbols('total_orbs below_orbs', integer=True)

    subsdict = {}       # for index substitutions
    fixed_inds = set()  # indices for left hand side
    for at in all_AT:

        for index in at.upper.args + at.lower.args:
            if not isinstance(index, Dummy):
                fixed_inds.add(index)
                label = Symbol(index.name, integer=True)
            else:
                label = Symbol("%s_%i" %(index.name, index.dummy_index),
                        integer=True)

            if index not in subsdict:
                # determine index range
                if index.assumptions0.get('above_fermi'):
                    start = below_orbs
                else:
                    start = 0
                if index.assumptions0.get('below_fermi'):
                    stop = below_orbs -1
                else:
                    stop = total_orbs -1

                idx = Idx(label, (start, stop))
                subsdict[index] = idx

    # substitute one AntiSymmetricTensor at a time
    for at in all_AT:
        inds = Tuple(*(at.upper.args + at.lower.args))
        inds = inds.subs(subsdict)

        # define arrays with shape because it cannot always be determined from
        # indices
        if at.symbol in symbols('t lhs'):
            shp = (total_orbs,)*(len(inds)%2) + (below_orbs,)*(len(inds)%2)
        elif at.symbol in symbols('fvwh'):
            shp = (total_orbs,)*len(inds)
        else:
            raise ValueError("Couldn't determine shape of %s" % at)

        base = IndexedBase("%s_%i" %(at.symbol, len(inds)), shape=shp)

        new_at = base[inds.args]
        expr = expr.subs(at, new_at)

        print "substituted %s for %s" %(new_at, at)

    return expr


if __name__ == "__main__":
    main()
