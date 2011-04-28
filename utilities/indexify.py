from sympy.physics.secondquant import AntiSymmetricTensor
from sympy.tensor import IndexedBase, Idx
from sympy.core.symbol import Symbol, Dummy
from sympy.core import Tuple, symbols


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

        elif at.symbol in symbols('r'):
            shp = (total_orbs,)*len(at.upper) + (below_orbs,)*len(at.lower)
        elif at.symbol in symbols('l'):
            # shp = (below_orbs,)*len(at.upper) + (total_orbs,)*len(at.lower)

            # change orbit_order for speed
            shp = (total_orbs,)*len(at.lower) + (below_orbs,)*len(at.upper)
            inds = Tuple(*(at.lower.args + at.upper.args))
            inds = inds.subs(subsdict)
        else:
            raise ValueError("Couldn't determine shape of %s" % at)

        base = IndexedBase("%s_%i" %(at.symbol, (len(inds)+1)/2), shape=shp)

        new_at = base[inds.args]
        expr = expr.subs(at, new_at)

        print "substituted %s for %s" %(new_at, at)

    return expr
