from sympy.core.symbol import Symbol, Dummy

class AmbiguousSymbolException(Exception):
    pass

def undummy(expr):
    """replaces all Dummy by corresponding Symbol

    We take care that the corresponding symbol is not already
    present in the expression.  In that case we raise an exception.
    """

    symbols = expr.atoms(Symbol)
    dummies = expr.atoms(Dummy)
    subslist = []
    for d in dummies:
        assum = d.assumptions0
        s = Symbol(d.name, **assum)
        if expr.has(s):
            raise AmbiguousSymbolException(
                    "Symbol `%s' already present in expression"%s
                    )
        else:
            subslist.append((d, s))
    return expr.subs(subslist)




