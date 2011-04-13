"""Module that defines common objects to spherical ccm"""

from sympy import Symbol, Assume, Q, S
from sympy.physics.braket import (
        braket_assumptions, DualSphericalTensorOperator, SphericalTensorOperator
        )

J_A = Symbol('J_A', nonnegative=True)
M_A = Symbol('M_A')
braket_assumptions.add(Assume(J_A, Q.integer))
braket_assumptions.add(Assume(M_A, Q.integer))

J_Ap1 = Symbol('J_Ap1', nonnegative=True)
M_Ap1 = Symbol('M_Ap1')
braket_assumptions.add(Assume(J_Ap1, 'half_integer'))
braket_assumptions.add(Assume(M_Ap1, 'half_integer'))

J_Am1 = Symbol('J_Am1', nonnegative=True)
M_Am1 = Symbol('M_Am1')
braket_assumptions.add(Assume(J_Am1, 'half_integer'))
braket_assumptions.add(Assume(M_Am1, 'half_integer'))

LA = DualSphericalTensorOperator('L', J_A, M_A)
LAp1 = DualSphericalTensorOperator('L', J_Ap1, M_Ap1)
LAm1 = DualSphericalTensorOperator('L', J_Am1, M_Am1)
RA = SphericalTensorOperator('R', J_A, M_A)
RAp1 = SphericalTensorOperator('R', J_Ap1, M_Ap1)
RAm1 = SphericalTensorOperator('R', J_Am1, M_Am1)
V = SphericalTensorOperator('V', S.Zero, S.Zero)
T = SphericalTensorOperator('T', S.Zero, S.Zero)
