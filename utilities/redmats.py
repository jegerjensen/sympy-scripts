from sympy.physics.braket import (
        ReducedMatrixElement, ClebschGordanCoefficient, ASigma
        )

class ReducedMatrixElement_right(ReducedMatrixElement):
    _definition = "right"

    def _eval_reduction_factor(self):
        left, op, right = self.args
        c_ket, t_ket = right.as_coeff_tensor()
        j1, m1 = t_ket.get_rank_proj()
        c_op, t_op = op.as_coeff_tensor()
        j2, m2 = t_op.get_rank_proj()
        J, M = left._j, left._m
        return c_op*c_ket*ClebschGordanCoefficient(j2, m2, j1, m1, J, M)
    def _eval_inverse_reduction_factor(self):
        left,op,right = self.args
        c_ket, t_ket = right.as_coeff_tensor()
        j1, m1 = t_ket.get_rank_proj()
        c_op, t_op = op.as_coeff_tensor()
        j2, m2 = t_op.get_rank_proj()
        J, M = left._j, left._m
        return ASigma(m1, m2)*ClebschGordanCoefficient(j2, m2, j1, m1, J, M)/c_op/c_ket

class ReducedMatrixElement_left(ReducedMatrixElement):
    _definition = "left"

    def _eval_reduction_factor(self):
        left, op, right = self.args
        c_bra, t_bra = left.as_coeff_tensor()
        j1, m1 = t_bra.get_rank_proj()
        c_op, t_op = op.as_coeff_tensor()
        j2, m2 = t_op.get_rank_proj()
        J, M = right._j, right._m
        return c_op*c_bra*ClebschGordanCoefficient(j1, m1, j2, m2, J, -M)*(-1)**(J-M)

    def _eval_inverse_reduction_factor(self):
        left,op,right = self.args
        left, op, right = self.args
        c_bra, t_bra = left.as_coeff_tensor()
        j1, m1 = t_bra.get_rank_proj()
        c_op, t_op = op.as_coeff_tensor()
        j2, m2 = t_op.get_rank_proj()
        J, M = right._j, right._m
        return ASigma(m1, m2)*ClebschGordanCoefficient(j1, m1, j2, m2, J, -M)/c_op/c_bra/(-1)**(J-M)
import sympy
sympy.physics.braket.ReducedMatrixElement_left = ReducedMatrixElement_left
sympy.physics.braket.ReducedMatrixElement_right = ReducedMatrixElement_right

