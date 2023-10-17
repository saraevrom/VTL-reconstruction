import pymc

from vtl_common.common_GUI.tk_forms_assist import FormNode, AlternatingNode, LabelNode
from fixed_rotator.slowmath import Vector2,Quaternion

from ..form_prototypes.distribution_parameter import DistBuilder


class FixedCaller(object):
    def __init__(self, data):
        self.data = data

    def call(self,R0_vec, U0_vec, R_quat):
        raise NotImplementedError

    def __call__(self, R0_vec, U0_vec, R_quat):
        return self.call(R0_vec, U0_vec, R_quat)


class TauCaller(FixedCaller):
    def call(self,R0_vec, U0_vec, R_quat):
        tau = self.data("tau")
        return pymc.Deterministic("nu_", 1 / tau)


class NuCaller(FixedCaller):
    def call(self, R0_vec, U0_vec, R_quat):
        nu = self.data("nu")
        return pymc.Deterministic("nu_",nu)


class DeterministicCaller(FixedCaller):
    def call(self, R0_vec, U0_vec, R_quat):
        R_mat = R_quat.to_matrix()
        rho_col = R_mat[:,2]
        rho = Vector2(rho_col[0],rho_col[1])/rho_col[2]
        nu = pymc.Deterministic("nu",U0_vec.dot(rho)/(1+R0_vec.dot(rho)))
        return pymc.Deterministic("nu_",nu)


def create_z_alter():
    class TauField(DistBuilder):
        DISPLAY_NAME = "tau"
        def get_data(self):
            data = super().get_data()
            return TauCaller(data)

    class NuField(DistBuilder):
        DISPLAY_NAME = "nu"
        def get_data(self):
            data = super().get_data()
            return NuCaller(data)

    class DeterministicField(LabelNode):
        DISPLAY_NAME = "Horizontal Fix"

        def get_data(self):
            return DeterministicCaller(None)

    class ZAlter(AlternatingNode):
        DISPLAY_NAME = "z_correction"
        SEL__nu = NuField
        SEL__tau = TauField
        SEL__horizontal = DeterministicField

    return ZAlter