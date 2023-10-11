import pymc

from vtl_common.common_GUI.tk_forms_assist import FormNode, AlternatingNode


from ..form_prototypes.distribution_parameter import DistBuilder


class FixedCaller(object):
    def __init__(self, data):
        self.data = data

    def call(self):
        raise NotImplementedError

    def __call__(self):
        return self.call()


class TauCaller(FixedCaller):
    def call(self):
        tau = self.data("tau")
        return pymc.Deterministic("u_z_", -1 / tau)


class UzCaller(FixedCaller):
    def call(self):
        u_z = self.data("u_z")
        return pymc.Deterministic("u_z_",u_z)


def create_z_alter():
    class TauField(DistBuilder):
        DISPLAY_NAME = "tau"
        def get_data(self):
            data = super().get_data()
            return TauCaller(data)

    class UzField(DistBuilder):
        DISPLAY_NAME = "u_z"
        def get_data(self):
            data = super().get_data()
            return UzCaller(data)

    class ZAlter(AlternatingNode):
        DISPLAY_NAME = "z_correction"
        SEL__u_z = UzField
        SEL__tau = TauField

    return ZAlter