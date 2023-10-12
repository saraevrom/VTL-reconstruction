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
        return pymc.Deterministic("nu_", 1 / tau)


class NuCaller(FixedCaller):
    def call(self):
        nu = self.data("nu")
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

    class ZAlter(AlternatingNode):
        DISPLAY_NAME = "z_correction"
        SEL__nu = NuField
        SEL__tau = TauField

    return ZAlter