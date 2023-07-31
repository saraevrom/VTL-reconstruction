import numpy as np
import pymc as pm

from ..model_base import FormPrototype
from ..form_prototypes import DistributionField
from vtl_common.common_GUI.tk_forms_assist import AlternatingNode

LC_COLOR = "#536878"

def E0_field():
    return DistributionField("uniform", lower=10.0, upper=60.0)

def tau_field():
    return DistributionField("halfnormal", sigma=1.0)

HLINE_STYLES = {
    "A" : dict(linestyles="solid", color=LC_COLOR, label="LC_A"),
    "B" : dict(linestyles="dashed", color=LC_COLOR, label="LC_B"),
    "C" : dict(linestyles="dashdot", color=LC_COLOR, label="LC_C"),
    "D" : dict(linestyles="dotted", color=LC_COLOR, label="LC_D"),
    "M" : dict(linestyles="solid", color="#555555", label="LC_all")
}

def _mod_style(s):
    s1 = s.copy()
    s1["linestyle"] = s1["linestyles"]
    del s1["linestyles"]
    return s1


PLOT_STYLES = {k: _mod_style(HLINE_STYLES[k]) for k in HLINE_STYLES.keys()}
del _mod_style

class LightCurve(FormPrototype):

    def get_lc(self, delta_k, k0, const_storage):
        raise NotImplementedError("This light curve is not implemented")

    def postprocess_plot(self, delta_k, k0, ax, trace, pmt, actual_x=None):
        if actual_x is None:
            actual_x = delta_k+k0
        self._postprocess(delta_k, k0, ax, trace, pmt, actual_x)

    def _postprocess(self, delta_k, k0, ax, trace, pmt, actual_x):
        pass

class ConstantLC(LightCurve):
    E0 = E0_field()
    def get_lc(self, delta_k, k0, const_storage):
        return self.E0("E0", const_storage)

    def _postprocess(self, delta_k, k0, ax, trace, pmt,actual_x):
        print(trace)
        e0 = trace.posterior["E0"].median()
        ax.hlines(e0, actual_x[0], actual_x[-1], **HLINE_STYLES[pmt[0]])

class LinearLC(LightCurve):
    TAU = DistributionField("normal", mu=0, sigma=1.0)
    E0 = E0_field()

    def get_lc(self, delta_k, k0, const_storage):
        tau = self.TAU("τ_LC", const_storage)
        e0 = self.E0("E0", const_storage)
        return e0*(1+delta_k/tau)

    def _postprocess(self, delta_k, k0, ax, trace, pmt, actual_x):
        tau = float(trace.posterior["τ_LC"].median())
        e0 = float(trace.posterior["E0"].median())
        # coeff = float(trace.posterior["K_LC"].median())
        # offset = float(trace.posterior["B_LC"].median())

        ax.plot(actual_x, e0*(1+delta_k/tau), **PLOT_STYLES[pmt[0]])

class GaussianLC(LightCurve):
    E0 = E0_field()
    mu_k0 = DistributionField("normal", mu=0, sigma=1.0)
    tau = tau_field()

    def get_lc(self, delta_k, k0, const_storage):
        e0 = self.E0("E0", const_storage)
        mu_k0 = self.mu_k0("mu_LC_k0", const_storage)
        mu = pm.Deterministic("mu_LC", k0+mu_k0)  # Maximum position from start of frame
        tau = self.tau("τ_LC", const_storage)
        return e0 * pm.math.exp(-(delta_k-mu_k0)**2/(2*tau**2))

    def _postprocess(self, delta_k, k0, ax, trace, pmt, actual_x):
        print(trace)
        e0 = float(trace.posterior["E0"].median())
        mu_k0 = float(trace.posterior["mu_LC_k0"].median())
        tau = float(trace.posterior["τ_LC"].median())
        ax.plot(actual_x, e0*np.exp(-(delta_k-mu_k0)**2/(2*tau**2)), **PLOT_STYLES[pmt[0]])


class ExponentialLC(LightCurve):
    E0 = E0_field()
    TAU = tau_field()

    def get_lc(self, delta_k, k0, const_storage):
        tau = self.TAU("τ_LC", const_storage)
        e0 = self.E0("E0", const_storage)
        return e0*pm.math.exp(delta_k/tau)

    def _postprocess(self, delta_k, k0, ax, trace, pmt, actual_x):
        e0 = float(trace.posterior["E0"].median())
        tau = float(trace.posterior["τ_LC"].median())
        ax.plot(actual_x, e0*np.exp(delta_k/tau), **PLOT_STYLES[pmt[0]])


# class TriangularLC(LightCurve):
#     E_max = E0_field()
#     E_start = E0_field()
#     E_end = E0_field()
#     x_0 = DistributionField("normal", mu=0, sigma=1.0)
#
#     def get_lc(self, delta_k, k0):
#         e_max = self.E_max("E_peak")
#         e_start = self.E_start("E_start")
#         e_end = self.E_start("E_end")
#         mu = pm.Deterministic("mu_LC", k0+mu_k0)  # Maximum position from start of frame
#         sigma = self.sigma("sigma_LC")
#         return e0 * pm.math.exp(-(delta_k-mu_k0)**2/(2*sigma**2))


class LC_Alter(AlternatingNode):
    DISPLAY_NAME = "LC"
    SEL__const = ConstantLC().generate_subform()
    SEL__linear = LinearLC().generate_subform()
    SEL__gauss = GaussianLC().generate_subform()
    SEL__exp = ExponentialLC().generate_subform()