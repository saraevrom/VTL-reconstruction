import pymc as pm
import numpy as np
import pytensor.tensor as pt
from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode, ComboNode, BoolNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
from .common import ensquared_energy_avg, track_threshold

def heaviside(x):
    return (pt.sgn(x) + 1) / 2


def sqrt_sigma(x):
    return x/pt.sqrt(1+x**2)

def lc_profile_invlogit(t, t1, t2, w):
    return pm.math.invlogit(w*(t-t1)) * pm.math.invlogit(w*(t2-t))



def lc_profile_sqrt(t, t1, t2, w):
    return sqrt_sigma(w*(t-t1)) * sqrt_sigma(w*(t2-t))



def lc_profile_heaviside(t, t1, t2):
    return heaviside(t-t1)*heaviside(t2-t)






def linear_track_model_alt(track_points, front_mul_mean=10.0, envelope_data=None, discrete_time=False, vmax=1.0,
                           threshold=6.0):


    if envelope_data is None:
        envelope_func = lc_profile_invlogit
        envelope_uses_width = True
    else:
        envelope_func, envelope_uses_width = envelope_data

    k_start, k_end = track_threshold(track_points, threshold)
    k_sigma = (k_end - k_start)/10.0
    k_sigma = max(k_sigma, 1.0)

    with pm.Model() as model:
        # image parameters
        e0 = pm.HalfNormal('E0', 1.)
        sigma0 = pm.HalfNormal('Sigma0', 1.)
        # b0 = pm.HalfNormal('b_0', 1.)
        sigma_psf = pm.HalfNormal('SigmaPSF', 1.)
        #front_mul = pm.Uniform("FrontMul", 0.1, 10.0)
        #front_mul = 100.0

        T = track_points.shape[0]

        if discrete_time:
            t1 = pm.DiscreteUniform('T1', 0, T)  # start time
            t2 = pm.DiscreteUniform('T2', t1, T)  # end time
        else:
            t1 = pm.TruncatedNormal("T1", mu=k_start, sigma=k_sigma, lower=0.0, upper=1.0*T)
            t2 = pm.TruncatedNormal("T2", mu=k_end, sigma=k_sigma, lower=t1, upper=1.0*T)

        #max_distance = delta_t*vmax

        # start point
        x1 = pm.Uniform('X1', -4.5, 4.5)
        y1 = pm.Uniform('Y1', -4.5, 4.5)
        # end point
        x2 = pm.Uniform('X2', -4.5, 4.5)
        y2 = pm.Uniform('Y2', -4.5, 4.5)

        #x2 = pm.Truncated("X2", lower=-4.5, upper=4.5, dist=pm.Uniform.dist(x1 - max_distance, x1 + max_distance))
        #y2 = pm.Truncated("Y2", lower=-4.5, upper=4.5, dist=pm.Uniform.dist(y1 - max_distance, y1 + max_distance))

        u0x = pm.Deterministic("u0_x", (x2 - x1)/(t2 - t1))
        u0y = pm.Deterministic("u0_y", (y2 - y1)/(t2 - t1))
        u0 = pm.Deterministic("u0", (u0x**2 + u0y**2)**0.5)

        ts = np.arange(T)

        x = x1 + u0x * (ts-t1)
        y = y1 + u0y * (ts-t1)

        x = pt.expand_dims(x, (1, 2))
        y = pt.expand_dims(y, (1, 2))

        #envelope = pm.math.invlogit(front_mul*(ts-t1)) * pm.math.invlogit(front_mul*(t2-ts))
        if envelope_uses_width:
            front_mul = pm.HalfNormal("FrontMul", front_mul_mean)
            profile = envelope_func(ts, t1, t2, front_mul)
        else:
            profile = envelope_func(ts, t1, t2)
        #envelope = (pt.math.switch(ts > t1, 0, 1) * pt.math.switch(ts < t2, 0, 1))
        profile = pt.expand_dims(profile, (1, 2))
        mu = e0 * profile * ensquared_energy_avg(x, y, u0x, u0y, sigma_psf, T)
        A = pm.Normal('A', mu=mu, sigma=sigma0,
                       observed=track_points, shape=(T, 8, 8))
        # A = pm.Laplace('A', mu=mu, b=b0,
        #                observed=track_points, shape=(T, 8, 8))

    return model

class EnvelopeSelector(ComboNode):
    SELECTION_READONLY = True
    DEFAULT_VALUE = "invlogit"
    VALUES = ["invlogit", "heaviside", "sqrt"]
    DISPLAY_NAME = get_locale("reconstruction_model.linear_alt.envelope")

    def get_data(self):
        data = super().get_data()
        if data == "invlogit":
            return lc_profile_invlogit,True
        elif data == "sqrt":
            return lc_profile_sqrt, True
        elif data == "heaviside":
            return lc_profile_heaviside,False

        return lc_profile_invlogit,True

class LinearTrackHardModelForm(FormNode):
    DISPLAY_NAME = get_locale("reconstruction_model.linear_alt")
    FIELD__front_mul_mean = create_value_field(FloatNode, get_locale("reconstruction_model.linear_alt.front_mul"), 10.0)
    FIELD__envelope_data = EnvelopeSelector
    FIELD__discrete_time = create_value_field(BoolNode, get_locale("reconstruction_model.linear_alt.discrete_time"),
                                              False)

    def get_data(self):
        data = super().get_data()
        return lambda x: linear_track_model_alt(x, **data)
