import pymc as pm
import numpy as np
import pytensor.tensor as pt
import numba as nb
from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
from .common import ensquared_energy_avg, track_threshold


def linear_track_model(track_points, minv=0.01, maxv = 0.45, min_e=0.0001, max_e=1.0):
    with pm.Model() as model:
        #k_start, k_end = track_threshold(track_points, threshold)
        k_start = 0
        k_end = track_points.shape[0]

        x0 = pm.Uniform('X0', -4.5, 4.5)
        y0 = pm.Uniform('Y0', -4.5, 4.5)
        phi0 = pm.Uniform('Phi0', -np.pi, np.pi)
        u0 = pm.Uniform('U0', minv, maxv)
        e0 = pm.Uniform('E0', min_e, max_e)
        #e0 = pm.HalfNormal('E0', 1.)

        sigma0 = pm.HalfNormal('Sigma0', 1.)
        sigmaPSF = pm.HalfNormal('SigmaPSF', 1.)

        kk = np.arange(k_start, k_end)
        k0 = (k_start+k_end)/2
        X = x0 + u0 * pm.math.cos(phi0) * (kk-k0)
        Y = y0 + u0 * pm.math.sin(phi0) * (kk-k0)
        dX = u0 * pm.math.cos(phi0) * 1.
        dY = u0 * pm.math.sin(phi0) * 1.

        X = pt.expand_dims(X,(1,2))
        Y = pt.expand_dims(Y,(1,2))

        mu = e0 * ensquared_energy_avg(X, Y, dX, dY, sigmaPSF, k_end - k_start)
        A = pm.Normal('A', mu=mu, sigma=sigma0,
                      observed=track_points[k_start:k_end], shape=(k_end - k_start, 8, 8))  # A = A[k,i,j]
    return model


class LinearTrackModelForm(FormNode):
    USE_SCROLLVIEW = False
    DISPLAY_NAME = get_locale("reconstruction_model.linear")
    #FIELD__threshold = create_value_field(FloatNode, get_locale("reconstruction_model.linear.threshold"), 6.0)
    FIELD__minv = create_value_field(FloatNode, get_locale("reconstruction_model.linear.minv"), 0.01)
    FIELD__maxv = create_value_field(FloatNode, get_locale("reconstruction_model.linear.maxv"), 0.45)
    FIELD__min_e = create_value_field(FloatNode, get_locale("reconstruction_model.linear.min_e"), 10.0)
    FIELD__max_e = create_value_field(FloatNode, get_locale("reconstruction_model.linear.max_e"), 60.0)

    def get_data(self):
        data = super().get_data()
        return lambda x: linear_track_model(x, **data)




