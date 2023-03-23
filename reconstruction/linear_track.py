import pymc as pm
import numpy as np
import pytensor.tensor as pt
import numba as nb
from vtl_common.common_GUI.tk_forms_assist import FormNode, LabelNode, FloatNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale

# coordinates of the pixel's center (from -3.5 till +3.5 in pixel size)
Xp = np.arange(-3.5, 4.5, 1.)
Yp = np.arange(-3.5, 4.5,  1.)

# pytensor EnsquaredEnergy function
def d_erf(dX, scale):
    return pm.math.erf((dX + 0.5) / scale) - pm.math.erf((dX - 0.5) / scale)


@nb.njit()
def find_track(track_data):
    first = -1
    last = -1
    triggered = False
    for i in range(len(track_data)):
        data = track_data[i]
        if data:
            triggered = True
            if first<0:
                first = i
            last = i
        elif triggered:
            if last==first:
                triggered = False
            else:
                break
    return first, last


@nb.njit()
def track_delta_matrix(T):
    '''
    Returns array with (T,8,8). (Something like meshgrid)
    :param X:
    :param Y:
    :return:
    '''
    Xp = np.arange(-3.5, 4.5, 1.)
    Yp = np.arange(-3.5, 4.5,  1.)
    delta_x = np.zeros((T, Xp.shape[0], Yp.shape[0]))
    delta_y = np.zeros((T, Xp.shape[0], Yp.shape[0]))
    for k in range(T):
        for i in range(8):
            for j in range(8):
                delta_x[k, i, j] = Xp[i]
                delta_y[k, i, j] = Yp[j]
    return delta_x, delta_y

def EnsquaredEnergy(X, Y, sigmaPSF, T):
    scale = sigmaPSF * np.sqrt(2.)
    x0, y0 = track_delta_matrix(T)
    a = d_erf(x0-X, scale) * d_erf(y0-Y, scale)
    # a_matr = []
    # for i in range(8):
    #     a_row = []
    #     for j in range(8):
    #         a_row.append(d_erf(X - Xp[i], scale) * d_erf(Y - Yp[j], scale))
    #     a_matr.append(pt.stack(a_row))
    # a0 = pt.stack(a_matr)
    # a = pt.moveaxis(a0, [0, 1, 2], [1, 2, 0])
    return 0.25 * a


def EnsquaredEnergyAvg(X, Y, dX, dY, sigmaPSF, T):
    return 0.2 * (EnsquaredEnergy(X, Y, sigmaPSF, T) + EnsquaredEnergy(X - 0.2 * dX, Y - 0.2 * dY,
                                                                    sigmaPSF, T) + EnsquaredEnergy(X - 0.4 * dX,
                                                                                                Y - 0.4 * dY,
                                                                                                sigmaPSF, T) + EnsquaredEnergy(
        X - 0.6 * dX, Y - 0.6 * dY, sigmaPSF, T) + EnsquaredEnergy(X - 0.8 * dX, Y - 0.8 * dY, sigmaPSF, T))

# class LinearTrackModel(pm.Model):
#     def __init__(self, track_points):
#         super().__init__(name="")
#
#         stds = np.std(track_points, axis=0)
#         tracks = track_points > (stds*3)
#         tracks = np.logical_or.reduce(tracks, axis=(1, 2))
#         k_start, k_end = find_track(tracks)
#         assert k_end>k_start>=0
#
#         x0 = pm.Uniform('X0', -4.5, 4.5)
#         y0 = pm.Uniform('Y0', -4.5, 4.5)
#         phi0 = pm.Uniform('Phi0', 0, 2 * np.pi)
#         u0 = pm.Uniform('U0', 0.05, 0.45)
#         e0 = pm.Uniform('E0', 0.5, 1.5)
#
#
#         sigma0 = pm.HalfNormal('Sigma0', 1.)
#         sigmaPSF = pm.HalfNormal('SigmaPSF', 1.)
#
#
#         kk = np.arange(k_start, k_end)
#         X = x0 + u0 * pm.math.cos(phi0) * (kk)
#         Y = y0 + u0 * pm.math.sin(phi0) * (kk)
#         dX = u0 * pm.math.cos(phi0) * 1.
#         dY = u0 * pm.math.sin(phi0) * 1.
#
#         #envelope = (pt.math.switch(start_t < kk , 0, 1) * pt.math.switch(end_t > kk , 0, 1)).flatten()
#         mu = e0 * EnsquaredEnergyAvg(X, Y, dX, dY, sigmaPSF)#*envelope
#         A = pm.Normal('A', mu=mu, sigma=sigma0,
#                       observed=track_points[k_start:k_end], shape=(k_end-k_start, 8, 8))  # A = A[k,i,j]

def linear_track_model(track_points, threshold=6.0):
    with pm.Model() as model:
        stds = np.std(track_points, axis=0)
        tracks = track_points > threshold
        tracks = np.logical_or.reduce(tracks, axis=(1, 2))
        #print(tracks)
        k_start, k_end = find_track(tracks)
        print("ASSUMING TRACK IN", k_start, "--", k_end)
        assert k_end > k_start >= 0

        x0 = pm.Uniform('X0', -4.5, 4.5)
        y0 = pm.Uniform('Y0', -4.5, 4.5)
        phi0 = pm.Uniform('Phi0', -np.pi, np.pi)
        u0 = pm.Uniform('U0', 0.05, 0.45)
        e0 = pm.Uniform('E0', 10.5, 60.0)
        #e0 = pm.HalfNormal('E0', 1.)

        sigma0 = pm.HalfNormal('Sigma0', 1.)
        sigmaPSF = pm.HalfNormal('SigmaPSF', 1.)

        kk = np.arange(k_start, k_end)
        k0 =(k_start+k_end)/2
        X = x0 + u0 * pm.math.cos(phi0) * (kk-k0)
        Y = y0 + u0 * pm.math.sin(phi0) * (kk-k0)
        dX = u0 * pm.math.cos(phi0) * 1.
        dY = u0 * pm.math.sin(phi0) * 1.

        X = pt.expand_dims(X,(1,2))
        Y = pt.expand_dims(Y,(1,2))

        # envelope = (pt.math.switch(start_t < kk , 0, 1) * pt.math.switch(end_t > kk , 0, 1)).flatten()
        mu = e0 * EnsquaredEnergyAvg(X, Y, dX, dY, sigmaPSF, k_end-k_start)  # *envelope
        A = pm.Normal('A', mu=mu, sigma=sigma0,
                      observed=track_points[k_start:k_end], shape=(k_end - k_start, 8, 8))  # A = A[k,i,j]
    return model


class LinearTrackModelForm(FormNode):
    USE_SCROLLVIEW = False
    DISPLAY_NAME = get_locale("reconstruction_model.linear")
    FIELD__threshold = create_value_field(FloatNode, get_locale("reconstruction_model.linear.threshold"), 6.0)

    def get_data(self):
        data = super().get_data()
        return lambda x: linear_track_model(x, **data)

def linear_track_model_alt(track_points):
    with pm.Model() as model:
        # image parameters
        # e0 = pm.Uniform('E0', 0.0, np.max(track_points))
        e0 = pm.HalfNormal('E0', 1.)
        sigma0 = pm.HalfNormal('Sigma0', 1.)
        # b0 = pm.HalfNormal('b_0', 1.)
        sigma_psf = pm.HalfNormal('SigmaPSF', 1.)
        front_mul = pm.HalfNormal("FrontMul", 10.0)

        # start point
        x1 = pm.Uniform('X1', -4.5, 4.5)
        y1 = pm.Uniform('Y1', -4.5, 4.5)
        # end point
        x2 = pm.Uniform('X2', -4.5, 4.5)
        y2 = pm.Uniform('Y2', -4.5, 4.5)

        T = track_points.shape[0]
        t1 = pm.Uniform('T1', 0.0, 1.0*T)  # start time
        t2 = pm.Uniform('T2', t1, 1.0*T)   # end time

        u0x = pm.Deterministic("u0_x", (x2 - x1)/(t2 - t1))
        u0y = pm.Deterministic("u0_y", (y2 - y1)/(t2 - t1))
        u0 = pm.Deterministic("u0", (u0x**2 + u0y**2)**0.5)

        ts = np.arange(T)

        x = x1 + u0x * (ts-t1)
        y = y1 + u0y * (ts-t1)

        x = pt.expand_dims(x, (1, 2))
        y = pt.expand_dims(y, (1, 2))

        envelope = pm.math.invlogit(front_mul*(ts-t1)) * pm.math.invlogit(front_mul*(t2-ts))
        #envelope = (pt.math.switch(ts > t1, 0, 1) * pt.math.switch(ts < t2, 0, 1))
        envelope = pt.expand_dims(envelope, (1, 2))
        mu = e0 * EnsquaredEnergyAvg(x, y, u0x, u0y, sigma_psf, T) * envelope
        A = pm.Normal('A', mu=mu, sigma=sigma0,
                       observed=track_points, shape=(T, 8, 8))
        # A = pm.Laplace('A', mu=mu, b=b0,
        #                observed=track_points, shape=(T, 8, 8))

    return model

class LinearTrackAltModelForm(LabelNode):
    DISPLAY_NAME = get_locale("reconstruction_model.linear_alt")

    def get_data(self):
        return linear_track_model_alt
