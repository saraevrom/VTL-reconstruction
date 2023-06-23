import numba as nb
import numpy as np
import pymc as pm

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
            if (last - first) < 20:
                triggered = False
            else:
                break
    return first, last


def track_threshold(track_data, threshold):
    print("THRESH =", threshold)
    amplitudes = np.max(track_data, axis=0)
    tracks_matrix = amplitudes > threshold

    tracks = (track_data > threshold)
    tracks = np.logical_or.reduce(tracks, axis=(1, 2))
    # print(tracks)
    k_start, k_end = find_track(tracks)
    print("ASSUMING TRACK IN", k_start, "--", k_end)
    assert k_end > k_start >= 0
    return k_start, k_end

# coordinates of the pixel's center (from -3.5 till +3.5 in pixel size)
Xp = np.arange(-3.5, 4.5, 1.)
Yp = np.arange(-3.5, 4.5,  1.)

# pytensor EnsquaredEnergy function
def d_erf(dX, scale, pixel_size=1.0):
    return pm.math.erf((dX + 0.5*pixel_size) / scale) - pm.math.erf((dX - 0.5*pixel_size) / scale)

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

def ensquared_energy(X, Y, sigmaPSF, T):
    scale = sigmaPSF * np.sqrt(2.)
    x0, y0 = track_delta_matrix(T)
    a = d_erf(x0-X, scale) * d_erf(y0-Y, scale)
    return 0.25 * a



def ensquared_energy_avg(X, Y, dX, dY, sigmaPSF, T):
    return 0.2 * (ensquared_energy(X, Y, sigmaPSF, T) +
                  ensquared_energy(X - 0.2 * dX, Y - 0.2 * dY, sigmaPSF, T) +
                  ensquared_energy(X - 0.4 * dX, Y - 0.4 * dY, sigmaPSF, T) +
                  ensquared_energy(X + 0.2 * dX, Y + 0.2 * dY, sigmaPSF, T) +
                  ensquared_energy(X + 0.4 * dX, Y + 0.4 * dY, sigmaPSF, T))

# def ensquared_energy_avg(X, Y, dX, dY, sigmaPSF, T):
#     return 0.2 * (ensquared_energy(X, Y, sigmaPSF, T) +
#                   ensquared_energy(X - 0.2 * dX, Y - 0.2 * dY, sigmaPSF, T) +
#                   ensquared_energy(X - 0.4 * dX, Y - 0.4 * dY, sigmaPSF, T) +
#                   ensquared_energy(X - 0.6 * dX, Y - 0.6 * dY, sigmaPSF, T) +
#                   ensquared_energy(X - 0.8 * dX, Y - 0.8 * dY, sigmaPSF, T))