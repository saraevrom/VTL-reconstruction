import numba as nb
import numpy as np
import pymc as pm

from vtl_common.parameters import PIXEL_SIZE
from vtl_common.localized_GUI.plotter import LOWER_EDGES
PIXEL_POSITIONS = LOWER_EDGES+PIXEL_SIZE/2

@nb.njit(cache=True)
def create_coord_mesh(T, t_off=0):
    result_x = np.zeros(shape=(T, 16, 16))
    result_y = np.zeros(shape=(T, 16, 16))
    result_t = np.zeros(shape=(T, 16, 16))

    for k in range(T):
        for i in range(16):
            for j in range(16):
                result_x[k,i,j] = PIXEL_POSITIONS[i]
                result_y[k,i,j] = PIXEL_POSITIONS[j]
                result_t[k,i,j] = k+t_off
    return result_x, result_y, result_t

@nb.njit(cache=True)
def create_temporal_part(T, t_start, src):
    result_t = np.zeros(shape=(T, 16, 16))

    for k in range(T):
        for i in range(16):
            for j in range(16):
                result_t[k,i,j] = src[k+t_start]
    return result_t

def d_erf(dX, scale, pixel_size=1.0):
    return pm.math.erf((dX + 0.5*pixel_size) / scale) - pm.math.erf((dX - 0.5*pixel_size) / scale)

def ensquared_energy_full(x_mesh, y_mesh, x0, y0, psf):
    scale = np.sqrt(2)*psf
    a = d_erf(x0 - x_mesh, scale, pixel_size=PIXEL_SIZE) * d_erf(y0 - y_mesh, scale, pixel_size=PIXEL_SIZE)/4
    return a

def ensquared_energy_avg(x, y, dx, dy, x0, y0, psf, steps=5):
    if steps==1:
        return ensquared_energy_full(x, y, x0, y0, psf)
    s = 0
    r = 0

    ds = np.linspace(-0.5,0.5,steps)

    for d in ds:
        s += ensquared_energy_full(x, y, x0 + d*dx, y0 + d*dy, psf)
        r += 1
    return s/r
