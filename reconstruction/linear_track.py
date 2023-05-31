import pymc as pm
import numpy as np
import pytensor.tensor as pt

from .common import ensquared_energy_avg
from .base_model import ReconstructionModelWrapper
from .form_prototypes import FloatField, DistributionField
from vtl_common.parameters import PIXEL_SIZE, HALF_PIXELS



def rect_raycast(x0, y0, phi, edge):
    dx = np.cos(phi)
    dy = np.sin(phi)
    if dx>0:
        tx = (edge-x0)/dx
    else:
        tx = (-edge-x0)/dx
    if dy > 0:
        ty = (edge - y0) / dy
    else:
        ty = (-edge - y0) / dy
    t = min(tx,ty)
    return x0+t*dx, y0+t*dy


class LinearTrackModel(ReconstructionModelWrapper):

    # Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    # One more static field is already included in ReconstructionModelWrapper:
    # It gets transformed into pymc.<Distribution>. And it automatically adds error parameters (sigma, nu, etc...).

    # final_distribution = FinalDistributionField()

    U0 = DistributionField("uniform", lower=0.05, upper=0.5)
    E0 = DistributionField("uniform", lower=10.0, upper=60.0)
    Phi0 = DistributionField("uniform", lower=-180.0, upper=180.0)
    # min_v = FloatField(0.01)   # Will be transformed to float
    # max_v = FloatField(0.45)
    # min_e = FloatField(10.0)
    # max_e = FloatField(60.0)
    # min_phi = FloatField(-180)
    # max_phi = FloatField(180)

    def get_pymc_model(self,observed):
        with pm.Model() as model:
            # k_start, k_end = track_threshold(track_points, threshold)
            k_start = 0
            k_end = observed.shape[0]

            x0 = pm.Uniform('X0', -4.5, 4.5)
            y0 = pm.Uniform('Y0', -4.5, 4.5)

            phi0_deg = self.Phi0("Phi0")
            #phi0_deg = pm.Uniform('Phi0', self.min_phi, self.max_phi)
            phi0 = phi0_deg * np.pi/180.0
            # print(self.min_v, self.max_v)
            # u0 = pm.Uniform('U0', self.min_v, self.max_v)
            # e0 = pm.Uniform('E0', self.min_e, self.max_e)
            u0 = self.U0("U0")
            e0 = self.E0("E0")

            #sigma0 = pm.HalfNormal('Sigma0', 1.)
            sigmaPSF = pm.HalfNormal('SigmaPSF', 1.)

            kk = np.arange(k_start, k_end)
            k0 = (k_start + k_end) / 2
            X = x0 + u0 * pm.math.cos(phi0) * (kk - k0)
            Y = y0 + u0 * pm.math.sin(phi0) * (kk - k0)
            dX = u0 * pm.math.cos(phi0) * 1.
            dY = u0 * pm.math.sin(phi0) * 1.

            X = pt.expand_dims(X, (1, 2))
            Y = pt.expand_dims(Y, (1, 2))

            mu = e0 * ensquared_energy_avg(X, Y, dX, dY, sigmaPSF, k_end - k_start)
            A = self.final_distribution('A', mu=mu,
                          observed=observed[k_start:k_end], shape=(k_end - k_start, 8, 8))  # A = A[k,i,j]
        return model

    def reconstruction_overlay(self, plotter, i_trace, offset):
        print("DRAW_CALL")
        x_off, y_off = offset
        post = i_trace.posterior
        x0 = post["X0"].median()*PIXEL_SIZE
        y0 = post["Y0"].median()*PIXEL_SIZE
        r = post["SigmaPSF"].median()
        phi = post["Phi0"].median()*np.pi/180
        u0 = post["U0"].median()

        x_start, y_start = rect_raycast(x0, y0, phi+np.pi, HALF_PIXELS/2*PIXEL_SIZE)
        x_end, y_end = rect_raycast(x0, y0, phi, HALF_PIXELS/2*PIXEL_SIZE)
        dx = x_end-x_start
        dy = y_end-y_start

        plotter.plot_arrow(x_start+x_off, y_start+y_off, dx, dy, color="red", width=r, length_includes_head=True)


LINEAR_TRACK = LinearTrackModel()
