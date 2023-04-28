import pymc as pm
import numpy as np
import pytensor.tensor as pt

from .common import ensquared_energy_avg
from .base_model import ReconstructionModelWrapper
from .form_prototypes import FloatField

class LinearTrackModel(ReconstructionModelWrapper):
    '''
    Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    '''
    # One more static field is already included in ReconstructionModelWrapper:
    # It gets transformed into pymc.<Distribution>. And it automatically adds error parameters (sigma, nu, etc...).

    # final_distribution = FinalDistributionField()
    min_v = FloatField(0.01)   # Will be transformed to float
    max_v = FloatField(0.45)
    min_e = FloatField(10.0)
    max_e = FloatField(60.0)

    def get_pymc_model(self,observed):
        with pm.Model() as model:
            # k_start, k_end = track_threshold(track_points, threshold)
            k_start = 0
            k_end = observed.shape[0]

            x0 = pm.Uniform('X0', -4.5, 4.5)
            y0 = pm.Uniform('Y0', -4.5, 4.5)
            phi0_deg = pm.Uniform('Phi0', -180.0, 180.0)
            phi0 = phi0_deg * np.pi/180.0
            print(self.min_v, self.max_v)
            u0 = pm.Uniform('U0', self.min_v, self.max_v)
            e0 = pm.Uniform('E0', self.min_e, self.max_e)

            sigma0 = pm.HalfNormal('Sigma0', 1.)
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
            A = self.final_distribution('A', mu=mu, sigma=sigma0,
                          observed=observed[k_start:k_end], shape=(k_end - k_start, 8, 8))  # A = A[k,i,j]
        return model


LINEAR_TRACK = LinearTrackModel()
