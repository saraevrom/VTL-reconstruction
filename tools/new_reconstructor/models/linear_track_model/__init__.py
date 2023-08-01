import numpy as np
import pymc as pm
import pytensor.tensor as pt

from ..model_base import ModelWithParameters, ReconstructionModelWrapper
from ..form_prototypes import DistributionField, PassthroughField
from .light_curves import create_lc_alter
from vtl_common.parameters import PIXEL_SIZE, HALF_GAP_SIZE, HALF_PIXELS
from common_functions import create_coord_mesh, ensquared_energy_avg

SPAN = HALF_GAP_SIZE + HALF_PIXELS*PIXEL_SIZE


def x0_from_pmt(pmt):
    low = HALF_GAP_SIZE
    high = -HALF_GAP_SIZE
    if ("A" in pmt) or ("C" in pmt):
        low = -SPAN
    if ("B" in pmt) or ("D" in pmt):
        high = SPAN
    assert high >= low
    return pm.Uniform('X0', low, high)

def y0_from_pmt(pmt):
    low = HALF_GAP_SIZE
    high = -HALF_GAP_SIZE
    if ("C" in pmt) or ("D" in pmt):
        low = -SPAN
    if ("A" in pmt) or ("B" in pmt):
        high = SPAN
    assert high >= low
    return pm.Uniform('Y0', low, high)

class LinearTrackModel(ReconstructionModelWrapper):

    # Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    # One more static field is already included in ReconstructionModelWrapper:
    # It gets transformed into pymc.<Distribution>. And it automatically adds error parameters (sigma, nu, etc...).

    # final_distribution = FinalDistributionField()

    accel = DistributionField("const", 0.0)
    U0 = DistributionField("uniform", lower=0.05, upper=0.5)
    Phi0 = DistributionField("uniform", lower=-180.0, upper=180.0)
    LC = PassthroughField(create_lc_alter)

    def generate_pymc_model(self, observed, cut_start, cut_end, broken, pmt, reconstructor_main) -> ModelWithParameters:
        with pm.Model() as model:
            consts = dict()
            k_start = cut_start
            k_end = cut_end
            T = k_end - k_start
            k0 = (k_start + k_end) / 2

            mesh_x, mesh_y, mesh_t = create_coord_mesh(T, k_start)
            alive = np.logical_not(broken)
            pixel_xs = mesh_x[:, alive]
            pixel_ys = mesh_y[:, alive]
            pixel_ts = mesh_t[:, alive]

            delta_k = pixel_ts - k0

            x0 = x0_from_pmt(pmt)
            y0 = y0_from_pmt(pmt)

            phi0_deg = self.Phi0("Phi0", consts)
            phi0 = phi0_deg * np.pi / 180.0
            u0 = self.U0("U0", consts)
            a = self.accel("accel", consts)

            sigmaPSF = pm.HalfNormal('SigmaPSF', 1.) * PIXEL_SIZE

            u = (u0 + a * delta_k) * PIXEL_SIZE
            lc = self.LC.get_lc(delta_k, k0, consts)

            u_int = (u0 * delta_k + a * delta_k ** 2 / 2)*PIXEL_SIZE
            X = x0 + u_int * pm.math.cos(phi0)
            Y = y0 + u_int * pm.math.sin(phi0)
            dX = u * pm.math.cos(phi0)
            dY = u * pm.math.sin(phi0)
            intensity = lc * ensquared_energy_avg(pixel_xs, pixel_ys, dX, dY, X, Y, sigmaPSF)
            obs = observed[0][k_start:k_end]
            observed_var = self.final_distribution('OBSERVED', mu=intensity,
                                                   observed=obs[:, alive])
        return ModelWithParameters(model, {
            "k0": k0,
            "k_start": k_start,
            "k_end": k_end,
        }, consts)

    def reconstruction_overlay(self, plotter, model_params: ModelWithParameters):
        print("DRAW_CALL")
        i_trace = model_params.idata
        params = model_params.parameters
        post = i_trace.posterior
        x0 = model_params.get_estimation("X0")
        y0 = model_params.get_estimation("Y0")
        r = model_params.get_estimation("SigmaPSF")*PIXEL_SIZE
        phi = model_params.get_estimation("Phi0")*np.pi/180
        u0 = model_params.get_estimation("U0")
        a = model_params.get_estimation("accel")


        k0 = params["k0"]
        k_start = params["k_start"]
        k_end = params["k_end"]

        u_int_start = (u0 * (k_start-k0) + a * (k_start-k0) ** 2 / 2) * PIXEL_SIZE
        u_int_end = (u0 * (k_end - k0) + a * (k_end - k0) ** 2 / 2) * PIXEL_SIZE


        x_start = x0+np.cos(phi)*u_int_start
        y_start = y0+np.sin(phi)*u_int_start
        x_end = x0+np.cos(phi)*u_int_end
        y_end = y0+np.sin(phi)*u_int_end

        #x_start, y_start = rect_raycast(x0, y0, phi+np.pi, HALF_PIXELS/2*PIXEL_SIZE)
        #x_end, y_end = rect_raycast(x0, y0, phi, HALF_PIXELS/2*PIXEL_SIZE)
        dx = x_end-x_start
        dy = y_end-y_start

        plotter.plot_arrow(x_start, y_start, dx, dy, color="red", width=r, length_includes_head=True)
        plotter.set_origin(x0, y0, k0)

    def postprocess(self, ax, model_params: ModelWithParameters):
        trace = model_params.idata
        k_start = model_params.parameters["k_start"]
        k_end = model_params.parameters["k_end"]
        kk = np.arange(k_start, k_end)
        k0 = model_params.parameters["k0"]
        delta_k = kk - k0
        self.LC.postprocess_plot(delta_k, k0, ax, model_params, model_params.pmt)

LINEAR_TRACK_FORM = LinearTrackModel()
