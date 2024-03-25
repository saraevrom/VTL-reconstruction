import numpy as np
import pymc as pm

from .model_base import ModelWithParameters, ReconstructionModelWrapper
from .form_prototypes import DistributionField, PassthroughField, IntField
from tools.new_reconstructor.models.light_curves import create_lc_alter
from vtl_common.parameters import PIXEL_SIZE, HALF_GAP_SIZE, HALF_PIXELS
from common_functions import create_coord_mesh, ensquared_energy_avg

SPAN = HALF_GAP_SIZE + HALF_PIXELS*PIXEL_SIZE


def resample(src:np.ndarray,change_mul):
    old_len = src.shape[0]
    new_len = int(old_len*change_mul)
    sampler = np.linspace(0, old_len-1, new_len)
    #print(sampler)
    return np.interp(sampler, np.arange(old_len), src)

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

class BaseLinearPlanarTrackModel(ReconstructionModelWrapper):

    # Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    # One more static field is already included in ReconstructionModelWrapper:
    # It gets transformed into pymc.<Distribution>. And it automatically adds error parameters (sigma, nu, etc...).

    # final_distribution = FinalDistributionField()

    # accel = DistributionField("const", 0.0)
    #U0 = DistributionField("uniform", lower=0.05, upper=0.5)
    #Phi0 = DistributionField("uniform", lower=-180.0, upper=180.0)

    SigmaPSF = DistributionField("exponential", lam={"selection_type": "lambda", "value":1.0})
    SigmaCoeff = DistributionField("const", 0.0)
    LC = PassthroughField(create_lc_alter)
    ee_steps = IntField(default_value=5)

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

            # phi0_deg = self.Phi0("Phi0", consts)
            # phi0 = phi0_deg * np.pi / 180.0
            #u0 = self.U0("U0", consts)
            #a = self.accel("accel", consts)

            # sigmaPSF = pm.HalfNormal('SigmaPSF', 1.) * PIXEL_SIZE
            sigmaPSF = self.SigmaPSF('SigmaPSF',consts) * PIXEL_SIZE
            sigmaCOEFF = self.SigmaCoeff('SigmaCoeff',consts) 

            #u = (u0 + a * delta_k) * PIXEL_SIZE
            #t_params = self.get_track_parameters(consts)

            #u = self.get_speed(delta_k,t_params)*PIXEL_SIZE
            lc = self.LC.get_lc(delta_k, k0, consts)

            #u_int = (u0 * delta_k + a * delta_k ** 2 / 2)*PIXEL_SIZE
            #u_int = self.get_displacement(delta_k,t_params)*PIXEL_SIZE

            orientation_form = reconstructor_main.orientation_form
            orientation = orientation_form.get_values()


            X,Y,dX,dY = self.get_kinematics(consts,delta_k,x0,y0,orientation)

            centerdist = (X**2+Y**2)**0.5
            psf = sigmaPSF+sigmaCOEFF*centerdist

            # X = x0 + u_int * pm.math.cos(phi0)
            # Y = y0 + u_int * pm.math.sin(phi0)
            # dX = u * pm.math.cos(phi0)
            # dY = u * pm.math.sin(phi0)

            intensity = lc * ensquared_energy_avg(pixel_xs, pixel_ys, dX, dY, X, Y, psf,steps=self.ee_steps)
            # intensity = lc * ensquared_energy_avg(pixel_xs, pixel_ys, dX, dY, X, Y, sigmaPSF)
            obs = observed[0][k_start:k_end]
            observed_var = self.final_distribution('OBSERVED', mu=intensity,
                                                   observed=obs[:, alive])
        return ModelWithParameters(model, {
            "k0": k0,
            "k_start": k_start,
            "k_end": k_end,
        }, consts)

    def get_kinematics(self,consts,delta_k,x0,y0,orientation):
        raise NotImplementedError

    def get_overlay_two_points(self,model_params: ModelWithParameters):
        raise NotImplementedError

    # def get_track_parameters(self, consts):
    #     raise NotImplementedError
    #
    # def get_displacement(self,delta_k,t_params):
    #     raise NotImplementedError
    #
    # def get_speed(self,delta_k,t_params):
    #     raise NotImplementedError

    def ask_k0(self, model_params: ModelWithParameters):
        return model_params.parameters["k0"]

    def ask_x0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("X0")

    def ask_y0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("Y0")

    def ask_sigma_psf(self, model_params:ModelWithParameters):
        return model_params.get_estimation('SigmaPSF')

    def reconstruction_overlay(self, plotter, model_params: ModelWithParameters):
        print("DRAW_CALL")
        params = model_params.parameters
        r = model_params.get_estimation("SigmaPSF")*PIXEL_SIZE
        x0 = model_params.get_estimation("X0")
        y0 = model_params.get_estimation("Y0")
        k0 = params["k0"]
        # params = model_params.parameters
        # x0 = model_params.get_estimation("X0")
        # y0 = model_params.get_estimation("Y0")
        # phi = model_params.get_estimation("Phi0")*np.pi/180
        # u0 = model_params.get_estimation("U0")
        # a = model_params.get_estimation("accel")
        #
        #
        # k0 = params["k0"]
        # k_start = params["k_start"]
        # k_end = params["k_end"]
        #
        # u_int_start = (u0 * (k_start-k0) + a * (k_start-k0) ** 2 / 2) * PIXEL_SIZE
        # u_int_end = (u0 * (k_end - k0) + a * (k_end - k0) ** 2 / 2) * PIXEL_SIZE
        #
        #
        # x_start = x0+np.cos(phi)*u_int_start
        # y_start = y0+np.sin(phi)*u_int_start
        # x_end = x0+np.cos(phi)*u_int_end
        # y_end = y0+np.sin(phi)*u_int_end
        x_start,x_end, y_start,y_end = self.get_overlay_two_points(model_params)

        #x_start, y_start = rect_raycast(x0, y0, phi+np.pi, HALF_PIXELS/2*PIXEL_SIZE)
        #x_end, y_end = rect_raycast(x0, y0, phi, HALF_PIXELS/2*PIXEL_SIZE)
        dx = x_end-x_start
        dy = y_end-y_start

        plotter.plot_arrow(x_start, y_start, dx, dy, color="red", width=r, length_includes_head=True)
        plotter.set_origin(x0, y0, k0)

    def postprocess(self, ax, model_params: ModelWithParameters, actual_x):
        trace = model_params.idata
        k_start = model_params.parameters["k_start"]
        k_end = model_params.parameters["k_end"]
        kk = np.arange(k_start, k_end)
        k0 = model_params.parameters["k0"]
        delta_k = kk - k0
        actual_x1 = resample(actual_x[kk], 10)
        delta_k1 = resample(delta_k,10)
        return self.LC.postprocess_plot(delta_k1, k0, ax, model_params, model_params.pmt, actual_x=actual_x1)
