import pymc as pm
import numpy as np
from ..base_flat_model import BaseLinearPlanarTrackModel
from ..form_prototypes import DistributionField
from ..model_base import ModelWithParameters
from common_functions.hor_to_dev import hor_to_dev
from vtl_common.parameters import PIXEL_SIZE, HALF_GAP_SIZE, HALF_PIXELS

class LinearTrackModel(BaseLinearPlanarTrackModel):
    accel = DistributionField("const", 0.0)
    U0 = DistributionField("uniform", lower=0.05, upper=0.5)
    Phi0 = DistributionField("uniform", lower=-180.0, upper=180.0)

    def get_kinematics(self,consts,delta_k,x0,y0,orientation):
        u0 = self.U0("U0", consts)
        a = self.accel("accel", consts)
        phi0_deg = self.Phi0("Phi0", consts)
        phi0 = phi0_deg * np.pi / 180.0
        u = (u0 + a * delta_k)*PIXEL_SIZE
        u_int = (u0*delta_k + a*delta_k**2/2)*PIXEL_SIZE

        X = x0 + u_int * pm.math.cos(phi0)
        Y = y0 + u_int * pm.math.sin(phi0)
        dX = u * pm.math.cos(phi0)
        dY = u * pm.math.sin(phi0)
        return X,Y,dX,dY

    def get_overlay_two_points(self,model_params: ModelWithParameters):
        params = model_params.parameters
        x0 = model_params.get_estimation("X0")
        y0 = model_params.get_estimation("Y0")
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
        return x_start,x_end,y_start,y_end

    def ask_nu(selfself, model_params: ModelWithParameters):
        return 0.0

    def ask_a(self, model_params: ModelWithParameters):
        return model_params.get_estimation("accel")

    def ask_u0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("U0")

    def ask_phi0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("Phi0")


LINEAR_TRACK_FORM = LinearTrackModel()