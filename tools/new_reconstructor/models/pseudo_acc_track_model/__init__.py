import pymc as pm
import numpy as np
from ..base_flat_model import BaseLinearPlanarTrackModel
from ..form_prototypes import DistributionField, PassthroughField
from ..model_base import ModelWithParameters
from vtl_common.parameters import PIXEL_SIZE, HALF_GAP_SIZE, HALF_PIXELS
from .velocity_alter import create_z_alter

def calculate_position(x0,y0,u0,u_z,phi0,delta_k):
    X = x0 + ( u0 * delta_k * np.cos(phi0)) / (1 + u_z * delta_k)
    Y = y0 + ( u0 * delta_k * np.sin(phi0)) / (1 + u_z * delta_k)
    return X,Y


class LinearTrackModelPseudoAcceleration(BaseLinearPlanarTrackModel):
    #u_z = DistributionField("const", 0.0)
    z_correction = PassthroughField(create_z_alter)
    U0 = DistributionField("uniform", lower=0.05, upper=0.5)
    Phi0 = DistributionField("uniform", lower=-180.0, upper=180.0)

    def get_kinematics(self,consts,delta_k,x0,y0):
        u0 = self.U0("U0", consts)*PIXEL_SIZE
        #u_z = self.u_z("u_z", consts)
        nu = self.z_correction()
        u_z = -nu
        phi0_deg = self.Phi0("Phi0", consts)
        phi0 = phi0_deg * np.pi / 180.0

        X = x0+(u0*delta_k*pm.math.cos(phi0))/(1+u_z*delta_k)
        Y = y0+(u0*delta_k*pm.math.sin(phi0))/(1+u_z*delta_k)
        dX = u0*pm.math.cos(phi0)/(1+u_z*delta_k)**2
        dY = u0*pm.math.sin(phi0)/(1+u_z*delta_k)**2
        #dX =(u0*pm.math.cos(phi0)-x0*u_z)/(1+u_z*delta_k)**2
        #dY =(u0*pm.math.sin(phi0)-y0*u_z)/(1+u_z*delta_k)**2
        return X,Y,dX,dY

    def get_overlay_two_points(self,model_params: ModelWithParameters):
        params = model_params.parameters
        x0 = model_params.get_estimation("X0")
        y0 = model_params.get_estimation("Y0")
        phi = model_params.get_estimation("Phi0")*np.pi/180
        u0 = model_params.get_estimation("U0")*PIXEL_SIZE
        nu = model_params.get_estimation("nu_")
        u_z = -nu

        k0 = params["k0"]
        k_start = params["k_start"]
        k_end = params["k_end"]


        x_start, y_start = calculate_position(x0,y0,u0,u_z,phi,k_start-k0)
        x_end, y_end = calculate_position(x0,y0,u0,u_z,phi,k_end-k0)

        return x_start,x_end,y_start,y_end

    def ask_nu(self, model_params: ModelWithParameters):
        return model_params.get_estimation("nu_")

    def ask_a(self, model_params: ModelWithParameters):
        return 0.0

    def ask_u0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("U0")

    def ask_phi0(self, model_params: ModelWithParameters):
        return model_params.get_estimation("Phi0")


LINEAR_TRACK_FORM_PSEUDO_A = LinearTrackModelPseudoAcceleration()
