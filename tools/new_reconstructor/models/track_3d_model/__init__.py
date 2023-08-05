import numpy as np
import pymc as pm
import pytensor.tensor as pt
from utils import binsearch_tgt
from vtl_common.common_GUI.tk_forms import TkDictForm

from ..model_base import ModelWithParameters, ReconstructionModelWrapper
from ..form_prototypes import DistributionField, PassthroughField
from ..linear_track_model.light_curves import create_lc_alter
from ..linear_track_model import x0_from_pmt, y0_from_pmt
from vtl_common.parameters import PIXEL_SIZE, HALF_GAP_SIZE, HALF_PIXELS, MAIN_LATITUDE,MAIN_LONGITUDE
from common_functions import create_coord_mesh, ensquared_energy_full, create_temporal_part

from fixed_rotator import Vector3, radec_to_ocef, unixtime_to_era, eci_to_ocef, ocef_to_detector_plane, detector_plane_to_ocef_f
from fixed_rotator import Vector2

SPAN = HALF_GAP_SIZE + HALF_PIXELS*PIXEL_SIZE

def get_interpolated_time(time_src:np.ndarray, float_index:float):
    start_index = int(float_index)
    assert start_index < time_src.shape[0]-1
    t1 = time_src[start_index]
    t2 = time_src[start_index+1]
    return t1 + (t2-t1)*(float_index-start_index)

class SpatialTrackModel(ReconstructionModelWrapper):

    # Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    # One more static field is already included in ReconstructionModelWrapper:
    # It gets transformed into pymc.<Distribution>. And it automatically adds error parameters (sigma, nu, etc...).

    # final_distribution = FinalDistributionField()

    accel = DistributionField("const", 0.0)
    V0 = DistributionField("normal", mu=0.05, sigma=5)
    dec_deg = DistributionField("const", 33.0)
    ra_deg = DistributionField("const", 112.0)
    LC = PassthroughField(create_lc_alter)
    distance0 = DistributionField("uniform", lower=70, upper=100)
    # x_pdm0 = DistributionField("uniform", lower=-SPAN, upper=SPAN)
    # y_pdm0 = DistributionField("uniform", lower=-SPAN, upper=SPAN)


    def generate_pymc_model(self, observed, cut_start, cut_end, broken, pmt, reconstructor_main) -> ModelWithParameters:
        print("CUTTER:", cut_start, cut_end)
        orientation_form: TkDictForm = reconstructor_main.orientation_form
        orientation = orientation_form.get_values()
        v_lat = orientation["VIEW_LATITUDE"] * np.pi / 180
        v_lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
        self_rotation = orientation["SELF_ROTATION"] * np.pi / 180
        main_lat = MAIN_LATITUDE * np.pi / 180
        main_lon = MAIN_LONGITUDE * np.pi / 180
        f = orientation["FOCAL_DISTANCE"]
        observed_signal, observed_time = observed
        with pm.Model() as model:
            consts = dict()
            k_start = cut_start
            k_end = cut_end
            T = k_end - k_start
            k0 = (k_start + k_end) / 2
            UT0 = get_interpolated_time(observed_time,k0)

            mesh_x, mesh_y, mesh_t = create_coord_mesh(T, k_start)
            alive = np.logical_not(broken)
            pixel_xs = mesh_x[:, alive]
            pixel_ys = mesh_y[:, alive]
            pixel_ks = mesh_t[:, alive]
            pixel_uts = create_temporal_part(T,k_start,observed_time)[:,alive]
            delta_k = pixel_ks - k0

            #pixel_era = unixtime_to_era(pixel_uts)
            dist = self.distance0("distance", consts)
            x0p = x0_from_pmt(pmt)
            y0p = y0_from_pmt(pmt)
            # x0p = self.x_pdm0("x0_pdm", consts)
            # y0p = self.y_pdm0("y0_pdm", consts)
            pdm0 = Vector2(x0p,y0p)
            point0 = detector_plane_to_ocef_f(pdm0, f)*dist

            x0,y0,z0 = point0.unpack()
            pm.Deterministic("z0_dev", x0)
            pm.Deterministic("x0_dev", y0)
            pm.Deterministic("y0_dev", z0)

            # x0 = self.x0("x0", consts)
            # y0 = self.y0("y0", consts)
            # z0 = self.z0("z0", consts)
            ra = self.ra_deg("RA", consts)*np.pi/180
            dec = self.dec_deg("DEC", consts)*np.pi/180
            v0 = self.V0("V0", consts)
            delta_t = (pixel_uts - UT0)



            # point_0_proj, vis = ocef_to_detector_plane(point0, f)
            # x0p, y0p = point_0_proj.unpack()
            # pm.Deterministic("x0_pdm", x0p)
            # pm.Deterministic("y0_pdm", y0p)

            vel_dev = -radec_to_ocef(ra,dec,v_lat, v_lon, self_rotation,unixtime_to_era(pixel_uts[0]))*v0
            dev_3d = point0 + vel_dev*delta_t

            eci2surf = eci_to_ocef(0.0, main_lat, main_lon, 0.0)
            eci2dev = eci_to_ocef(0.0, v_lat, v_lon, self_rotation)

            #surf2dev = eci2surf.inverse()*eci2dev

            dev2surf = eci2dev.inverse()*eci2surf
            surf_point0 = dev2surf*point0
            s0x, s0y,s0z = surf_point0.unpack()
            pm.Deterministic("z0_surf", s0x)
            pm.Deterministic("x0_surf", s0y)
            pm.Deterministic("y0_surf", s0z)
            #surf_3d = dev2surf*dev_3d

            # dev_3d = surf2dev*surf_3d
            # vel_dev = surf2dev*vel_surf

            fs_2d, visible = ocef_to_detector_plane(dev_3d, f)

            x_fs, y_fs = fs_2d.unpack()

            # XFS_VAR = pm.Deterministic("X_fs_mm", x_fs)
            # YFS_VAR = pm.Deterministic("Y_fs_mm", y_fs)


            sigmaPSF = pm.HalfNormal('SigmaPSF', 1.) * PIXEL_SIZE
            lc = self.LC.get_lc(delta_t, UT0, consts)
            intensity = lc * ensquared_energy_full(pixel_xs, pixel_ys, x_fs, y_fs, sigmaPSF)
            obs = observed_signal[k_start:k_end]
            observed_var = self.final_distribution('OBSERVED', mu=intensity,
                                                   observed=obs[:, alive])
        return ModelWithParameters(model, {
            "k0": k0,
            "k_start": k_start,
            "k_end": k_end,
            "v_lat":v_lat,
            "v_lon":v_lon,
            "f": f,
            "self_rotation":self_rotation,
            "UT0": UT0,
            "times": observed_time[k_start:k_end],
        }, consts)

    def reconstruction_overlay(self, plotter, model_params: ModelWithParameters):
        params = model_params.parameters
        # k0 = params["k0"]
        k_start = params["k_start"]
        k_end = params["k_end"]
        v_lat = params["v_lat"]
        v_lon = params["v_lon"]
        UT0 = params["UT0"]
        f = params["f"]
        self_rotation = params["self_rotation"]
        x0 = float(model_params.get_estimation("z0_dev"))
        y0 = float(model_params.get_estimation("x0_dev"))
        z0 = float(model_params.get_estimation("y0_dev"))
        ra = float(model_params.get_estimation("RA"))*np.pi/180
        dec = float(model_params.get_estimation("DEC"))*np.pi/180
        v0 = float(model_params.get_estimation("V0"))
        pixel_uts = params["times"]

        point_0 = Vector3(x0, y0, z0)
        vel_dev = -radec_to_ocef(ra, dec, v_lat, v_lon, self_rotation, unixtime_to_era(pixel_uts)) * v0
        dev_3d = point_0 + vel_dev * (pixel_uts - UT0)

        fs_2d, visible = ocef_to_detector_plane(dev_3d, f)

        X0 = float(model_params.get_estimation("X0"))
        Y0 = float(model_params.get_estimation("Y0"))

        k0 = model_params.parameters["k0"]
        plotter.set_origin(X0, Y0, k0)

        fs_x, fs_y = fs_2d.unpack()
        print("OVERLAY DIM",fs_x.shape[0])
        print("OVERLAY RANGE",k_start, k_end)
        r = model_params.get_estimation("SigmaPSF")*PIXEL_SIZE
        plotter.plot_lines(fs_x, fs_y, "-o", color="red", linewidth=r*2)

    def postprocess(self, ax, model_params: ModelWithParameters, actual_x):
        trace = model_params.idata
        k_start = model_params.parameters["k_start"]
        k_end = model_params.parameters["k_end"]
        t = model_params.parameters["times"]
        UT0 = model_params.parameters["UT0"]

        kk = np.arange(k_start, k_end)
        k0 = model_params.parameters["k0"]
        self.LC.postprocess_plot(t-UT0, UT0, ax, model_params, model_params.pmt, actual_x=actual_x[kk])

TRACK_3d_FORM = SpatialTrackModel()
