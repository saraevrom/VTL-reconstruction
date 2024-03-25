import numpy as np
import pymc as pm
import pytensor.tensor as pt
import h5py
import numba as nb

# from .stellar_tensor_math import eci_to_ocef_pt, ocef_to_detector_plane_pt, rotate_yz_pt, ecef_to_ocef_pt
# from .stellar_tensor_math import ocef_to_altaz_pt
from fixed_rotator.astro_math_z_aligned import eci_to_ocef, ocef_to_detector_plane, Quaternion, ecef_to_ocef, Vector3, altaz_to_ocef
from fixed_rotator.astro_math_z_aligned import ocef_to_altaz, unixtime_to_era

from utils import binsearch_tgt
from .time_interval import TimeRange
from .database_reader import StarEntry
from vtl_common.parameters import PIXEL_SIZE
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE

from common_functions import create_coord_mesh
from common_functions import ensquared_energy_full

def get_dual(h5file:h5py.File, primary, secondary):
    if primary in h5file.keys():
        return h5file[primary]
    else:
        return h5file[secondary]

def get_signal(h5file):
    return get_dual(h5file,"data0","pdm_2d_rot_global")

def get_time(h5file):
    return get_dual(h5file,"UT0","unixtime_dbl_global")


@nb.njit()
def mean_axis0(data):
    res = np.zeros(shape=data.shape[1:])
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            res[i,j] = np.mean(data[:,i,j])
    return res

@nb.njit()
def collapse(data, stride):
    src_len = data.shape[0]
    tgt_len = src_len//stride
    #print("COLLAPSE", src_len, tgt_len)
    res = np.zeros(shape=(tgt_len,data.shape[1],data.shape[2]))
    for i in range(tgt_len):
        res[i] = mean_axis0(data[i*stride:(i+1)*stride])
    return res

@nb.njit()
def get_times(time_src,stride):
    src_len = time_src.shape[0]
    tgt_len = src_len//stride
    res = np.zeros(shape=(tgt_len,))
    print("TIMES", src_len,tgt_len)
    for i in range(tgt_len):
        res[i] = np.mean(time_src[i*stride:(i+1)*stride])
    return res


def mux_error(uniform,constructor,*args,**kwargs):
    if uniform:
        return constructor(*args, **kwargs)
    else:
        return constructor(*args, shape=(16, 16), **kwargs)


ERROR_DICT = {
    "laplace":"lapl_b",
    "student":"Sigma0",
    "gauss":"σ"
}



def create_model(datafile, intervals, stars, known_params, unixtime, tuner, broken, ffmodel):
    era = unixtime_to_era(unixtime)
    ut0 = np.array(get_time(datafile))
    times = []
    observed = []
    for interval in intervals:
        interval: TimeRange
        ut_start, ut_end = interval.unixtime_intervals()
        i_start = binsearch_tgt(ut0, ut_start)
        i_end = binsearch_tgt(ut0, ut_end)
        #times.append(ut0[i_start:i_end:interval.stride])
        new_time = get_times(ut0[i_start:i_end], stride=interval.stride)
        #new_time = list(map())
        times.append(new_time)
        #observed.append(get_signal(datafile)[i_start:i_end:interval.stride])
        observed.append(collapse(get_signal(datafile)[i_start:i_end],stride=interval.stride))
        print("INTERVAL SRC", interval.name())
        print(f"INTERVAL {i_start} - {i_end}")
    times = np.concatenate(times)
    eras = unixtime_to_era(times)
    observed = np.concatenate(observed)
    if ffmodel is not None:
        observed = ffmodel.apply(observed)
    assert times.shape[0] == observed.shape[0]
    assert observed.shape[1] == 16
    assert observed.shape[2] == 16

    break_matrix = np.expand_dims(broken, 0)
    print("BROKEN SHAPE (PRE)", break_matrix.shape)
    break_matrix = np.repeat(break_matrix, observed.shape[0], axis=0)
    print("OBSERVED SHAPE", observed.shape)
    print("BROKEN SHAPE", break_matrix.shape)
    assert break_matrix.shape == observed.shape
    print("USED PIXELS:", np.sum(broken))
    #
    # observed[break_matrix] = 0.0

    T = times.shape[0]
    x_mesh, y_mesh, _ = create_coord_mesh(T)
    with pm.Model() as model:
        tune_lat = tuner["tune_lat"]
        tune_lon = tuner["tune_lon"]
        tune_rot = tuner["tune_rot"]
        tune_f = tuner["tune_f"]
        tune_psf = tuner["tune_psf"]
        #use_laplace = tuner["use_laplace"]
        final_dist = tuner["final_dist"]
        tune_kappa = tuner["tune_kappa"]
        tune_b = tuner["tune_b"]
        tune_b_auto_assume = tuner["tune_b_auto_assume"]

        view_latitude = tune_lat("lat", known_params["VIEW_LATITUDE"]) #DEC
        view_longitude = tune_lon("lon", known_params["VIEW_LONGITUDE"]) #HOUR ANGLE
        self_rotation = tune_rot("Ω", known_params["SELF_ROTATION"])
        focal_distance = tune_f("f", known_params["FOCAL_DISTANCE"])
        psf = tune_psf("PSF", known_params["PSF"])*PIXEL_SIZE

        alive_matrix = np.logical_not(break_matrix)
        off_mu = np.median(observed[alive_matrix])
        off_sigma = np.mean((observed[alive_matrix] - off_mu) ** 2) ** 0.5

        max_energy = np.max([star.energy() for star in stars])
        max_value = np.max(observed[alive_matrix])
        max_light = max_value - off_mu

        #integral_mul = erf(1/(2**1.5*known_params["PSF"]))**2
        #assumed_energy = max_light/integral_mul

        mul_mu = max_light/max_energy

        print("ASSUMPTIONS:")
        #print("MAX_VALUE", max_value)
        #print("MAX_AMPL", max_light)
        #print("PSF CORRECTION", integral_mul)
        #print("MAX_ASSUMED_LIGHT", assumed_energy)
        print("BASE", off_mu)
        print("MUL ~", max_light, "/", max_energy, "=", mul_mu)

        #mul = pm.TruncatedNormal("A", mu=mul_mu, sigma=tune_a_std, lower=0.0)
        #mul = pm.HalfNormal("A", sigma=mul_mu*np.sqrt(np.pi/2))
        mul = tune_kappa("κ", known_params["MULTIPLIER"])
        if tune_b_auto_assume:
            print("OFFSET ~", off_mu, "±", off_sigma)
            off = pm.Normal("B", sigma=off_sigma, mu=off_mu)
        else:
            print("OFFSET ~", known_params["OFFSET"])
            off = tune_b("B", known_params["OFFSET"])


        # ["VIEW_LATITUDE", FloatNode, MAIN_LATITUDE],
        # ["VIEW_LONGITUDE", FloatNode, MAIN_LONGITUDE],
        # ["SELF_ROTATION", FloatNode, 0.0],
        # ["FOCAL_DISTANCE", FloatNode, 165.0],
        # ["PSF", FloatNode, 0.25],
        # ["MULTIPLIER", FloatNode, 1.0],
        # ["OFFSET", FloatNode, 0.0],

        modelled = None

        for star in stars:
            star:StarEntry
            #x_eci, y_eci, z_eci = star.get_eci()
            eci = Vector3(*star.get_eci())
            energy = star.energy_u()
            eci2ocef = eci_to_ocef(eras, dec=view_latitude * np.pi / 180, gha=view_longitude * np.pi / 180,
                                   self_rot=self_rotation * np.pi / 180, backend=pt)
            #eci2ocef = Quaternion.rotate_yz(self_rotation * np.pi / 180, backend=pt)*eci2ocef
            pdf_vec, v = ocef_to_detector_plane(eci2ocef*eci, focal_distance)
            x_pdf = pt.expand_dims(pdf_vec.x, (1, 2))
            y_pdf = pt.expand_dims(pdf_vec.y, (1, 2))
            v = pt.expand_dims(v, (1, 2))

            energy_dist = energy*v*ensquared_energy_full(x_mesh, y_mesh, x_pdf, y_pdf, psf)
            if modelled is None:
                modelled = energy_dist
                print("Started with", star.primary_name())
            else:
                modelled = modelled + energy_dist
                print("Added", star.primary_name())

        nonbroken_model = mul*modelled+off
        #broken_observed = pt.switch(break_matrix, 0.0, nonbroken_observed)  # Turning off broken pixels

        # Deterministic auxiliary parameters
        view_lat_rad = view_latitude * np.pi / 180
        view_lon_rad = view_longitude * np.pi / 180

        x_view_ecef = pt.cos(view_lon_rad)*pt.cos(view_lat_rad)
        y_view_ecef = pt.sin(view_lon_rad)*pt.cos(view_lat_rad)
        z_view_ecef = pt.sin(view_lat_rad)
        view_ecef = Vector3(x_view_ecef, y_view_ecef, z_view_ecef)

        ecef2ocef = ecef_to_ocef(MAIN_LATITUDE * np.pi / 180, MAIN_LONGITUDE * np.pi / 180)
        view_alt, view_az = ocef_to_altaz(ecef2ocef*view_ecef, backend=pt)
        view_azimuth = pm.Deterministic("AZ", view_az*180/np.pi)
        view_altangle = pm.Deterministic("ALT", view_alt*180/np.pi)
        alive = np.logical_not(broken)
        broken_model = nonbroken_model[:,alive]
        broken_observed = observed[:,alive]

        if final_dist=="laplace":
            b = pm.HalfNormal("lapl_b", sigma=1.0)
            intensity = pm.Laplace("I", mu=broken_model, b=b, observed=broken_observed)
        elif final_dist=="student":
            sigma0 = pm.HalfNormal('Sigma0', 1.0)
            nu = pm.Exponential('nu', 1.0  )
            intensity = pm.StudentT("I", mu=broken_model, sigma=sigma0, nu=nu, observed=broken_observed)
        else:
            # Normal distributed variable may save some time
            sigma = pm.HalfNormal("σ", sigma=1.0)
            intensity = pm.Normal("I", mu=broken_model, sigma=sigma, observed=broken_observed)

    return model
