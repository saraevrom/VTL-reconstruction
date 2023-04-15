import numpy as np
import pymc as pm
import pytensor.tensor as pt

from orientation.stellar_math import unixtime_to_era
from reconstruction.common import d_erf
from .stellar_tensor_math import eci_to_ocef_pt, ocef_to_detector_plane_pt, rotate_yz_pt
from utils import binsearch_tgt
from .time_interval import TimeRange
from .database_reader import StarEntry
from vtl_common.parameters import PIXEL_SIZE
from vtl_common.localized_GUI.plotter import LOWER_EDGES
import numba as nb

PIXEL_POSITIONS = LOWER_EDGES+PIXEL_SIZE/2


@nb.njit(cache=True)
def create_coord_mesh(T):
    result_x = np.zeros(shape=(T, 16, 16))
    result_y = np.zeros(shape=(T, 16, 16))
    for k in range(T):
        for i in range(16):
            for j in range(16):
                result_x[k,i,j] = PIXEL_POSITIONS[i]
                result_y[k,i,j] = PIXEL_POSITIONS[j]
    return result_x, result_y

def ensquared_energy_full(x_mesh, y_mesh, x0, y0, psf):
    scale = np.sqrt(2)*psf
    a = d_erf(x0 - x_mesh, scale, pixel_size=PIXEL_SIZE) * d_erf(y0 - y_mesh, scale, pixel_size=PIXEL_SIZE)
    return a

def create_model(datafile, intervals, stars, known_params, unixtime, tuner, broken, ffmodel):
    era = unixtime_to_era(unixtime)
    ut0 = np.array(datafile["UT0"])
    times = []
    observed = []
    for interval in intervals:
        interval: TimeRange
        ut_start, ut_end = interval.unixtime_intervals()
        i_start = binsearch_tgt(ut0, ut_start)
        i_end = binsearch_tgt(ut0, ut_end)
        times.append(ut0[i_start:i_end:interval.stride])
        observed.append(datafile["data0"][i_start:i_end:interval.stride])
        print("INTERVAL SRC", interval.name())
        print(f"INTERVAL {i_start} - {i_end}")
    times = np.concatenate(times)
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

    observed[break_matrix] = 0.0

    T = times.shape[0]
    x_mesh, y_mesh = create_coord_mesh(T)
    with pm.Model() as model:
        tune_lat = tuner["tune_lat"]
        tune_lon = tuner["tune_lon"]
        tune_rot = tuner["tune_rot"]
        tune_f = tuner["tune_f"]
        tune_psf = tuner["tune_psf"]
        use_laplace = tuner["use_laplace"]

        view_latitude = tune_lat("lat", known_params["VIEW_LATITUDE"])
        view_longitude = tune_lon("lon", known_params["VIEW_LONGITUDE"])
        self_rotation = tune_rot("Ω", known_params["SELF_ROTATION"])
        focal_distance = tune_f("f", known_params["FOCAL_DISTANCE"])
        psf = tune_psf("PSF", known_params["PSF"])*PIXEL_SIZE

        # view_latitude = pm.TruncatedNormal("lat",
        #                                    lower=clamp_lat_min,
        #                                    upper = clamp_lat_max,
        #                                    mu=known_params["VIEW_LATITUDE"],
        #                                    sigma=tune_lat,
        #                                    )
        #
        # view_longitude = pm.TruncatedNormal("lon",
        #                                     lower=clamp_lon_min,
        #                                     upper=clamp_lon_max,
        #                                     mu=known_params["VIEW_LONGITUDE"],
        #                                     sigma=tune_lon,
        #                                     )

        #self_rotation = pm.Normal("Ω", mu=known_params["SELF_ROTATION"], sigma=tuner["tune_rot"])

        #focal_distance = pm.TruncatedNormal("f", lower=100, upper=200, mu=known_params["FOCAL_DISTANCE"], sigma=tuner["tune_f"])
        # if tuner["fix_f"]:
        #     focal_distance = tuner["tune_f"]
        # else:
        #     focal_distance = pm.Normal("f", mu=known_params["FOCAL_DISTANCE"], sigma=tuner["tune_f"])
        #psf = pm.TruncatedNormal("PSF", lower=0.01*PIXEL_SIZE, upper=10.0*PIXEL_SIZE, mu=known_params["PSF"], sigma=10.0)
        #psf = pm.Uniform("PSF", lower=0.01*PIXEL_SIZE, upper=10.0*PIXEL_SIZE)

        off_mu = np.median(datafile["data0"])
        off_sigma = np.mean((datafile["data0"] - off_mu) ** 2) ** 0.5

        max_energy = np.max([star.energy() for star in stars])
        max_light = np.max(datafile["data0"] - off_mu)

        mul_mu = max_light/max_energy

        print("ASSUMPTIONS:")
        print("OFF ~", off_mu)
        print("MUL ~", mul_mu)

        mul = pm.TruncatedNormal("A", mu=mul_mu, sigma=100.0, lower=0.0)
        off = pm.Normal("B", sigma=off_sigma, mu=off_mu)


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
            x_eci, y_eci, z_eci = star.get_eci()
            energy = star.energy()
            x_ocef, y_ocef, z_ocef = eci_to_ocef_pt(x_eci, y_eci, z_eci, era,
                                                    lat=view_latitude * np.pi / 180,
                                                    lon=view_longitude * np.pi / 180)
            x_ocef, y_ocef, z_ocef = rotate_yz_pt(x_ocef, y_ocef, z_ocef, self_rotation * np.pi / 180)
            x_pdf, y_pdf, v = ocef_to_detector_plane_pt(x_ocef, y_ocef, z_ocef, focal_distance)
            energy_dist = energy*v*ensquared_energy_full(x_mesh, y_mesh, x_pdf, y_pdf, psf)
            if modelled is None:
                modelled = energy_dist
                print("Started with", star.primary_name())
            else:
                modelled = modelled + energy_dist
                print("Added", star.primary_name())

        nonbroken_observed = mul*modelled+off
        broken_observed = pt.switch(break_matrix, 0.0, nonbroken_observed)  # Turning off broken pixels

        if use_laplace:
            b = pm.HalfNormal("lapl_b", sigma=1.0)
            intensity = pm.Laplace("I", mu=broken_observed, b=b, observed=observed)
        else:
            # Normal distributed variable may save some time
            sigma = pm.HalfNormal("σ", sigma=1.0)
            intensity = pm.Normal("I", mu=broken_observed, sigma=sigma, observed=observed)

    return model