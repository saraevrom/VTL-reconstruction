import numpy as np
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE
from fixed_rotator.slowmath import Quaternion
from fixed_rotator.astro_math_z_aligned import latlon_quaternion


def hor_to_dev(orientation):
    self_rot = orientation["SELF_ROTATION"] * np.pi / 180
    dev_dec = orientation["VIEW_LATITUDE"] * np.pi / 180
    dev_gha = orientation["VIEW_LONGITUDE"] * np.pi / 180
    dev_lat = MAIN_LATITUDE * np.pi / 180
    dev_lon = MAIN_LONGITUDE * np.pi / 180
    # R matrix from article
    R = Quaternion.rotate_xy(self_rot) * \
        latlon_quaternion(dev_dec, dev_gha).inverse() * \
        latlon_quaternion(dev_lat, dev_lon)
    #assert R*R.inverse() == Quaternion(w=1,x=0,y=0,z=0)
    return R