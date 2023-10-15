import numpy as np

from vtl_common.common_GUI.tk_forms_assist import FormNode, AlternatingNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from tools.new_reconstructor.orientation_pointing_form import OrientedPoint, Declination
from fixed_rotator.astro_math_z_aligned import Vector2, Vector3, Quaternion, latlon_quaternion
from fixed_rotator.astro_math_z_aligned import radec_to_eci, eci_to_ocef
from vtl_common.localization import get_locale
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE


class DirectionDef(object):
    def calculate(self,orientation,era):
        raise NotImplementedError

class RaDecDef(DirectionDef):
    def __init__(self,ra,dec):
        self.ra = ra
        self.dec = dec

    def calculate(self,orientation,era):
        ra = self.ra
        dec = self.dec
        self_rot = orientation["SELF_ROTATION"] * np.pi / 180
        dev_dec = orientation["VIEW_LATITUDE"] * np.pi / 180
        dev_gha = orientation["VIEW_LONGITUDE"] * np.pi / 180
        # dev_lat = MAIN_LATITUDE * np.pi / 180
        # dev_lon = MAIN_LONGITUDE * np.pi / 180
        v_eci = - radec_to_eci(ra, dec)

        # Velocity in detector frame
        v_dev = eci_to_ocef(era, dev_dec, dev_gha, self_rot) * v_eci
        return v_dev

class HorizontalDef(DirectionDef):
    def __init__(self,alt,az):
        self.alt = alt
        self.az = az

    def calculate(self,orientation,era):
        self_rot = orientation["SELF_ROTATION"] * np.pi / 180
        dev_dec = orientation["VIEW_LATITUDE"] * np.pi / 180
        dev_gha = orientation["VIEW_LONGITUDE"] * np.pi / 180
        dev_lat = MAIN_LATITUDE * np.pi / 180
        dev_lon = MAIN_LONGITUDE * np.pi / 180

        R_quat = Quaternion.rotate_xy(self_rot) * \
                 latlon_quaternion(dev_dec, dev_gha).inverse() * \
                 latlon_quaternion(dev_lat, dev_lon)

        z = np.sin(self.alt)
        y = np.cos(self.alt)*np.cos(self.az)
        x = np.cos(self.alt)*np.sin(self.az)
        v_hor = Vector3(x,y,z)
        v_dev = R_quat*v_hor
        return v_dev


class RaDecEntry(OrientedPoint):
    DISPLAY_NAME = get_locale("tools.trajectory_calculator.radiant")
    def get_data(self):
        data = super().get_data()
        ra = data["ra"]
        dec = data["dec"]
        return RaDecDef(ra=ra,dec=dec)


class HorEntry(FormNode):
    DISPLAY_NAME = get_locale("tools.trajectory_calculator.horizontal")
    FIELD__alt = create_value_field(Declination, "alt [°]", 0.0)
    FIELD__az = create_value_field(Declination, "az [°]", 0.0)

    def get_data(self):
        data = super().get_data()
        alt = data["alt"]
        az = data["az"]
        return HorizontalDef(alt=alt,az=az)


class ExtendedDirection(AlternatingNode):
    DISPLAY_NAME = get_locale("tools.trajectory_calculator.direction")
    SEL__sky = RaDecEntry
    SEL__horizon = HorEntry
