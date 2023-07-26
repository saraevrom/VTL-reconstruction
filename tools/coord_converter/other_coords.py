import datetime
import tkinter as tk

import numpy as np

from vtl_common.common_GUI.settings_frame import SettingMenu, DoubleValue
from vtl_common.localization import get_locale
from tools.orientation.orientation.stellar_math import eci_to_ocef, ocef_to_detector_plane
from tools.orientation.orientation.stellar_math import ocef_to_altaz, radec_to_eci, datetime_to_era
from tools.orientation.orientation.stellar_math import ocef_to_eci, detector_plane_to_ocef_f
from tools.orientation.orientation.stellar_math import altaz_to_ocef, eci_to_radec, rotate_yz


class CoordSystem(object):
    def __init__(self, settings_frame: SettingMenu, schedule):
        self.settings_frame = settings_frame
        self.schedule = schedule
        self.build_menu()

    def build_menu(self):
        pass

    def get_eci(self, location:dict, orientation:dict, dt:datetime.datetime):
        raise NotImplementedError()

    def set_eci(self, eci, location:dict, orientation:dict, dt:datetime.datetime):
        raise NotImplementedError()

    def part_changed(self):
        self.schedule.append(self)

    def add_tracer(self, key):
        self.settings_frame.add_tracer(key, self.part_changed)


class AstronomicalSystem(CoordSystem):

    def build_menu(self):
        self.settings_frame.add_separator(get_locale("converter.astronomy"))
        self.ra_hours = self.settings_frame.add_setting(DoubleValue, "ra_h", "RA [h]", 0.0)
        self.ra_deg = self.settings_frame.add_setting(DoubleValue, "ra_deg", "RA [°]", 0.0)
        self.dec_deg = self.settings_frame.add_setting(DoubleValue, "dec_deg", "Dec [°]", 0.0)
        self.settings_frame.add_tracer("ra_h", lambda: self.ra_h_to_deg(True))
        self.settings_frame.add_tracer("ra_deg", lambda: self.ra_deg_to_h(True))

        self.add_tracer("dec_deg")

    def ra_h_to_deg(self, auto=False):
        d = self.ra_hours.get_value()*180/12
        d = round(d, 6)
        self.ra_deg.set_value(d)
        if auto:
            self.part_changed()

    def ra_deg_to_h(self, auto=False):
        h = self.ra_deg.get_value()*12/180
        h = round(h, 6)
        self.ra_hours.set_value(h)
        if auto:
            self.part_changed()

    def get_eci(self, location, orientation, dt):
        ra = self.ra_deg.get_value()*np.pi/180
        dec = self.dec_deg.get_value()*np.pi/180
        x,y,z = radec_to_eci(ra, dec)
        return x,y,z

    def set_eci(self, eci, location, orientation, dt):
        x, y, z = eci
        ra, dec = eci_to_radec(x, y, z)
        dec_deg = dec*180/np.pi
        ra_deg = ra*180/np.pi
        dec_deg = round(dec_deg,3)
        ra_deg = round(ra_deg,3)

        self.dec_deg.set_value(dec_deg)
        self.ra_deg.set_value(ra_deg)
        self.ra_deg_to_h()

class GeoSystem(CoordSystem):
    def build_menu(self):
        self.settings_frame.add_separator(get_locale("converter.altaz"))
        self.zang_deg = self.settings_frame.add_setting(DoubleValue, "alt_deg", "θ [°]", 0.0)
        self.az_deg = self.settings_frame.add_setting(DoubleValue, "az_deg", "φ [°]", 0.0)
        self.add_tracer("alt_deg")
        self.add_tracer("az_deg")

    def get_eci(self, location, orientation, dt):
        lat = location["lat"]*np.pi/180
        lon = location["lon"]*np.pi/180
        alt = (90-self.zang_deg.get_value())*np.pi/180
        az = self.az_deg.get_value()*np.pi/180

        x0, y0, z0 = altaz_to_ocef(alt, az)
        era = datetime_to_era(dt)
        return ocef_to_eci(x0,y0,z0,era,lat,lon)

    def set_eci(self, eci, location:dict, orientation:dict, dt:datetime.datetime):
        lat = location["lat"] * np.pi / 180
        lon = location["lon"] * np.pi / 180
        era = datetime_to_era(dt)
        x_eci, y_eci, z_eci = eci
        x0, y0, z0 = eci_to_ocef(x_eci, y_eci, z_eci, era, lat, lon)
        alt, az = ocef_to_altaz(x0, y0, z0)

        az = az*180/np.pi
        az = round(az, 6)
        alt = alt*180/np.pi
        zang = 90-alt
        zang = round(zang, 6)
        self.az_deg.set_value(az)
        self.zang_deg.set_value(zang)


class DeviceSystem(CoordSystem):
    def build_menu(self):
        self.settings_frame.add_separator(get_locale("converter.device"))
        self.gamma_deg = self.settings_frame.add_setting(DoubleValue, "ang_sep", "γ [°]", 0.0)
        self.psi_deg = self.settings_frame.add_setting(DoubleValue, "psi", "ψ [°]", 0.0)

        self.add_tracer("psi_deg")
        self.add_tracer("gamma_deg")

    def get_eci(self, location, orientation, dt):
        lat = orientation["VIEW_LATITUDE"]*np.pi/180
        lon = orientation["VIEW_LONGITUDE"]*np.pi/180
        self_rot = orientation["SELF_ROTATION"]*np.pi/180
        gamma = (90-self.gamma_deg.get_value())*np.pi/180
        psi = (self.psi_deg.get_value() - 90)*np.pi/180

        x0, y0, z0 = altaz_to_ocef(gamma, psi)
        x0, y0, z0 = rotate_yz(x0, y0, z0, -self_rot)
        era = datetime_to_era(dt)
        return ocef_to_eci(x0,y0,z0,era,lat,lon)

    def set_eci(self, eci, location:dict, orientation:dict, dt:datetime.datetime):
        lat = orientation["VIEW_LATITUDE"] * np.pi / 180
        lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
        self_rot = orientation["SELF_ROTATION"] * np.pi / 180
        era = datetime_to_era(dt)
        x_eci, y_eci, z_eci = eci
        x0, y0, z0 = eci_to_ocef(x_eci, y_eci, z_eci, era, lat, lon)
        x0, y0, z0 = rotate_yz(x0, y0, z0, self_rot)
        gamma, psi = ocef_to_altaz(x0, y0, z0)

        gamma = gamma * 180 / np.pi
        gamma = 90-gamma
        gamma = round(gamma, 6)

        psi = (psi*180/np.pi+90) % 360
        if psi>180:
            psi = psi-360

        psi = round(psi, 6)

        self.gamma_deg.set_value(gamma)
        self.psi_deg.set_value(psi)

class PlanarSystem(CoordSystem):
    def build_menu(self):
        self.settings_frame.add_separator(get_locale("converter.planar"))
        self.x = self.settings_frame.add_setting(DoubleValue, "x_mm", "X [mm]", 0.0)
        self.y = self.settings_frame.add_setting(DoubleValue, "y_mm", "Y [mm]", 0.0)

        self.add_tracer("x_mm")
        self.add_tracer("y_mm")

    def get_eci(self, location, orientation, dt):
        lat = orientation["VIEW_LATITUDE"]*np.pi/180
        lon = orientation["VIEW_LONGITUDE"]*np.pi/180
        self_rot = orientation["SELF_ROTATION"]*np.pi/180
        f = orientation["FOCAL_DISTANCE"]
        x = self.x.get_value()
        y = self.y.get_value()

        x0, y0, z0 = detector_plane_to_ocef_f(x, y, f)
        x0, y0, z0 = rotate_yz(x0, y0, z0, -self_rot)
        era = datetime_to_era(dt)
        return ocef_to_eci(x0,y0,z0,era,lat,lon)

    def set_eci(self, eci, location:dict, orientation:dict, dt:datetime.datetime):
        lat = orientation["VIEW_LATITUDE"] * np.pi / 180
        lon = orientation["VIEW_LONGITUDE"] * np.pi / 180
        self_rot = orientation["SELF_ROTATION"] * np.pi / 180
        f = orientation["FOCAL_DISTANCE"]
        era = datetime_to_era(dt)
        x_eci, y_eci, z_eci = eci
        x0, y0, z0 = eci_to_ocef(x_eci, y_eci, z_eci, era, lat, lon)
        x0, y0, z0 = rotate_yz(x0, y0, z0, self_rot)
        #gamma, psi = ocef_to_altaz(x0, y0, z0)
        x, y, v = ocef_to_detector_plane(x0, y0, z0, f)
        if not v:
            x = 0
            y = 0

        x = round(x, 2)
        y = round(y, 2)
        self.x.set_value(x)
        self.y.set_value(y)

class CoordList(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._formdata = dict()
        self._schedule = []
        self.settings_frame = SettingMenu(self, autocommit=True)
        self.settings_frame.commit_action = self.on_form_commit
        self.settings_frame.pack(fill="both", expand=True)
        self.systems = [
            AstronomicalSystem(self.settings_frame,self._schedule),
            GeoSystem(self.settings_frame,self._schedule),
            DeviceSystem(self.settings_frame,self._schedule),
            PlanarSystem(self.settings_frame,self._schedule),
        ]
        self._orientation = None
        self._location = None
        self._time = None

    def set_fixed(self, orientation, location, time):
        self._orientation = orientation
        self._location = location
        self._time = time
        self.syncto(self.systems[0])

    def syncto(self, src):
        self.settings_frame.notify_disable()
        if self._time is not None:
            for dst in self.systems:
                if dst != src:
                    print(src, ">>>", dst)
                    eci = src.get_eci(self._location, self._orientation, self._time)
                    dst.set_eci(eci, self._location, self._orientation, self._time)
        self.settings_frame.notify_enable()

    def on_form_commit(self):
        self.settings_frame.push_settings_dict(self._formdata)
        while self._schedule:
            src = self._schedule.pop(0)
            self.syncto(src)
