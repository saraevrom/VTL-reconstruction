import datetime
import json
import tkinter as tk
import tkinter.messagebox as messagebox

import numpy as np

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.workspace_manager import Workspace
from vtl_common.datetime_parser import parse_datetimes_dt
from ..tool_base import ToolBase
from .input_form import make_form
from tools.orientation.parameters import ParametersForm
from specific_ui import DataOutput
from vtl_common.parameters import PIXEL_SIZE, MAIN_LATITUDE, MAIN_LONGITUDE
from fixed_rotator.astro_math_z_aligned import Vector2, Vector3, Quaternion, latlon_quaternion
from fixed_rotator.astro_math_z_aligned import radec_to_eci, eci_to_ocef, ocef_to_altaz
from fixed_rotator import datetime_to_era
from vtl_common.localization import get_locale
from vtl_common.common_GUI.button_panel import ButtonPanel
from .input_form import M_S, M_MRAD, M_MM, M_KM, DT_DEFAULT
from .mode_asker import OptionDialog

W_DETECTOR = get_locale("tools.trajectory_calculator.DETECTOR")
W_HORIZON = get_locale("tools.trajectory_calculator.HORIZON")

ORIENTATION_WORKSPACE = Workspace("orientation")
SKY_CATALOG = Workspace("sky_catalog")

def calculate_omega_proj(tres,F,X,Y, u_x,u_y):
    S = Vector3(-X, Y, 0)
    ez = Vector3.back()
    V = Vector3(-u_x, u_y, 0)
    omega = F ** 2 / (F ** 2 + S.sqr_len()) * (1 / F ** 2 * S.cross(V) + 1 / F * ez.cross(V)) / tres
    return omega

def calculate_dev_vector(F, x, y, omega, v_dev):
    direction = Vector3(-x,y,F).normalized()
    length = (direction.cross(v_dev)).length()/omega.length()
    return direction*length
    # S_0_3d = Vector3(-x, y, F)
    # omega_sqr = omega.sqr_len()
    # z_dev = F / (omega_sqr * S_0_3d.sqr_len()) * omega.dot(S_0_3d.cross(v_dev))
    # r_dev = z_dev / F * S_0_3d
    # return r_dev

def hor_to_dev(orientation):
    self_rot = orientation["SELF_ROTATION"] * np.pi / 180
    dev_dec = orientation["VIEW_LATITUDE"] * np.pi / 180
    dev_gha = orientation["VIEW_LONGITUDE"] * np.pi / 180
    dev_lat = MAIN_LATITUDE * np.pi / 180
    dev_lon = MAIN_LONGITUDE * np.pi / 180
    # R matrix from article
    return Quaternion.rotate_xy(self_rot) * \
             latlon_quaternion(dev_dec, dev_gha).conj() * \
             latlon_quaternion(dev_lat, dev_lon)

def calculate_dev_vector_bypass(U,R,F,v_dev):
    print("BYPASS",U,R,F,v_dev)
    z_dev_1 = (-F*v_dev.x-v_dev.z*R.x)/U.x
    z_dev_2 = (F*v_dev.y - v_dev.z*R.y)/U.y
    S_0_3d = Vector3(-R.x, R.y, F)
    return z_dev_1/F*S_0_3d, z_dev_2/F*S_0_3d

class LinearEstimator(ToolBase):
    TOOL_KEY = "tools.trajectory_calculator"

    def __init__(self, master):
        super().__init__(master)
        self._raw_parameters = None
        rpanel = tk.Frame(self)
        rpanel.pack(side="right",fill="both", expand=True)

        bpanel = ButtonPanel(rpanel)
        bpanel.add_button(get_locale("tools.trajectory_calculator.import_reco"), self.on_import_reco, 0)
        bpanel.add_button(get_locale("tools.trajectory_calculator.align_az"), self.on_balance_azimuth, 1)
        bpanel.add_button(get_locale("tools.trajectory_calculator.align_hor"), self.on_horizontal_align, 1)
        bpanel.add_button(get_locale("form.point.load"), self.on_point_load, 2)
        bpanel.add_button(get_locale("form.point.save"), self.on_point_save, 2)
        bpanel.pack(side="top", fill="x")

        self.orientation_form = SaveableTkDictForm(rpanel,
                                                   ParametersForm().get_configuration_root(),
                                                   save_label="form.orientation.save",
                                                   load_label="form.orientation.load",
                                                   file_asker=ORIENTATION_WORKSPACE,
                                                   optional=True,
                                                   option_default_visible=False
                                                   )
        self.orientation_form.pack(side="bottom",fill="x")

        self.parameter_form, self.parameter_parser = make_form(rpanel)
        self.parameter_form.pack(side="bottom",fill="both",expand=True)
        self._parameters = None
        self._orientation = None
        self.parameter_form.on_commit = self.on_commit
        self.orientation_form.on_commit = self.on_commit

        self.output_panel = DataOutput(self)
        self.output_panel.pack(side="left", fill="both", expand=True)
        self.on_commit()

    def on_point_save(self):
        #radiant = self._raw_parameters["radiant"]
        direction = self._raw_parameters["direction"]
        if not direction["selection_type"]=="sky":
            return

        radiant = direction["value"]

        filename = SKY_CATALOG.asksaveasfilename(title=get_locale("app.filedialog.save_settings.title"),
                                                     filetypes=[
                                                         (get_locale("app.filedialog_formats.form_json"), "*.json")],
                                                     initialdir=".",
                                                     parent=self)
        if filename:
            with open(filename, "w") as fp:
                json.dump(radiant, fp, indent=4)

    def on_point_load(self):
        filename = SKY_CATALOG.askopenfilename(title=get_locale("app.filedialog.open_settings.title"),
                                                   filetypes=[
                                                       (get_locale("app.filedialog_formats.form_json"), "*.json")],
                                                   initialdir=".",
                                                   parent=self)
        if filename:
            with open(filename, "r") as fp:
                jsd = json.load(fp)

            r = {
                "direction": {
                    "selection_type":"sky",
                    "value":jsd,
                }
            }
            self.parameter_form.set_values(r, force=True)
            self.parameter_form.trigger_change()

    def on_import_reco(self):
        reco, dt = self.ask_from_tool(0,"reco_result")
        if reco is None or len(reco.keys()) == 0:
           return

        reco:dict
        if len(reco.keys()) == 1:
            selected_mode = list(reco.keys())[0]
        else:
            selected_mode = OptionDialog(self, list(reco.keys())).result
        if selected_mode is None:
            return

        applier = dict()
        applier.update(reco[selected_mode])

        if dt is not None:
            applier["k0_datetime"] = dt.strftime(DT_DEFAULT)
        self.parameter_form.set_values(applier)
        self.on_commit()

    def on_commit(self):
        orientation = self.orientation_form.get_values()
        parameters = self.parameter_form.get_values()
        self._raw_parameters = parameters
        self.parameter_parser.parse_formdata(parameters)
        self._parameters = self.parameter_parser.get_data()
        self._orientation = orientation

        u0 = self._parameters["u0"]*PIXEL_SIZE
        phi0 = self._parameters["phi0"]*np.pi/180
        x0 = self._parameters["x0"]
        y0 = self._parameters["y0"]
        k0 = self._parameters["k0"]
        k = self._parameters["k"]
        tres = self._parameters["tres"] # Temporal resolution
        v = self._parameters["v"]
        F = self._orientation["FOCAL_DISTANCE"]

        base_dt = datetime.datetime.now()
        base_dt = base_dt.replace(microsecond=0)
        dt = parse_datetimes_dt(self._parameters["k0_datetime"], base_dt)
        era = datetime_to_era(dt+datetime.timedelta(seconds=(k-k0)*tres))

        X,Y,u_x, u_y = self.get_kinematics()

        self.output_panel.clear()
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.focal_plane"))
        self.output_panel.add_entry("x",f"{X:.3f}")
        self.output_panel.add_entry("y",f"{Y:.3f}")
        self.output_panel.add_entry("U",f"{((u_x**2+u_y**2)**0.5)/PIXEL_SIZE:.3f}")

        # omega = calculate_omega(tres, F, X, Y, u, phi0)
        omega = calculate_omega_proj(tres, F, X, Y, u_x, u_y)

        # R matrix from article
        R_quat = hor_to_dev(self._orientation)
        Rt_quat = R_quat.conj()

        omega_hor = Rt_quat*omega

        # self.output_panel.add_entry("omega vector (DETECTOR VIEW) [mrad/t]", f"[{wx:.3f}, {wy:.3f}, {wz:.3f}]")
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.omega"))
        self.output_panel.add_entry(f"ω ({W_DETECTOR}), [{M_MRAD}/{M_S}]", str(1000*omega))
        self.output_panel.add_entry(f"ω ({W_HORIZON}), [{M_MRAD}/{M_S}]", str(1000*omega_hor))
        self.output_panel.add_entry(f"ω [{M_MRAD}/{M_S}]", f"{1000*omega.length():.3f}")

        v_dev = v*self._parameters["direction"].calculate(orientation,era)
        v_hor = Rt_quat*v_dev

        u_k0_x = u0*np.cos(phi0)
        u_k0_y = u0*np.sin(phi0)

        omega_0 = calculate_omega_proj(tres, F, x0, y0, u_k0_x, u_k0_y)

        r_dev = calculate_dev_vector(F, x0, y0, omega_0, v_dev)
        r_dev = r_dev + v_dev*(k-k0)*tres
        r_hor = Rt_quat*r_dev
        h = r_hor.z
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.distances"))
        self.output_panel.add_entry(f"v_dev [{M_KM}/{M_S}]", str(v_dev))
        self.output_panel.add_entry(f"v_hor [{M_KM}/{M_S}]", str(v_hor))
        self.output_panel.add_entry(f"dev [{M_KM}]", str(r_dev))
        self.output_panel.add_entry(f"H [{M_KM}]", h)

        # ALT

        U = Vector2(u_x, u_y)/tres
        R = Vector2(X,Y)
        dev_1, dev_2 = calculate_dev_vector_bypass(U,R,F,v_dev)
        r_hor_1 = Rt_quat * dev_1
        r_hor_2 = Rt_quat * dev_2
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.distances_alt"))
        self.output_panel.add_entry("U [mm/s]", str(U))
        self.output_panel.add_entry("R [mm]", str(R))
        self.output_panel.add_entry("dev1 [km]", str(dev_1))
        self.output_panel.add_entry("dev2 [km]", str(dev_2))
        h1 = r_hor_1.z
        h2 = r_hor_2.z
        self.output_panel.add_entry("h [km]", f'{h1:.3f}; {h2:.3f}')

    def pseudo_acc_solve_dev(self):
        tres = self._parameters["tres"]  # Temporal resolution
        F = self._orientation["FOCAL_DISTANCE"]
        u0 = self._parameters["u0"] * PIXEL_SIZE / tres
        phi0 = self._parameters["phi0"] * np.pi / 180
        v = self._parameters["v"]
        nu = self._parameters["nu"]
        x0 = self._parameters["x0"]
        y0 = self._parameters["y0"]

        z0 = F*v/np.sqrt((u0*np.cos(phi0)-nu*x0)**2 + (u0*np.sin(phi0)-nu*y0)**2)
        v_z = -nu*z0
        v_x = -z0/F*(u0*np.cos(phi0)-nu*x0)
        v_y = z0/F*(u0*np.sin(phi0)-nu*y0)
        v_dev = Vector3(v_x,v_y,v_z)
        return v_dev,z0

    def on_horizontal_align(self):
        v_dev, z0 = self.pseudo_acc_solve_dev()
        print("ESTIMATED V_DEV", v_dev)
        Rt_quat = hor_to_dev(self._orientation).conj()
        v_hor = Rt_quat*v_dev
        alt, az = ocef_to_altaz(v_hor, allow_neg=True)
        dat = {
            "direction": {
                "selection_type": "horizon",
                "value": {
                    "az":   az*180/np.pi,
                    "alt": alt*180/np.pi
                },
            }
        }
        self.parameter_form.set_values(dat)
        self.on_commit()

    def on_balance_azimuth(self):
        direction = self._raw_parameters["direction"]
        if direction["selection_type"] != "horizon":
            messagebox.showwarning(title=get_locale("tools.trajectory_calculator.warning.no_hor.title"),
                                   message=get_locale("tools.trajectory_calculator.warning.no_hor.message"),
                                   )
        tres = self._parameters["tres"]  # Temporal resolution
        F = self._orientation["FOCAL_DISTANCE"]
        alt = direction["value"]["alt"]*np.pi/180
        X, Y, u_x, u_y = self.get_kinematics()
        u_x /= tres
        u_y /= tres
        Ax = F*u_y
        Ay = F*u_x
        Az = X*u_y - Y*u_x
        Rt_quat = hor_to_dev(self._orientation).conj()
        A = Vector3(Ax,Ay,Az)
        B = Rt_quat*A
        C1 = B.x*np.cos(alt)
        C2 = B.y*np.cos(alt)
        C3 = B.z*np.sin(alt)
        if C2==C3:
            if C1 == 0:
                return
            t = -C2/C1
        else:
            s = np.sqrt(C1**2-C3**2+C2**2)
            if s<0:
                return
            t = (s-C1)/(C3-C2)
        az = 2*np.arctan(t)*180/np.pi
        dat = {
            "direction": {
                "selection_type": "horizon",
                "value": {
                    "az": az
                },
            }
        }
        self.parameter_form.set_values(dat)
        self.on_commit()




    def get_kinematics(self):
        a = self._parameters["a"] * PIXEL_SIZE
        u0 = self._parameters["u0"] * PIXEL_SIZE
        phi0 = self._parameters["phi0"] * np.pi / 180
        x0 = self._parameters["x0"]
        y0 = self._parameters["y0"]
        k0 = self._parameters["k0"]
        k = self._parameters["k"]
        u_z = -self._parameters["nu"]


        delta_k = k - k0
        X = x0 + (np.cos(phi0) * (u0 * delta_k)) / (1 + u_z * delta_k) + np.cos(phi0) * a * (delta_k ** 2) / 2
        Y = y0 + (np.sin(phi0) * (u0 * delta_k)) / (1 + u_z * delta_k) + np.sin(phi0) * a * (delta_k ** 2) / 2
        u_x = u0 * np.cos(phi0) / (1 + u_z * delta_k) ** 2 + a * delta_k * np.cos(phi0)
        u_y = u0 * np.sin(phi0) / (1 + u_z * delta_k) ** 2 + a * delta_k * np.sin(phi0)
        return X,Y,u_x, u_y