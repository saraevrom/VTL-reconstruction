import datetime
import json
import tkinter as tk

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
from fixed_rotator.astro_math_z_aligned import radec_to_eci, eci_to_ocef
from fixed_rotator import datetime_to_era
from vtl_common.localization import get_locale
from vtl_common.common_GUI.button_panel import ButtonPanel
from .input_form import M_S, M_MRAD, M_MM, M_KM, DT_DEFAULT
from .mode_asker import OptionDialog

W_DETECTOR = get_locale("tools.trajectory_calculator.DETECTOR")

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
        bpanel.add_button(get_locale("form.point.load"), self.on_point_load, 1)
        bpanel.add_button(get_locale("form.point.save"), self.on_point_save, 1)
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
        radiant = self._raw_parameters["radiant"]
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
                "radiant": jsd,
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

        a = self._parameters["a"]*PIXEL_SIZE
        u0 = self._parameters["u0"]*PIXEL_SIZE
        phi0 = self._parameters["phi0"]*np.pi/180
        x0 = self._parameters["x0"]
        y0 = self._parameters["y0"]
        k0 = self._parameters["k0"]
        k = self._parameters["k"]
        u_z = self._parameters["u_z"]
        tres = self._parameters["tres"] # Temporal resolution
        v = self._parameters["v"]

        base_dt = datetime.datetime.now()
        base_dt = base_dt.replace(microsecond=0)
        dt = parse_datetimes_dt(self._parameters["k0_datetime"], base_dt)
        era = datetime_to_era(dt+datetime.timedelta(seconds=(k-k0)*tres))

        delta_k = k-k0
        X = (x0+np.cos(phi0)*(u0*delta_k+a*(delta_k**2)/2))/(1+u_z*delta_k)
        Y = (y0+np.sin(phi0)*(u0*delta_k+a*(delta_k**2)/2))/(1+u_z*delta_k)
        #u = u0+a*delta_k
        u_x = ((u0+a*delta_k+(a*u_z*delta_k**2)/2)*np.cos(phi0)-x0*u_z)/(1+u_z*delta_k)**2
        u_y = ((u0+a*delta_k+(a*u_z*delta_k**2)/2)*np.sin(phi0)-y0*u_z)/(1+u_z*delta_k)**2

        self.output_panel.clear()
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.focal_plane"))
        self.output_panel.add_entry("x",f"{X:.3f}")
        self.output_panel.add_entry("y",f"{Y:.3f}")
        self.output_panel.add_entry("U",f"{((u_x**2+u_y**2)**0.5)/PIXEL_SIZE:.3f}")

        F = self._orientation["FOCAL_DISTANCE"]
        # omega = calculate_omega(tres, F, X, Y, u, phi0)
        omega = calculate_omega_proj(tres, F, X, Y, u_x, u_y)
        wx,wy,wz = (1000*omega).unpack() #mrad/s
        # self.output_panel.add_entry("omega vector (DETECTOR VIEW) [mrad/t]", f"[{wx:.3f}, {wy:.3f}, {wz:.3f}]")
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.omega"))
        self.output_panel.add_entry(f"ω ({W_DETECTOR}), [{M_MRAD}/{M_S}]", f"[{wx:.3f}, {wy:.3f}, {wz:.3f}]")
        self.output_panel.add_entry(f"ω [{M_MRAD}/{M_S}]", f"{1000*omega.length():.3f}")

        self_rot = self._orientation["SELF_ROTATION"]*np.pi/180
        dev_dec = self._orientation["VIEW_LATITUDE"]*np.pi/180
        dev_gha = self._orientation["VIEW_LONGITUDE"]*np.pi/180
        dev_lat = MAIN_LATITUDE*np.pi/180
        dev_lon = MAIN_LONGITUDE*np.pi/180

        # R matrix from article
        R_quat = Quaternion.rotate_xy(self_rot) * \
                 latlon_quaternion(dev_dec,dev_gha).conj() * \
                 latlon_quaternion(dev_lat,dev_lon)

        Rt_quat = R_quat.conj()

        ra = self._parameters["radiant"]["ra"]
        dec = self._parameters["radiant"]["dec"]
        v_eci = -v*radec_to_eci(ra,dec)

        # Velocity in detector frame
        v_dev = eci_to_ocef(era,dev_dec,dev_gha,self_rot) * v_eci

        u_k0_x = u0*np.cos(phi0)-x0*u_z
        u_k0_y = u0*np.sin(phi0)-y0*u_z

        omega_0 = calculate_omega_proj(tres, F, x0, y0, u_k0_x, u_k0_y)

        r_dev = calculate_dev_vector(F, x0, y0, omega_0, v_dev)
        r_dev = r_dev + v_dev*(k-k0)*tres
        r_hor = Rt_quat*r_dev
        h = r_hor.z
        self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.distances"))
        self.output_panel.add_entry(f"v_dev [{M_KM}/{M_S}]", str(v_dev))
        self.output_panel.add_entry(f"dev [{M_KM}]", str(r_dev))
        self.output_panel.add_entry(f"H [{M_KM}]", h)

        # r_dev = calculate_dev_vector(F, X, Y, omega, v_dev)
        # r_hor = Rt_quat * r_dev
        # h = r_hor.z
        # self.output_panel.add_separator(get_locale("tools.trajectory_calculator.section.distances_alt"))
        # self.output_panel.add_entry(f"dev [{M_KM}]", str(r_dev))
        # self.output_panel.add_entry(f"H [{M_KM}]", h)

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