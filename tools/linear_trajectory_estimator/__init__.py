import tkinter as tk

import numpy as np

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.workspace_manager import Workspace
from ..tool_base import ToolBase
from .input_form import make_form
from tools.orientation.parameters import ParametersForm
from specific_ui import DataOutput
from vtl_common.parameters import PIXEL_SIZE
from fixed_rotator import Vector2, Vector3
from vtl_common.localization import get_locale
from .input_form import M_S, M_MRAD, M_MM

W_DETECTOR = get_locale("tools.trajectory_calculator.DETECTOR")

ORIENTATION_WORKSPACE = Workspace("orientation")
SKY_CATALOG = Workspace("sky_catalog")

class LinearEstimator(ToolBase):
    TOOL_KEY = "tools.trajectory_calculator"

    def __init__(self, master):
        super().__init__(master)
        rpanel = tk.Frame(self)
        rpanel.pack(side="right",fill="y")

        self.orientation_form = SaveableTkDictForm(rpanel,
                                                   ParametersForm().get_configuration_root(),
                                                   save_label="form.point.save",
                                                   load_label="form.point.load",
                                                   file_asker=ORIENTATION_WORKSPACE
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

    def on_commit(self):
        orientation = self.orientation_form.get_values()
        parameters = self.parameter_form.get_values()
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
        tres = self._parameters["tres"]
        delta_k = k-k0
        X = x0+np.cos(phi0)*(u0*delta_k+a*(delta_k**2)/2)
        Y = y0+np.sin(phi0)*(u0*delta_k+a*(delta_k**2)/2)
        u = u0+a*delta_k

        self.output_panel.clear()
        self.output_panel.add_entry("x",f"{X:.3f}")
        self.output_panel.add_entry("y",f"{Y:.3f}")
        self.output_panel.add_entry("U",f"{u/PIXEL_SIZE:.3f}")

        F = self._orientation["FOCAL_DISTANCE"]
        S = Vector3(-X, Y,0)
        V = Vector3(-u*np.cos(phi0), u*np.sin(phi0),0)
        ez = Vector3(0,0,1)
        omega = F**2/(F**2+S.dot(S))*(1/F**2*S.cross(V)+1/F*ez.cross(V))/tres
        wx,wy,wz = (1000*omega).unpack() #mrad/s
        # self.output_panel.add_entry("omega vector (DETECTOR VIEW) [mrad/t]", f"[{wx:.3f}, {wy:.3f}, {wz:.3f}]")
        self.output_panel.add_entry(get_locale(f"ω ({W_DETECTOR}), [{M_MRAD}/{M_S}]"), f"[{wx:.3f}, {wy:.3f}, {wz:.3f}]")
        self.output_panel.add_entry(get_locale(f"ω [{M_MRAD}/{M_S}]"), f"{1000*omega.length():.3f}")
