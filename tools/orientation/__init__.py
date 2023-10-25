import tkinter as tk
from tkinter import ttk
from tkinter.simpledialog import askinteger

import numpy as np
import pymc as pm
import arviz as az
from pandastable import Table

from .datetime_picker import DatetimePicker
from .astronomy_display import SkyPlotter
from .source_explorer import SourceExplorer
from .parameters import ParametersForm
from .form import OrientationForm, MasterCoeffTuner
from ..tool_base import ToolBase
from vtl_common.common_GUI.button_panel import ButtonPanel
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localization import get_locale
from .gui_lists import StarList, TimeList
from .gui_lists.time_list import TimeRange
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.workspace_manager import Workspace
from tools.orientation.orientation.model import create_model
from specific_ui.data_output import DataOutput
from fixed_rotator.astro_math_z_aligned import ocef_to_altaz, Vector3
from common_functions.hor_to_dev import hor_to_dev

ORIENTATION_WORKSPACE = Workspace("orientation")


def modify_gauss_sigma(src, sigma):

    src["sigma"] = sigma
    # if "range" in src.keys():
    #     range_ = src["range"]
    #

def modify_parameter(posterior, post_key, key, key_alt, param_formdata, formdata, inside_alter=False):
    '''

    :param posterior:
    :param post_key: key of parameter in Model
    :param key: key of parameter in ParametersForm
    :param key_alt:  key of parameter in OrientationForm
    :param param_formdata:
    :param formdata:
    :return:
    '''
    if post_key in posterior.keys() and formdata[key_alt] is not None:
        mu, std = get_stats(posterior, post_key)
        print(f"SET {post_key}={mu}±{std}")
        param_formdata[key] = mu
        modify_gauss_sigma(formdata[key_alt], std)
    else:
        print(f"{post_key} is unchanged")

def get_stats(xarr, key):
    mu = xarr[key].median()
    mad = np.median(np.abs(xarr[key] - mu))
    std = mad*np.sqrt(np.pi/2)
    return float(mu), std

class OrientationTool(ToolBase):
    TOOL_KEY = "tools.orientation"

    def __init__(self, master):
        super().__init__(master)
        self._orientation = None
        self._trace = None
        self.datetime_picker = DatetimePicker(self)
        self.datetime_picker.pack(side=tk.TOP, fill=tk.X)
        self.datetime_picker.on_commit = self.on_date_select

        main_panel = tk.Frame(self)
        main_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        main_panel.rowconfigure(0,weight=1)

        result_panel = tk.Frame(main_panel)
        result_panel.grid(row=1, column=0, sticky="nsew",columnspan=3)

        result_notebook = ttk.Notebook(result_panel)
        result_notebook.pack(fill="both",expand=True)

        table_tab = tk.Frame(result_notebook)
        table_tab.pack(fill="both",expand=True)

        self.brief_data = DataOutput(result_notebook)
        self.brief_data.pack(fill="both",expand=True)

        result_notebook.add(table_tab,text="reco table")
        result_notebook.add(self.brief_data,text="brief")

        #result_panel.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        self.result_table = Table(table_tab)
        self.result_table.show()

        self.source_explorer = SourceExplorer(main_panel)
        self.source_explorer.grid(row=0, column=1, sticky="nsew")

        self.sky_plotter = SkyPlotter(main_panel)
        self.sky_plotter.grid(row=0, column=2, sticky="nsew")


        main_panel.columnconfigure(0, weight=1)
        main_panel.columnconfigure(1, weight=1)
        main_panel.columnconfigure(2, weight=1)

        main_panel.rowconfigure(0, weight=1)

        right_panel = tk.Frame(self)
        right_panel.pack(side=tk.RIGHT,fill=tk.Y)
        #right_panel.grid(row=0, column=3, sticky="nsew")
        self.control_panel = ButtonPanel(right_panel)
        self.control_panel.add_button(text=get_locale("orientation.btn.load_data"),
                                      command=self.on_load_file,
                                      row=0)
        self.control_panel.add_button(text=get_locale("orientation.btn.align"),
                                      command=self.on_align,
                                      row=1)
        self.control_panel.add_button(text=get_locale("orientation.btn.accept"),
                                      command=self.on_accept,
                                      row=2)
        self.control_panel.add_button(text=get_locale("orientation.btn.attach_master"),
                                      command=self.on_attach_master,
                                      row=3)
        self.control_panel.pack(side=tk.TOP,fill=tk.X)

        #self.brief_panel = DataOutput(right_panel)
        #self.brief_panel.pack(side=tk.BOTTOM,fill=tk.X)


        self.orientation_form_parser = OrientationForm()
        self.form = TkDictForm(right_panel, self.orientation_form_parser.get_configuration_root(), use_scrollview=True)
        self.form.on_commit = self.on_form_commit

        self.parameters_form = SaveableTkDictForm(right_panel, ParametersForm().get_configuration_root(),
                                                  use_scrollview=True, file_asker=ORIENTATION_WORKSPACE)
        self.parameters_form.pack(side=tk.BOTTOM, fill=tk.X)
        self.parameters_form.on_commit = self.on_parameters_change

        self.form.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)


        left_panel = tk.Frame(main_panel)
        left_panel.grid(row=0, column=0, sticky="nsew")

        self.starlist = StarList(left_panel)
        self.starlist.on_list_changed = self.on_stars_changed
        self.starlist.grid(row=0, column=0, sticky="nsew")

        self.time_list = TimeList(left_panel)
        self.time_list.grid(row=1, column=0, sticky="nsew")
        self.time_list.on_list_changed = self.on_intervals_changed

        left_panel.columnconfigure(0,weight=1)
        left_panel.rowconfigure(0,weight=1)
        left_panel.rowconfigure(1,weight=1)

        self._formdata = None
        self._sync_form()
        self.on_parameters_change()

    def on_stars_changed(self):
        self.source_explorer.set_stars(self.starlist.get_items())
        self.on_form_commit()

    def on_intervals_changed(self):
        iv = self.time_list.get_items()
        self.source_explorer.set_intervals(iv)
        self.on_form_commit()

    def on_load_file(self):
        if self.source_explorer.on_load_file():
            t_start, t_end = self.source_explorer.get_unixtime_interval()
            self.datetime_picker.set_limits(t_start, t_end)
            self.time_list.clear()
            self.time_list.limits = TimeRange.from_unixtime(t_start, t_end)
            self.on_form_commit()

    def on_align(self):
        self._sync_form()
        file = self.source_explorer.get_file()
        stars = self.starlist.get_items()
        times = self.time_list.get_items()
        if file and stars and times and self._orientation and self._formdata:
            model = create_model(file, times, stars, self._orientation,
                                 unixtime=self.datetime_picker.get_unixtime(), tuner=self._formdata["tuner"],
                                 broken=self.source_explorer.get_broken_pixels(),
                                 ffmodel=self.source_explorer.get_ffmodel())
            print("Got model")

            self._trace = self._formdata["sampler"].sample(model)

            df = az.summary(self._trace, stat_focus="median")
            df.insert(0, "parameters", df.index)
            df.reset_index(drop=True)
            self.result_table.model.df = df
            self.result_table.redraw()

            axs = az.plot_trace(self._trace).flatten()
            print(axs)
            fig = axs[0].get_figure()
            fig.tight_layout()
            fig.show()


    def on_attach_master(self):
        self._sync_form()
        self.on_parameters_change()
        param_formdata = self.parameters_form.get_values()
        a = param_formdata["MULTIPLIER"]
        b = param_formdata["OFFSET"]
        print("MASTER COEFFS:")
        coeff = 1/a
        off = -b/a
        print(f"K={coeff}")
        print(f"B={off}")
        self.source_explorer.attach_master_coeff_offset(coeff, off)
        print("MASTER COEFFS attached")

    def on_accept(self):
        if self._trace is None:
            return
        tgt_chain = askinteger("ACCEPT", "chain=")
        if tgt_chain is not None and tgt_chain<self._trace.posterior.chain.shape[0]:
            print("USING chain", tgt_chain)
            workon = self._trace.posterior.isel(chain=tgt_chain)
            outer_formdata = self.form.get_values()
            formdata = outer_formdata["tuner"]
            param_formdata = self.parameters_form.get_values()
            # No need to ensure target formdata != None. If it is None, this variable will be just float constant
            print("AVAILABLE FORMDATA:", formdata)
            if "A" in workon.keys():
                src = formdata["tune_a"]
                mu, std = get_stats(workon, "A")
                if src["selection_type"] != "norm":
                    src["selection_type"] = "norm"
                    src["value"] = {
                        "sigma": std
                    }
                    print("SW A to normal")
                print(f"SET A={mu}±{std}")
                param_formdata["MULTIPLIER"] = mu  # Set center
            else:
                print("A is unchanged")
            modify_parameter(workon, "PSF", "PSF", "tune_psf", param_formdata, formdata)
            modify_parameter(workon, "f", "FOCAL_DISTANCE", "tune_f", param_formdata, formdata)
            modify_parameter(workon, "Ω", "SELF_ROTATION", "tune_rot", param_formdata, formdata)
            modify_parameter(workon, "lat", "VIEW_LATITUDE", "tune_lat", param_formdata, formdata)
            modify_parameter(workon, "lon", "VIEW_LONGITUDE", "tune_lon", param_formdata, formdata)
            self.form.set_values(outer_formdata)
            self.parameters_form.set_values(param_formdata)
            self.on_form_commit()
            self.on_parameters_change() # Trigger event manually

    def _sync_form(self):
        formdata = self.form.get_values()
        self.orientation_form_parser.parse_formdata(formdata)
        self._formdata = self.orientation_form_parser.get_data()

    def _display_brief(self):
        Rt_mat = hor_to_dev(self._orientation).conj()
        hor_vec = Rt_mat*Vector3(0,0,1)
        alt, az = ocef_to_altaz(hor_vec,np)
        alt *= 180/np.pi
        az *= 180/np.pi
        self.brief_data.clear()
        self.brief_data.add_entry("Θ", f"{90-alt:.3f}")
        self.brief_data.add_entry("ALT", f"{alt:.3f}")
        self.brief_data.add_entry("AZ", f"{az:.3f}")


    def on_form_commit(self):
        self._sync_form()
        self.on_date_select()

    def on_date_select(self):
        self.time_list.reference_time = self.datetime_picker.get_datetime()
        ut = self.datetime_picker.get_unixtime()
        self.sky_plotter.plot_stars(ut)
        self.source_explorer.set_unixtime(ut)

    def on_parameters_change(self):
        params = self.parameters_form.get_values()
        self._orientation = params
        self._display_brief()
        self.sky_plotter.draw_fov(params)
        self.source_explorer.set_orientation(params)
        self.on_form_commit()