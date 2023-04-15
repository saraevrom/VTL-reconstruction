import tkinter as tk

import pymc as pm
import arviz as az
from pandastable import Table

from .datetime_picker import DatetimePicker
from .astronomy_display import SkyPlotter
from .source_explorer import SourceExplorer
from .parameters import ParametersForm
from .form import OrientationForm
from ..tool_base import ToolBase
from vtl_common.common_GUI.button_panel import ButtonPanel
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localization import get_locale
from .gui_lists import StarList, TimeList
from .gui_lists.time_list import TimeRange
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.workspace_manager import Workspace
from orientation.model import create_model

ORIENTATION_WORKSPACE = Workspace("orientation")


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
        #result_panel.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        self.result_table = Table(result_panel)
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
        self.control_panel.pack(side=tk.TOP,fill=tk.X)

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
            with model:
                self._trace = pm.sample(**self._formdata["sampler"])
                df = az.summary(self._trace)
                df.insert(0, "parameters", df.index)
                df.reset_index(drop=True)
                self.result_table.model.df = df
                self.result_table.redraw()

                axs = az.plot_trace(self._trace).flatten()
                print(axs)
                fig = axs[0].get_figure()
                fig.tight_layout()
                fig.show()


    def on_accept(self):
        pass

    def _sync_form(self):
        formdata = self.form.get_values()
        self.orientation_form_parser.parse_formdata(formdata)
        self._formdata = self.orientation_form_parser.get_data()

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
        self.sky_plotter.draw_fov(params)
        self.source_explorer.set_orientation(params)
        self.on_form_commit()