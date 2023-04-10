import tkinter as tk

from .datetime_picker import DatetimePicker
from .astronomy_display import SkyPlotter
from .source_explorer import SourceExplorer
from .form import OrientationForm
from ..tool_base import ToolBase
from vtl_common.common_GUI.button_panel import ButtonPanel
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localization import get_locale


class OrientationTool(ToolBase):
    TOOL_KEY = "tools.orientation"

    def __init__(self, master):
        super().__init__(master)
        self.datetime_picker = DatetimePicker(self)
        self.datetime_picker.pack(side=tk.TOP, fill=tk.X)
        self.datetime_picker.on_commit = self.on_date_select
        bottom_panel = tk.Frame(self)
        bottom_panel.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.sky_plotter = SkyPlotter(bottom_panel)
        self.sky_plotter.grid(row=0, column=1, sticky="nsew")

        self.source_explorer = SourceExplorer(bottom_panel)
        self.source_explorer.grid(row=0, column=0, sticky="nsew")

        bottom_panel.columnconfigure(0, weight=1)
        bottom_panel.columnconfigure(1, weight=1)
        bottom_panel.rowconfigure(0, weight=1)

        right_panel = tk.Frame(bottom_panel)
        right_panel.grid(row=0, column=2, sticky="nsew")
        self.control_panel = ButtonPanel(right_panel)
        self.control_panel.add_button(text=get_locale("orientation.btn.load_data"),
                                      command=self.on_load_file,
                                      row=0)
        self.control_panel.pack(side=tk.TOP,fill=tk.X)

        self.orientation_form_parser = OrientationForm()
        self.form = TkDictForm(right_panel, self.orientation_form_parser.get_configuration_root(), use_scrollview=True)
        self.form.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.form.on_commit = self.on_form_commit
        self._formdata = None
        self._sync_form()


    def on_load_file(self):
        if self.source_explorer.on_load_file():
            t_start, t_end = self.source_explorer.get_unixtime_interval()
            self.datetime_picker.set_limits(t_start, t_end)
            self.on_form_commit()

    def _sync_form(self):
        formdata = self.form.get_values()
        self.orientation_form_parser.parse_formdata(formdata)
        self._formdata = self.orientation_form_parser.get_data()

    def on_form_commit(self):
        self._sync_form()
        self.on_date_select()

    def on_date_select(self):
        ut = self.datetime_picker.get_unixtime()
        self.sky_plotter.plot_stars(ut)
        self.source_explorer.set_unixtime(ut)