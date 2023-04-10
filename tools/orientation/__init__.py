import tkinter as tk

from .datetime_picker import DatetimePicker
from .astronomy_display import SkyPlotter
from ..tool_base import ToolBase


class OrientationTool(ToolBase):
    TOOL_KEY = "tools.orientation"

    def __init__(self, master):
        super().__init__(master)
        self.datetime_picker = DatetimePicker(self)
        self.datetime_picker.pack(side=tk.TOP, fill=tk.X)
        bottom_panel = tk.Frame(self)
        bottom_panel.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.sky_plotter = SkyPlotter(bottom_panel)
        self.sky_plotter.grid(row=0, column=1, sticky="nsew")
        bottom_panel.columnconfigure(1, weight=1)
        bottom_panel.rowconfigure(0, weight=1)