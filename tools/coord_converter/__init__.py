import tkinter as tk

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localization import get_locale
from ..tool_base import ToolBase
from .orientation import Orientation, Location
from .date_time import DatetimeEntry
from .other_coords import CoordList



def create_frame(master, local_key):
    panel = tk.Frame(master)
    tk.Label(panel,text=get_locale(local_key), anchor="w").pack(side="top", fill="x")
    return panel

class CoordConverter(ToolBase):
    TOOL_KEY = "tools.converter"

    def __init__(self, master):
        super().__init__(master)

        fixed_panel = create_frame(self, "converter.fixed_params")
        fixed_panel.grid(row=0, column=0, sticky="nsew")

        mutable_panel = create_frame(self, "converter.coords")
        mutable_panel.grid(row=0, column=1, sticky="nsew")

        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        self.location = Location(fixed_panel)
        self.location.pack(side="top", fill="x")
        self.location.set_commit(self.propagate_fixed)

        self.orientation = Orientation(fixed_panel)
        self.orientation.pack(side="top", fill="both")
        self.orientation.set_commit(self.propagate_fixed)

        self.datetime_entry = DatetimeEntry(fixed_panel)
        self.datetime_entry.pack(side="top", fill="x")
        self.datetime_entry.commit = self.propagate_fixed

        self.coords = CoordList(mutable_panel)
        self.coords.pack(side="top", fill="both", expand=True)

        self.propagate_fixed()

    def propagate_fixed(self):
        orient = self.orientation.get_data()
        loc = self.location.get_data()
        t = self.datetime_entry.get_time()
        self.coords.set_fixed(orient, loc, t)
