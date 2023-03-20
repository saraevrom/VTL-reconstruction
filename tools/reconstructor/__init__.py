from ..tool_base import ToolBase
from .marked_plotter import HighlightingPlotter
from vtl_common.common_GUI.button_panel import ButtonPanel
from vtl_common.localization import get_locale
from vtl_common.workspace_manager import Workspace
from vtl_common.localized_GUI.signal_plotter import PopupPlotable
import tkinter as tk
import zipfile
import h5py
import numpy as np


TRACKS_WORKSPACE = Workspace("tracks")

class ReconstructorTool(ToolBase, PopupPlotable):
    TOOL_KEY = "tools.reconstruction"

    def __init__(self, master):
        super().__init__(master)
        self.track_plotter = HighlightingPlotter(self)
        self.track_plotter.pack(side="left",fill="both", expand=True)
        rpanel = tk.Frame(self)
        rpanel.pack(side="right", fill="y")
        self.control_panel = ButtonPanel(rpanel)
        self.control_panel.pack(side="top", fill="x")
        self.control_panel.add_button(get_locale("reconstruction.btn.load"), self.on_load, 0)
        self.control_panel.add_button(get_locale("reconstruction.btn.prev"), self.on_prev, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.next"), self.on_next, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.reconstruct"), self.on_reconstruct, 2)
        PopupPlotable.__init__(self, self.track_plotter)
        self.loaded_file = None
        self.pointer = 0
        self._filelist = []
        self._length = 0
        self._loaded_data0 = None
        self._loaded_ut0 = None

    def get_plot_data(self):
        if self._loaded_data0 is None:
            return None
        xs = np.arange(self._loaded_data0.shape[0])
        return xs, self._loaded_data0

    def on_load(self):
        filename = TRACKS_WORKSPACE.askopenfilename(auto_formats=["zip"])
        if filename:
            if self.loaded_file:
                self.loaded_file.close()
            self.loaded_file = zipfile.ZipFile(filename, "r")
            self._filelist = self.loaded_file.namelist()
            self._length = len(self._filelist)
            self.pointer = 0
            self.show_event()

    def show_event(self):
        if self.pointer < 0:
            self.pointer = 0
        if self.pointer >= self._length:
            self.pointer = self._length-1
        filename = self._filelist[self.pointer]
        print("Loading", filename)
        with self.loaded_file.open(filename) as fp:
            with h5py.File(fp, "r") as h5file:
                self._loaded_data0 = h5file["data0"][:]
                self._loaded_ut0 = h5file["UT0"][:]
                flattened = np.max(self._loaded_data0, axis=0)
                self.track_plotter.buffer_matrix = flattened
                self.track_plotter.update_matrix_plot(True)
                self.track_plotter.set_mask(h5file.attrs)
                self.track_plotter.draw()

    def on_reconstruct(self):
        if self.loaded_file:
            pass

    def on_prev(self):
        self.pointer -= 1
        self.show_event()

    def on_next(self):
        self.pointer += 1
        self.show_event()

