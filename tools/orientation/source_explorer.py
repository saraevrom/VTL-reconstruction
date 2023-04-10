import json

from vtl_common.localized_GUI import GridPlotter
import tkinter as tk
from vtl_common.workspace_manager import Workspace
import h5py
import numpy as np
from vtl_common.common_flatfielding.models import FlatFieldingModel
from .astronomy_display import DETECTOR_SPAN

DATA_FILES_WORKSPACE = Workspace("merged_data")
FF_WORKSPACE = Workspace("ff_calibration")

class SourceExplorer(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.frame_plotter = GridPlotter(self)
        self.frame_plotter.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.frame_plotter.axes.arrow(x=0.0, y=0.0, dx=DETECTOR_SPAN, dy=0.0, color="red")
        self.frame_plotter.axes.arrow(x=0.0, y=0.0, dx=0.0, dy=DETECTOR_SPAN, color="blue")
        self._loaded_file = None
        self._ffmodel = None
        self._index = 0
        self._unixtime = 0

    def get_unixtime_interval(self):
        if self._loaded_file:
            return self._ut0[0], self._ut0[-1]

    def set_unixtime(self, unixtime):
        self._index, self._unixtime = self.get_frame_by_unixtime(unixtime)
        self.update_plot()

    def on_load_file(self):
        filename = DATA_FILES_WORKSPACE.askopenfilename(auto_formats=["h5"])
        if filename:
            if self._loaded_file is not None:
                self._loaded_file.close()
            self._loaded_file = h5py.File(filename,"r")
            self._ut0 = np.array(self._loaded_file["UT0"])
            if "ffmodel" in self._loaded_file.attrs.keys():
                jsd = json.loads(self._loaded_file.attrs["ffmodel"])
                self._ffmodel = FlatFieldingModel.create_from_parameters(jsd)
                self.frame_plotter.set_broken(self._ffmodel.broken_query())
            return True
        return False

    def on_load_ff(self):
        filename = FF_WORKSPACE.askopenfilename(auto_formats=["h5"])
        if filename:
            with open(filename, "r") as fp:
                jsd = json.load(fp)
            self._ffmodel = FlatFieldingModel.create_from_parameters(jsd)
            self.frame_plotter.set_broken(self._ffmodel.broken_query())
            return True
        return False

    def update_plot(self, use_ff=True):
        if self._loaded_file and self._index>=0:
            frame = self._loaded_file["data0"][self._index]
            if use_ff and self._ffmodel:
                frame = self._ffmodel.apply(frame)
            self.frame_plotter.buffer_matrix = frame
            self.frame_plotter.update_matrix_plot(True)
            self.frame_plotter.draw()

    def get_frame_by_unixtime(self, unixtime):
        if self._loaded_file:
            index = np.argmin(np.abs(self._ut0-unixtime))
            return index, self._ut0[index]
        return -1, 0.0