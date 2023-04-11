import json

from vtl_common.localized_GUI import GridPlotter
import tkinter as tk
from vtl_common.workspace_manager import Workspace
import h5py
import numpy as np
from vtl_common.common_flatfielding.models import FlatFieldingModel
from .astronomy_display import DETECTOR_SPAN
from .gui_lists.star_list import StarEntry
import numba as nb
from orientation.stellar_math import unixtime_to_era
from orientation.database_reader import VEGA_LUM
import matplotlib.pyplot as plt

DATA_FILES_WORKSPACE = Workspace("merged_data")
FF_WORKSPACE = Workspace("ff_calibration")


@nb.njit()
def binsearch(array):
    n = array.shape[0]
    start = array[0]
    end = array[-1]
    if start>=0 and end>=0:
        if start<=end:
            return 0
        else:
            return n-1
    elif start<=0 and end<=0:
        if start>=end:
            return 0
        else:
            return n-1
    else:
        # Binary search itself
        start_i = 0
        end_i = n-1
        middle_i = (start_i+end_i)//2
        while start_i!=middle_i:
            if array[middle_i] == 0:
                return middle_i
            if array[start_i]*array[middle_i]<0:
                end_i = middle_i
            elif array[end_i]*array[middle_i]<0:
                start_i = middle_i
            else:
                raise RuntimeError("How did we get there?")
            middle_i = (start_i + end_i) // 2
        return middle_i

class StarDot(object):
    def __init__(self, star_entry: StarEntry):
        self.star_entry = star_entry
        self._circle = None
        self._text = None

    def draw_at(self, axes, params, unixtime):
        era = unixtime_to_era(unixtime)
        lum = self.star_entry.energy()*params["MULTIPLIER"]+params["OFFSET"]
        x, y, visible = self.star_entry.position_on_plane(params, era)


        radius = params["PSF"]*3

        if self._circle is None:
            if visible:
                self._circle = plt.Circle((x, y), radius, color="yellow")
                axes.add_patch(self._circle)
        else:
            self._circle.set_center((x,y))
            self._circle.set_radius(radius)
            self._circle.set_visible(visible)

        if self._text is None:
            if visible:
                self._text = axes.text(x, y, self.star_entry.primary_name())
        else:
            self._text.set_position((x,y))
            self._text.set_text(self.star_entry.primary_name())
            self._text.set_visible(visible)
    def __del__(self):
        print()
        if self._circle is not None:
            self._circle.remove()
        if self._text is not None:
            self._text.remove()


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
        self._stars = None
        self._star_list = []
        self._star_plot = None
        self._orientation = None

    def set_orientation(self, v):
        self._orientation = v

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
        if self._orientation is not None and self._star_list:
            for star in self._star_list:
                star.draw_at(self.frame_plotter.axes, self._orientation, self._unixtime)
        self.frame_plotter.draw()

    def get_frame_by_unixtime(self, unixtime):
        if self._loaded_file:
            #index = np.argmin(np.abs(self._ut0-unixtime))
            seeking = self._ut0-unixtime
            index = binsearch(seeking)
            print("FOUND DIFF", seeking[index])
            return index, self._ut0[index]
        return -1, unixtime

    def set_stars(self, new_list):
        self._star_list = [StarDot(item) for item in new_list]