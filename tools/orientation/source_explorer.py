import json
from pprint import pprint

from vtl_common.localized_GUI import GridPlotter
import tkinter as tk
from vtl_common.workspace_manager import Workspace
import h5py
import numpy as np
import matplotlib.pyplot as plt
from vtl_common.common_flatfielding.models import FlatFieldingModel
from .astronomy_display import DETECTOR_SPAN
from .gui_lists.star_list import StarEntry
from .orientation.stellar_math import unixtime_to_era
from vtl_common.parameters import PIXEL_SIZE
from vtl_common.localized_GUI.signal_plotter import PopupPlotable

from utils import binsearch_tgt

DATA_FILES_WORKSPACE = Workspace("merged_data")
FF_WORKSPACE = Workspace("ff_calibration")




class StarDot(object):
    def __init__(self, star_entry: StarEntry):
        self.star_entry = star_entry
        self._circle = None
        self._text = None

    def draw_at(self, axes, params, unixtime):
        era = unixtime_to_era(unixtime)
        lum = self.star_entry.energy_u()*params["MULTIPLIER"]+params["OFFSET"]
        x, y, visible = self.star_entry.position_on_plane(params, era)


        radius = params["PSF"]*PIXEL_SIZE

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


class SourceExplorer(tk.Frame, PopupPlotable):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
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
        self._intervals = None
        PopupPlotable.__init__(self, self.frame_plotter)


    def set_intervals(self,v):
        self._intervals = v

    def set_orientation(self, v):
        self._orientation = v

    def get_ffmodel(self):
        return self._ffmodel

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
            self._loaded_file = h5py.File(filename, "r+")
            self._ut0 = np.array(self._loaded_file["UT0"])
            if "ffmodel" in self._loaded_file.attrs.keys():
                jsd = json.loads(self._loaded_file.attrs["ffmodel"])
                self.set_ffmodel(jsd)
            return True
        return False

    def attach_master_coeff_offset(self, coeff, offset):
        if (self._ffmodel is not None) and (self._loaded_file is not None):
            jsd = self._ffmodel.dump(coeff, offset)
            pprint(jsd)
            self._loaded_file.attrs["ffmodel"] = json.dumps(jsd)

    def get_broken_pixels(self):
        return self.frame_plotter.get_broken()

    def on_load_ff(self):
        filename = FF_WORKSPACE.askopenfilename(auto_formats=["h5"])
        if filename:
            with open(filename, "r") as fp:
                jsd = json.load(fp)
            self.set_ffmodel(jsd)
            return True
        return False

    def set_ffmodel(self, jsd):
        self._ffmodel = FlatFieldingModel.create_from_parameters(jsd)
        self._ffmodel.master_coeff = 1.0  # Calibration of this coefficient is a purpose of this application.
        self._ffmodel.master_offset = 0.0  # Same
        self.frame_plotter.set_broken(self._ffmodel.broken_query())

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
            index = binsearch_tgt(self._ut0, unixtime)
            #print("FOUND DIFF", index, self._ut0[index] - unixtime)
            return index, self._ut0[index]
        return -1, unixtime

    def get_file(self):
        return self._loaded_file

    def set_stars(self, new_list):
        self._star_list = [StarDot(item) for item in new_list]

    def get_plot_data(self):
        if self._loaded_file and self._intervals:
            datafile = self._loaded_file
            intervals = self._intervals
            ut0 = np.array(datafile["UT0"])
            times = []
            observed = []
            for interval in intervals:
                ut_start, ut_end = interval.unixtime_intervals()
                i_start = binsearch_tgt(ut0, ut_start)
                i_end = binsearch_tgt(ut0, ut_end)
                times.append(ut0[i_start:i_end:interval.stride])
                observed.append(datafile["data0"][i_start:i_end:interval.stride])
                print("INTERVAL SRC", interval.name())
                print(f"INTERVAL {i_start} - {i_end}")
            times = np.concatenate(times)
            observed = np.concatenate(observed)
            if self._ffmodel:
                observed = self._ffmodel.apply(observed)
            return times, observed