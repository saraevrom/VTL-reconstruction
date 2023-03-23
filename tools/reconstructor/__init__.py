import arviz

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
from pandastable import Table
import pandas as pd
import pymc as pm
from vtl_common.parameters import NPROC
import arviz as az
import matplotlib.pyplot as plt
from .form import ControlForm
from vtl_common.common_GUI.tk_forms import TkDictForm


TRACKS_WORKSPACE = Workspace("tracks")
# USED_MODEL = linear_track_model_alt

STDEV_PREFIX = "σ_"

def create_records(*args):
    res = dict()
    res["SRC"] = []
    res["PMT"] = []
    res["Chain"] = []
    for arg in args:
        res[str(arg)] = []
        res[STDEV_PREFIX+str(arg)] = []
    return res

def reconstruct_event(src_file, pmt, plot_data, form_data):
    used_model = form_data["model"]
    sampler_params = form_data["sampler"]
    re_model = used_model(plot_data)
    chains = sampler_params["chains"]
    print("Got model")
    with re_model:
        idata_0 = pm.sample(return_inferencedata=True, progressbar=True, cores=NPROC, **sampler_params)
        # idata_0 = pm.sample(2000, chains=chains, tune=2000, random_seed=5, target_accept=0.95,
        #                     return_inferencedata=True, progressbar=True, cores=NPROC)
        # idata_0 = pm.sample(2000, chains=chains, tune=2000, random_seed=5,
        #                     return_inferencedata=True, progressbar=True, cores=NPROC)
        #idata_0 = pm.sample(2000, chains=chains, tune=2000, random_seed=5,
        #                      return_inferencedata=True, progressbar=True)
    # for chain in range(chains):
    #     data["SRC"].append(src_file)
    #     data["PMT"].append(pmt)
    #     data["Chain"].append(chain)
    #     for arg in re_model.unobserved_RVs:
    #         sarg = str(arg)
    #         try:
    #             data[sarg].append(np.mean(np.array(idata_0["posterior"][sarg][chain])))
    #             data[STDEV_PREFIX+sarg].append(np.std(np.array(idata_0["posterior"][sarg][chain])))
    #         except KeyError:
    #             data[sarg].append("-")
    #             data[STDEV_PREFIX + sarg].append("-")
    summary = arviz.summary(idata_0)
    # summary['parameter'] = summary.index
    summary.insert(0,'parameter', summary.index)
    summary = summary.reset_index(drop=True)
    summary.insert(0, 'PMT', pmt)
    summary.insert(0, 'SRC', src_file)
    # summary["SRC"] = src_file
    # summary["PMT"] = pmt
    return idata_0, summary


class ReconstructorTool(ToolBase, PopupPlotable):
    TOOL_KEY = "tools.reconstruction"

    def __init__(self, master):
        super().__init__(master)


        rpanel = tk.Frame(self)
        rpanel.pack(side="right", fill="y")

        self.track_plotter = HighlightingPlotter(self)
        self.track_plotter.pack(fill="both", expand=True)
        self.control_panel = ButtonPanel(rpanel)
        self.control_panel.pack(side="top", fill="x")
        self.control_panel.add_button(get_locale("reconstruction.btn.load"), self.on_load, 0)
        self.control_panel.add_button(get_locale("reconstruction.btn.prev"), self.on_prev, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.next"), self.on_next, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.reconstruct"), self.on_reconstruct, 2)
        self.control_panel.add_button(get_locale("reconstruction.btn.trace"), self.on_traces, 3)


        self.ctrl_form_parser = ControlForm()
        self.ctrl_form = TkDictForm(rpanel, self.ctrl_form_parser.get_configuration_root())
        self.ctrl_form.pack(side="bottom", fill="both")

        bottom_panel = tk.Frame(self)
        bottom_panel.pack(side="bottom", fill="x")

        self.result_table = Table(bottom_panel)
        #self.result_table.show()

        PopupPlotable.__init__(self, self.track_plotter)
        self.loaded_file = None
        self.pointer = 0
        self._filelist = []
        self._length = 0
        self._loaded_data0 = None
        self._loaded_ut0 = None
        self._bottom_left = False
        self._bottom_right = False
        self._top_left = False
        self._top_right = False
        self._bottom_visible = False
        self._traces = []

    def on_traces(self):
        for i, trace in enumerate(self._traces):
            axs = az.plot_trace(trace).flatten()
            print(axs)
            fig = axs[0].get_figure()
            fig.canvas.manager.set_window_title(f'Trace {i}')
            fig.show()

    def on_prev(self):
        self.pointer -= 1
        self.show_event()

    def on_next(self):
        self.pointer += 1
        self.show_event()

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
                self.track_plotter.axes.set_title(filename)
                self.track_plotter.draw()
                self._bottom_left = h5file.attrs["bottom_left"]
                self._bottom_right = h5file.attrs["bottom_right"]
                self._top_left = h5file.attrs["top_left"]
                self._top_right = h5file.attrs["top_right"]

    def on_reconstruct(self):
        formdata = self.ctrl_form.get_values()
        self.ctrl_form_parser.parse_formdata(formdata)
        formdata = self.ctrl_form_parser.get_data()
        used_model = formdata["model"]
        if self.loaded_file:
            # re_model = used_model(self._loaded_data0)
            # data = create_records(*re_model.unobserved_RVs)
            # print("KEYS:", data.keys())
            src_file = self._filelist[self.pointer]
            self._traces.clear()
            dfs = []
            keys = []
            if self._bottom_left:
                print("RECONSTRUCTING BL")
                trace, summary = reconstruct_event(src_file=src_file, pmt="BL", plot_data=self._loaded_data0[:, :8, :8],
                                          form_data=formdata)
                self._traces.append(trace)
                dfs.append(summary)
                keys.append("BL")
            if self._bottom_right:
                print("RECONSTRUCTING BR")
                trace, summary = reconstruct_event(src_file=src_file, pmt="BR", plot_data=self._loaded_data0[:, 8:, :8],
                                          form_data=formdata)
                self._traces.append(trace)
                dfs.append(summary)
                keys.append("BR")
            if self._top_left:
                print("RECONSTRUCTING TL")
                trace, summary = reconstruct_event(src_file=src_file, pmt="TL", plot_data=self._loaded_data0[:, :8, 8:],
                                          form_data=formdata)
                self._traces.append(trace)
                keys.append("TL")
                dfs.append(summary)
            if self._top_right:
                print("RECONSTRUCTING TR")
                trace, summary = reconstruct_event(src_file=src_file, pmt="TR", plot_data=self._loaded_data0[:, 8:, 8:],
                                          form_data=formdata)
                self._traces.append(trace)
                dfs.append(summary)
                keys.append("TR")

            df = pd.concat(dfs)
            #df = pd.DataFrame(data=data)
            self.result_table.model.df = df
            if not self._bottom_visible:
                self.result_table.show()
                self._bottom_visible = True
            else:
                self.result_table.redraw()


