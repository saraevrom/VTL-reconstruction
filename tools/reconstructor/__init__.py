import os
import tempfile

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
import json
import traceback
from xarray import Dataset
import io
import zipfile
import xarray as xr
from .selection_dialog import SelectionDialog
# az.style.use(["arviz-white", "arviz-redish"])

from tkinter import messagebox

import matplotlib.pyplot as plt
from .form import ControlForm
from vtl_common.common_GUI.tk_forms import TkDictForm
az.rcParams["data.load"] = "eager"


TRACKS_WORKSPACE = Workspace("tracks")
TRACES_WORKSPACE = Workspace("traces")
RECON_WORKSPACE = Workspace("reconstruction_results")
# USED_MODEL = linear_track_model_alt

STDEV_PREFIX = "σ_"

def noext(x):
    return  os.path.splitext(x)[0]

def create_records(*args):
    res = dict()
    res["SRC"] = []
    res["PMT"] = []
    res["Chain"] = []
    for arg in args:
        res[str(arg)] = []
        res[STDEV_PREFIX+str(arg)] = []
    return res


def find_trace_entry(traces, entry_name):
    print("Finding trace", entry_name)
    if entry_name in traces.keys():
        print("FOUND")
        return traces[entry_name]
    print("NOT FOUND")
    return None


def get_track_attr(fp, attr):
    if attr in fp.attrs.keys():
        return int(fp.attrs[attr])
    else:
        return 1

def reconstruct_event(form_data, measured_data):
    used_model = form_data["model"]
    sampler_params = form_data["sampler"]

    # some_matr = np.max(measured_data, axis=0)
    # plt.matshow(some_matr)
    # plt.title(f"To reconstruct:")
    # plt.show()

    re_model = used_model(np.array(measured_data))
    with re_model:
        print("Sampling")
        idata_0 = pm.sample(return_inferencedata=True, progressbar=True,
                            **sampler_params)
    return idata_0


def render_event(idata_0, form_data):
    split_chains = form_data["split_chains"]
    sampler_params = form_data["sampler"]
    #print(idata_0)
    post: Dataset = idata_0.posterior
    chains = sampler_params["chains"]
    if split_chains:
        summaries = []
        for chain in range(chains):

            chained = post.isel(chain=slice(chain, chain+1)) # Needed for summary
            subsum = arviz.summary(chained)
            subsum.insert(0, 'parameter', subsum.index)
            subsum = subsum.reset_index(drop=True)
            subsum.insert(0, 'chain', chain)
            summaries.append(subsum)

        whole_summary = pd.concat(summaries)
        whole_summary = whole_summary.reset_index(drop=True)

        # summary['parameter'] = summary.index
    else:
        print(post)
        posterior_collapsed = post.median(dim="chain")
        posterior_collapsed = posterior_collapsed.expand_dims(dim={"chain": 1})
        whole_summary = az.summary(posterior_collapsed)
        whole_summary.insert(0, 'parameter', whole_summary.index)
        whole_summary = whole_summary.reset_index(drop=True)


    return whole_summary


def create_checkbox(parent, label_key, variable):
    check = tk.Checkbutton(parent,text=get_locale(label_key), variable=variable, justify="left", anchor="w")
    check.pack(side="bottom", fill="x")

class ReconstructorTool(ToolBase, PopupPlotable):
    TOOL_KEY = "tools.reconstruction"

    def __init__(self, master):
        super().__init__(master)

        self._formdata = None
        rpanel = tk.Frame(self)
        rpanel.pack(side="right", fill="y")

        self.track_plotter = HighlightingPlotter(self)
        self.track_plotter.pack(fill="both", expand=True)
        self.control_panel = ButtonPanel(rpanel)
        self.control_panel.pack(side="top", fill="x")
        self.control_panel.add_button(get_locale("reconstruction.btn.load"), self.on_load_archive, 0)
        self.control_panel.add_button(get_locale("reconstruction.btn.prev"), self.on_prev, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.next"), self.on_next, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.reconstruct"), self.on_reconstruct, 2)
        self.control_panel.add_button(get_locale("reconstruction.btn.trace"), self.on_plot_arviz_traces, 3)
        self.control_panel.add_button(get_locale("reconstruction.btn.save_traces"), self.on_save_traces, 4)
        self.control_panel.add_button(get_locale("reconstruction.btn.save_df"), self.on_save_dataframe, 4)
        self.control_panel.add_button(get_locale("reconstruction.btn.load_traces"), self.on_load_traces, 5)
        self.control_panel.add_button(get_locale("reconstruction.btn.clear_traces"), self.on_clear_traces, 6)

        self.ctrl_form_parser = ControlForm()
        self.ctrl_form = TkDictForm(rpanel, self.ctrl_form_parser.get_configuration_root())
        self.ctrl_form.pack(side="bottom", fill="both",expand=True)

        bottom_panel = tk.Frame(self)
        bottom_panel.pack(side="bottom", fill="x")

        self.result_table = Table(bottom_panel)
        self.result_table.show()

        PopupPlotable.__init__(self, self.track_plotter)
        self.loaded_file = None
        self.pointer = 0
        self._filelist = []
        self._length = 0
        self._loaded_data0 = None
        self._loaded_ut0 = None
        self._bottom_left = tk.IntVar(self)
        self._bottom_left.trace("w", self.on_vars_change)
        self._bottom_right = tk.IntVar(self)
        self._bottom_right.trace("w", self.on_vars_change)
        self._top_left = tk.IntVar(self)
        self._top_left.trace("w", self.on_vars_change)
        self._top_right = tk.IntVar(self)
        self._top_right.trace("w", self.on_vars_change)
        self._traces = dict()
        self._is_archive = False

        create_checkbox(rpanel, "reconstruction.bl", self._bottom_left)
        create_checkbox(rpanel, "reconstruction.br", self._bottom_right)
        create_checkbox(rpanel, "reconstruction.tl", self._top_left)
        create_checkbox(rpanel, "reconstruction.tr", self._top_right)

    def on_plot_arviz_traces(self):
        options_to_select = self._traces.keys()
        if options_to_select:
            selector = SelectionDialog(self, options_to_select)
            label = selector.result
            if label:
                trace = self._traces[label]
                axs = az.plot_trace(trace).flatten()
                print(axs)
                fig = axs[0].get_figure()
                fig.tight_layout()
                fig.canvas.manager.set_window_title(f'Trace {label}')
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

    def on_load_archive(self):
        filename = TRACKS_WORKSPACE.askopenfilename(filetypes=[(get_locale("reconstruction.filetypes.data"), "*.zip *.h5")])
        if filename:
            if self.loaded_file:
                self.loaded_file.close()
            if filename.endswith(".zip"):
                self.loaded_file = zipfile.ZipFile(filename, "r")
                self._filelist = self.loaded_file.namelist()
                self._length = len(self._filelist)
                self._is_archive = True
                self.pointer = 0
            elif filename.endswith(".h5"):
                self.loaded_file = h5py.File(filename, "r")
                self._filelist = [os.path.basename(filename)]
                self._length = 1
                self._is_archive = False
                self.pointer = 0
            else:
                raise RuntimeError("File extension must be .zip or .h5")
            self.show_event()



    def on_load_traces(self):
        self._update_formdata()
        filename = TRACES_WORKSPACE.askopenfilename(auto_formats=["zip"])
        if filename:
            with zipfile.ZipFile(filename, "r", zipfile.ZIP_DEFLATED) as zipf:
                namelist = zipf.namelist()
                self._traces.clear()
                for name in namelist:
                    source_label = noext(name)
                    number, temp_filename = tempfile.mkstemp(suffix=".h5")
                    with open(temp_filename, "wb") as tmp_file:
                        data = zipf.read(name)
                        #print(data)
                        tmp_file.write(data)
                    self._traces[source_label] = az.from_netcdf(temp_filename, engine="h5netcdf")
                    print(self._traces[source_label])
                    os.remove(temp_filename)

                self.render_traces()

    def on_save_traces(self):
        if self._traces:
            filename = TRACES_WORKSPACE.asksaveasfilename(auto_formats=["zip"])
            if filename:
                with zipfile.ZipFile(filename, "w", zipfile.ZIP_DEFLATED) as zipf:
                    for label in self._traces.keys():
                        trace = self._traces[label]
                        trace_name = f"{label}.h5"
                        print("Saving", trace_name)
                        number, temp_filename = tempfile.mkstemp(prefix="saved_", suffix="_"+trace_name)
                        trace.to_netcdf(temp_filename, engine="h5netcdf")
                        print("Saved temporary as", temp_filename)
                        zipf.write(temp_filename, trace_name)
                        os.remove(temp_filename)

    def on_save_dataframe(self):
        filename = RECON_WORKSPACE.asksaveasfilename(auto_formats=["csv"])
        if filename:
            df:pd.DataFrame
            df = self.result_table.model.df
            df.to_csv(filename, sep=",")

    def update_track_locations(self):
        self.track_plotter.set_mask_vars(bl=self._bottom_left,
                                         br=self._bottom_right,
                                         tl=self._top_left,
                                         tr=self._top_right)

    def on_vars_change(self, *args):
        self.update_track_locations()
        self.track_plotter.draw()


    def _show_h5(self, h5file, filename):
        self._loaded_data0 = h5file["data0"][:]
        self._loaded_ut0 = h5file["UT0"][:]
        flattened = np.max(self._loaded_data0, axis=0)

        self._bottom_left.set(get_track_attr(h5file, "bottom_left"))  # h5file.attrs["bottom_left"]
        self._bottom_right.set(get_track_attr(h5file, "bottom_right"))  # h5file.attrs["bottom_right"]
        self._top_left.set(get_track_attr(h5file, "top_left"))  # h5file.attrs["top_left"]
        self._top_right.set(get_track_attr(h5file, "top_right"))  # h5file.attrs["top_right"]

        self.track_plotter.buffer_matrix = flattened
        self.track_plotter.update_matrix_plot(True)
        # self.update_track_locations()
        self.track_plotter.axes.set_title(filename)
        self.track_plotter.draw()

    def show_event(self):
        if self.pointer < 0:
            self.pointer = 0
        if self.pointer >= self._length:
            self.pointer = self._length-1
        filename = self._filelist[self.pointer]
        print("Loading", filename)
        if self._is_archive:
            with self.loaded_file.open(filename) as fp:
                with h5py.File(fp, "r") as h5file:
                    self._show_h5(h5file, filename)
        else:
            self._show_h5(self.loaded_file, filename)

    def _reconstruct_fill(self, pmt, i_slice, j_slice):
        src_file = self._filelist[self.pointer]
        upper_pmt = pmt.upper()
        formdata = self._formdata
        cutters = formdata["cutter"]
        try:
            print("RECONSTRUCTING", upper_pmt)
            # reco_data = cutters[pmt].cut(self._loaded_data0)
            # trace, summary = reconstruct_event_to_remove(self._traces, src_file=src_file, pmt=upper_pmt,
            #                                    plot_data=reco_data[:, i_slice, j_slice],
            #                                    form_data=formdata)
            # output_dfs.append(summary)

            trace_identifier = f"{noext(src_file)}_{pmt}"

            trace = find_trace_entry(self._traces, trace_identifier)
            if trace is None or formdata["overwrite"]:
                reconstruction_data = cutters[pmt].cut(self._loaded_data0)
                if reconstruction_data is not None:
                    trace = reconstruct_event(formdata, reconstruction_data[:, i_slice, j_slice])
                    self._traces[trace_identifier] = trace

            # summary = render_event(trace, formdata)
            # summary.insert(0, 'PMT', pmt)
            # summary.insert(0, 'SRC', src_file)
            # self._buffer_dataframes.append(summary)

        except ValueError as e:
            print("Skipped")
            print("Reason: ValueError:", e)
            print("Traceback:")
            print(traceback.format_exc())

        except KeyboardInterrupt:
            print("Skipped")


    def on_clear_traces(self):
        self._traces.clear()


    def render_traces(self):
        formdata = self._formdata
        buffer_dataframes = []
        for label in self._traces.keys():
            trace = self._traces[label]
            summary = render_event(trace, formdata)
            summary.insert(0, 'SRC', label)
            buffer_dataframes.append(summary)

        if buffer_dataframes:
            df = pd.concat(buffer_dataframes)
            self.result_table.model.df = df
            self.result_table.redraw()


    def _update_formdata(self):
        formdata = self.ctrl_form.get_values()
        print("USED PARAMETERS:", json.dumps(formdata, indent=4, sort_keys=True))
        self.ctrl_form_parser.parse_formdata(formdata)
        self._formdata = self.ctrl_form_parser.get_data()

    def on_reconstruct(self):
        self._update_formdata()
        if self.loaded_file:
            # re_model = used_model(self._loaded_data0)
            # data = create_records(*re_model.unobserved_RVs)
            # print("KEYS:", data.keys())

            lower_slice = slice(None, 8)
            upper_slice = slice(8, None)

            if self._bottom_left.get():
                self._reconstruct_fill("bl", lower_slice, lower_slice)
            if self._bottom_right.get():
                self._reconstruct_fill("br", upper_slice, lower_slice)
            if self._top_left.get():
                self._reconstruct_fill("tl", lower_slice, upper_slice)
            if self._top_right.get():
                self._reconstruct_fill("tr", upper_slice, upper_slice)

            self.render_traces()