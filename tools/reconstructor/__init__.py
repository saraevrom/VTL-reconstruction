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
import arviz as az
import json
import traceback
from xarray import Dataset
import io
import zipfile
import xarray as xr
from .selection_dialog import SelectionDialog
# az.style.use(["arviz-white", "arviz-redish"])
from vtl_common.parameters import DATETIME_FORMAT
# from vtl_common.parameters import NPROC


import matplotlib.pyplot as plt
from .form import ControlForm
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.parameters import HALF_GAP_SIZE, HALF_PIXELS, PIXEL_SIZE
from datetime import datetime
from ..orientation.parameters import ParametersForm as OrientationParametersForm
from .orientation_pointing_form import OrientedPoint

az.rcParams["data.load"] = "eager"


TRACKS_WORKSPACE = Workspace("tracks")
TRACES_WORKSPACE = Workspace("traces")
RECON_WORKSPACE = Workspace("reconstruction_results")
RECO_PARAMETERS_WORKSPACE = Workspace("reconstruction_parameters")
ORIENTATION_WORKSPACE = Workspace("orientation")
# USED_MODEL = linear_track_model_alt

STDEV_PREFIX = "σ_"


MAIN_OFFSET = HALF_PIXELS*PIXEL_SIZE/2+HALF_GAP_SIZE
OFFSETS = {
    "bl": (-MAIN_OFFSET, -MAIN_OFFSET),
    "br": (MAIN_OFFSET, -MAIN_OFFSET),
    "tl": (-MAIN_OFFSET, MAIN_OFFSET),
    "tr": (MAIN_OFFSET, MAIN_OFFSET)
}

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

def get_slice(pmt):
    lower_slice = slice(None, 8)
    upper_slice = slice(8, None)
    if pmt=="bl":
        return (lower_slice, lower_slice)
    if pmt=="br":
        return (upper_slice, lower_slice)
    if pmt=="tl":
        return (lower_slice, upper_slice)
    if pmt=="tr":
        return (upper_slice, upper_slice)

def get_track_attr(fp, attr):
    if attr in fp.attrs.keys():
        return int(fp.attrs[attr])
    else:
        return 1


PMT_MAPPING = {
    "bl":"C",
    "br":"D",
    "tl":"A",
    "tr":"B",
}

def reconstruct_event(form_data, measured_data, pmt, start,end, broken):
    used_model = form_data["model"][PMT_MAPPING[pmt]][1]
    sampler_params = form_data["sampler"]

    # some_matr = np.max(measured_data, axis=0)
    # plt.matshow(some_matr)
    # plt.title(f"To reconstruct:")
    # plt.show()

    re_model = used_model(np.array(measured_data),start,end, broken)
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
        whole_summary = az.summary(post,stat_focus="median")
        #posterior_collapsed = post.median(dim="chain")
        #posterior_collapsed = posterior_collapsed.expand_dims(dim={"chain": 1})
        #whole_summary = az.summary(posterior_collapsed)
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
        self._last_traces = {"bl":None, "br":None, "tl":None, "tr":None}
        rpanel = tk.Frame(self)
        rpanel.pack(side="right", fill="y")
        rpanel.config(width=500)
        rpanel.pack_propagate(0)

        lpanel = tk.Frame(self)
        lpanel.pack(side="left", fill="y")
        lpanel.config(width=500)
        lpanel.pack_propagate(0)



        self.point_parser = OrientedPoint()
        self.point_form = TkDictForm(lpanel, self.point_parser.get_configuration_root(), False)
        self.point_form.pack(side="top", fill="x")
        self.point_form.on_commit = self.on_orienatation_update

        self.orientation_parser = OrientationParametersForm()
        self.orientation_form = SaveableTkDictForm(lpanel, self.orientation_parser.get_configuration_root(), False,
                                                   file_asker=ORIENTATION_WORKSPACE)
        self.orientation_form.pack(side="top", fill="x")
        self.orientation_form.on_commit = self.on_orienatation_update

        self.track_plotter = HighlightingPlotter(self)
        self.track_plotter.pack(fill="both", expand=True)


        self.control_panel = ButtonPanel(rpanel)
        self.control_panel.pack(side="top", fill="x")
        self.control_panel.add_button(get_locale("reconstruction.btn.load"), self.on_load_archive, 0)
        self.control_panel.add_button(get_locale("reconstruction.btn.prev"), self.on_prev, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.next"), self.on_next, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.reconstruct"), self.on_reconstruct, 2)
        self.control_panel.add_button(get_locale("reconstruction.btn.analyze"), self.on_analyze, 2)
        self.control_panel.add_button(get_locale("reconstruction.btn.trace"), self.on_plot_arviz_traces, 3)
        self.control_panel.add_button(get_locale("reconstruction.btn.save_traces"), self.on_save_traces, 4)
        self.control_panel.add_button(get_locale("reconstruction.btn.load_traces"), self.on_load_traces, 4)
        self.control_panel.add_button(get_locale("reconstruction.btn.clear_traces"), self.on_clear_traces, 5)
        self.control_panel.add_button(get_locale("reconstruction.btn.save_df"), self.on_save_dataframe, 6)

        self.ctrl_form_parser = ControlForm()
        self.ctrl_form = SaveableTkDictForm(rpanel, self.ctrl_form_parser.get_configuration_root(),
                                            file_asker=RECO_PARAMETERS_WORKSPACE)
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
        self._orientation = None
        self.eager_reconstruction = False

        create_checkbox(rpanel, "reconstruction.br", self._bottom_right)
        create_checkbox(rpanel, "reconstruction.bl", self._bottom_left)
        create_checkbox(rpanel, "reconstruction.tr", self._top_right)
        create_checkbox(rpanel, "reconstruction.tl", self._top_left)
        self._update_formdata()

    def on_plot_arviz_traces(self):
        options_to_select = self._traces.keys()
        if options_to_select:
            selector = SelectionDialog(self, options_to_select)
            label = selector.result
            if label:
                trace = self._traces[label]
                axs = az.plot_trace(trace).flatten()
                print(axs)
                for ax in axs:
                    fig = ax.get_figure()
                    fig.tight_layout()
                fig = axs[0].get_figure()
                fig.canvas.manager.set_window_title(f'Trace {label}')
                fig.show()

    def on_orienatation_update(self):
        if self._loaded_ut0 is not None:
            formdata_orient = self.orientation_form.get_values()
            formdata = self.point_form.get_values()
            self.point_parser.parse_formdata(formdata)
            formdata = self.point_parser.get_data()
            formdata["orientation"] = formdata_orient
            self.track_plotter.set_point_direction(formdata, self._loaded_ut0)



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

    def postprocess_plot(self, axes):
        cutters = self._formdata["cutter"]
        for pmt in ["bl", "br", "tl", "tr"]:
            i_slice, j_slice = get_slice(pmt)
            if self._last_traces[pmt] is not None:
                start, end = cutters[pmt].cut(self._loaded_data0[:, i_slice, j_slice])
                model_wrapper = self._formdata["model"][PMT_MAPPING[pmt]][0]
                model_wrapper.postprocess(axes, start, end, pmt, self._last_traces[pmt])

    def _clear_last(self):
        for pmt in ["bl", "br", "tl", "tr"]:
            self._last_traces[pmt] = None

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
        self.track_plotter.axes.set_title(filename+"\n"+
                                          datetime.utcfromtimestamp(self._loaded_ut0[0]).strftime(DATETIME_FORMAT)+"\n"+
                                          f"Δt={self._loaded_ut0[1]-self._loaded_ut0[0]}s")

        self.track_plotter.draw()
        self.on_orienatation_update()
        self.track_plotter.set_origin(0, 0, 0)

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

    def _reconstruct_fill(self, pmt):
        i_slice, j_slice = get_slice(pmt)
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
            if trace is None or self.eager_reconstruction:
                reconstruction_data = cutters[pmt].cut(self._loaded_data0[:, i_slice, j_slice])
                if reconstruction_data is not None:
                    start, end = reconstruction_data
                    broken = self.track_plotter.get_broken()[i_slice, j_slice]
                    trace = reconstruct_event(formdata, self._loaded_data0[:, i_slice, j_slice], pmt, start, end, broken)
                    self._traces[trace_identifier] = trace
                else:
                    return
            self._last_traces[pmt] = trace
            model_wrapper = formdata["model"][PMT_MAPPING[pmt]][0]
            assert trace is not None
            model_wrapper.reconstruction_overlay(plotter=self.track_plotter, i_trace=trace, offset=OFFSETS[pmt])
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
        self._clear_last()


    def on_reconstruct(self):
        self.eager_reconstruction = True
        self.on_reconstruction_run()

    def on_analyze(self):
        self.eager_reconstruction = False
        self.on_reconstruction_run()

    def on_reconstruction_run(self):
        self._update_formdata()
        if self.loaded_file:
            # re_model = used_model(self._loaded_data0)
            # data = create_records(*re_model.unobserved_RVs)
            # print("KEYS:", data.keys())

            lower_slice = slice(None, 8)
            upper_slice = slice(8, None)
            self.track_plotter.clear_added_patches()
            self._clear_last()

            if self._bottom_left.get():
                self._reconstruct_fill("bl")
            if self._bottom_right.get():
                self._reconstruct_fill("br")
            if self._top_left.get():
                self._reconstruct_fill("tl")
            if self._top_right.get():
                self._reconstruct_fill("tr")

            self.render_traces()
            self.track_plotter.draw()