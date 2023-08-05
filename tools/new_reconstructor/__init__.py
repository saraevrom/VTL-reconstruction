import json
import os
import tkinter as tk
from tkinter import ttk
import zipfile
from datetime import datetime
import traceback

from pandastable import Table
import h5py
import numpy as np
import pandas as pd
import  arviz


from vtl_common.localized_GUI.signal_plotter import PopupPlotable
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.localization import get_locale
from vtl_common.workspace_manager import Workspace
from vtl_common.common_GUI.button_panel import ButtonPanel
from vtl_common.parameters import DATETIME_FORMAT

from .marked_plotter import HighlightingPlotter
from .orientation_pointing_form import OrientedPoint
from ..tool_base import ToolBase
from .orientation_pointing_form import OrientationParametersForm
from .marked_plotter import NA_TEXT, PlotProxy
from .multiform import Multiform
from .model_form import create_reco_params
from .models.model_base import ReconstructionModelWrapper, ModelWithParameters
from .cutters import Cutter
from .selection_dialog import SelectionDialog, CheckListDialog
from common_functions import ut0_to_copystr

ORIENTATION_WORKSPACE = Workspace("orientation")
TRACKS_WORKSPACE = Workspace("tracks")
RECO_PARAMETERS_WORKSPACE = Workspace("reconstruction_parameters")


def create_checkbox(parent, label_key, variable):
    check = tk.Checkbutton(parent,text=get_locale(label_key), variable=variable, justify="left", anchor="w")
    check.pack(side="bottom", fill="x")


def get_track_attr(fp, attr):
    if attr in fp.attrs.keys():
        return int(fp.attrs[attr])
    else:
        return 1

def render_event(idata_0, split):
    split_chains = split
    #print(idata_0)
    post = idata_0.posterior
    chains = idata_0.posterior["chain"]
    if split_chains:
        summaries = []
        for chain_ in chains:
            chain = int(chain_)
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
        whole_summary = arviz.summary(post,stat_focus="median")
        #posterior_collapsed = post.median(dim="chain")
        #posterior_collapsed = posterior_collapsed.expand_dims(dim={"chain": 1})
        #whole_summary = az.summary(posterior_collapsed)
        whole_summary.insert(0, 'parameter', whole_summary.index)
        whole_summary = whole_summary.reset_index(drop=True)


    return whole_summary


def cut_workon(data, pmt: str):
    if pmt.startswith("M"):
        return data
    if pmt.startswith("A"):
        return data[:, :8, 8:]
    if pmt.startswith("B"):
        return data[:, 8:, 8:]
    if pmt.startswith("C"):
        return data[:, :8, :8]
    if pmt.startswith("D"):
        return data[:, 8:, :8]

class NewReconstructorTool(ToolBase, PopupPlotable):
    TOOL_KEY = "tools.reconstruction"

    def __init__(self, master):
        ToolBase.__init__(self, master)

        self._is_archive = False
        self._length = 0
        self._filelist = []
        self.pointer = 0
        self.loaded_file = None
        self._loaded_data0 = None
        self._loaded_ut0 = None
        self._traces = dict()

        rpanel = tk.Frame(self)
        rpanel.pack(side="right", fill="y")
        rpanel.config(width=500)
        rpanel.pack_propagate(0)

        lpanel = tk.Frame(self)
        lpanel.pack(side="left", fill="y")
        lpanel.config(width=500)
        lpanel.pack_propagate(0)

        self._pointparams = tk.StringVar()
        self._pointparams.set(NA_TEXT)

        angshow = tk.Label(lpanel, textvariable=self._pointparams, anchor="w")
        angshow.pack(side="top", fill="both", anchor="nw")

        self.point_parser = OrientedPoint()
        self.point_form = TkDictForm(lpanel, self.point_parser.get_configuration_root(), False)
        self.point_form.pack(side="top", fill="x")
        self.point_form.on_commit = self.on_orientation_update

        self.orientation_parser = OrientationParametersForm()
        self.orientation_form = SaveableTkDictForm(lpanel, self.orientation_parser.get_configuration_root(), False,
                                                   file_asker=ORIENTATION_WORKSPACE,
                                                   save_label="form.orientation.save",
                                                   load_label="form.orientation.load")
        self.orientation_form.pack(side="top", fill="x")
        self.orientation_form.on_commit = self.on_orientation_update

        self.track_plotter = HighlightingPlotter(self)
        self.track_plotter.pack(fill="both", expand=True)
        PopupPlotable.__init__(self, self.track_plotter, True)

        self.mod_notebook = Multiform(rpanel, create_reco_params)
        self.control_panel = ButtonPanel(rpanel)
        self.control_panel.pack(side="top", fill="x")
        self.control_panel.add_button(get_locale("reconstruction.btn.load"), self.on_load, 0)
        self.control_panel.add_button(get_locale("reconstruction.btn.prev"), self.on_prev, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.next"), self.on_next, 1)
        self.control_panel.add_button(get_locale("reconstruction.btn.reconstruct"), self.on_reconstruct, 2)
        self.control_panel.add_button(get_locale("reconstruction.btn.clear_traces"), self.on_traces_clear, 3)
        self.control_panel.add_button(get_locale("reconstruction.btn.copy_datetime"), self.on_copy_time, 3)
        self.control_panel.add_button(get_locale("reconstruction.btn.trace"), self.on_show_trace, 4)
        self.control_panel.add_button(get_locale("reconstruction.btn.pairplot"), self.on_pair_plot, 4)
        # self.control_panel.add_button(get_locale("reconstruction.btn.reset_individual"),
        #                               self.mod_notebook.reset_individuals, 5)

        self._pmt_a = tk.IntVar(self)
        self._pmt_a.trace("w", self.on_vars_change)

        self._pmt_b = tk.IntVar(self)
        self._pmt_b.trace("w", self.on_vars_change)

        self._pmt_c = tk.IntVar(self)
        self._pmt_c.trace("w", self.on_vars_change)

        self._pmt_d = tk.IntVar(self)
        self._pmt_d.trace("w", self.on_vars_change)

        self._pmt_m = tk.IntVar(self)
        self._pmt_m.trace("w", self.on_vars_change)

        self._use_time = tk.IntVar(self)
        self._use_time.trace("w", self.on_var_invalidate)

        self._split_chains = tk.IntVar()
        self._split_chains.trace("w", self.on_vars_change)

        self._propagate_formchange = tk.IntVar(self)
        self._propagate_formchange.set(1)
        self._propagate_formchange.trace("w", self.on_vars_change)


        self.low_button_panel = ButtonPanel(rpanel)
        self.low_button_panel.add_button(get_locale("reconstruction.btn.load_params"), self.on_parameters_load, 0)
        self.low_button_panel.add_button(get_locale("reconstruction.btn.save_params"), self.on_parameters_save, 0)

        self.low_button_panel.pack(side="bottom", fill="x")

        self.mod_notebook.pack(side="bottom", fill="both", expand=True)

        create_checkbox(rpanel, "reconstruction.use_seconds", self._use_time)
        create_checkbox(rpanel, "reconstruction.form.split_chains", self._split_chains)
        #create_checkbox(rpanel, "reconstruction.propagate_formchange", self._propagate_formchange)
        create_checkbox(rpanel, "reconstruction.all", self._pmt_m)
        create_checkbox(rpanel, "reconstruction.br", self._pmt_d)
        create_checkbox(rpanel, "reconstruction.bl", self._pmt_c)
        create_checkbox(rpanel, "reconstruction.tr", self._pmt_b)
        create_checkbox(rpanel, "reconstruction.tl", self._pmt_a)

        bottom_panel = tk.Frame(self)
        bottom_panel.pack(side="bottom", fill="x")



        self.result_table = Table(bottom_panel)
        self.result_table.show()
        self.on_vars_change()


    def redraw_traces(self):
        buffer_dataframes = []
        for label in self._traces.keys():
            trace = self._traces[label].idata
            summary = render_event(trace, self._split_chains.get())
            summary.insert(0, 'SRC', label)
            buffer_dataframes.append(summary)

        if buffer_dataframes:
            df = pd.concat(buffer_dataframes)
            self.result_table.model.df = df
            self.result_table.redraw()

    def on_copy_time(self):
        if self._loaded_ut0 is not None:
            ut0 = self._loaded_ut0[0]
            s = ut0_to_copystr(ut0)
            self.winfo_toplevel().clipboard_clear()
            self.winfo_toplevel().clipboard_append(s)

    def on_traces_clear(self):
        self._traces.clear()
        self.redraw_traces()
        self.redraw_tracks()

    def on_orientation_update(self):
        if self._loaded_ut0 is not None:
            formdata_orient = self.orientation_form.get_values()
            formdata = self.point_form.get_values()
            self.point_parser.parse_formdata(formdata)
            formdata = self.point_parser.get_data()
            formdata["orientation"] = formdata_orient
            self._pointparams.set(self.track_plotter.set_point_direction(formdata, self._loaded_ut0))
        else:
            self._pointparams.set(NA_TEXT)

    def on_vars_change(self, *args):
        self.track_plotter.set_mask_vars(bl=self._pmt_c,
                                         br=self._pmt_d,
                                         tl=self._pmt_a,
                                         tr=self._pmt_b)
        self.track_plotter.draw()
        self.mod_notebook.propagate_master_change = bool(self._propagate_formchange.get())
        self.redraw_traces()
        self.redraw_tracks()

    def on_var_invalidate(self):
        self.invalidate_popup_plot()

    def get_pmt_modes(self):
        if self._pmt_m.get():
            mode = "M"
            if self._pmt_a.get():
                mode += "A"
            if self._pmt_b.get():
                mode += "B"
            if self._pmt_c.get():
                mode += "C"
            if self._pmt_d.get():
                mode += "D"
            return [mode]
        else:
            res = []
            if self._pmt_a.get():
                res.append("A")
            if self._pmt_b.get():
                res.append("B")
            if self._pmt_c.get():
                res.append("C")
            if self._pmt_d.get():
                res.append("D")
            return res

    def on_parameters_save(self):
        filename = RECO_PARAMETERS_WORKSPACE.asksaveasfilename(auto_formats=["json"])
        if filename:
            data = {
                "a_en": self._pmt_a.get(),
                "b_en": self._pmt_b.get(),
                "c_en": self._pmt_c.get(),
                "d_en": self._pmt_d.get(),
                "m_en": self._pmt_m.get(),
                "propagate": self._propagate_formchange.get(),
                "params": self.mod_notebook.get_values(),
            }
            with open(filename, "w") as fp:
                json.dump(data, fp)

    def on_parameters_load(self):
        filename = RECO_PARAMETERS_WORKSPACE.askopenfilename(auto_formats=["json"])
        if filename:
            with open(filename, "r") as fp:
                data = json.load(fp)
            self.mod_notebook.set_values(data["params"], force=True)
            self._pmt_a.set(data["a_en"])
            self._pmt_b.set(data["b_en"])
            self._pmt_c.set(data["c_en"])
            self._pmt_d.set(data["d_en"])
            self._pmt_m.set(data["m_en"])
            self._propagate_formchange.set(data["propagate"])


    def on_reconstruct(self):
        if self._loaded_data0 is None:
            return
        formdata = self.mod_notebook.get_data()
        for mode in self.get_pmt_modes():
            self._reconstruct(formdata, mode)
        self.on_orientation_update()
        self.invalidate_popup_plot()

    def _reconstruct(self, formdata, pmt):
        try:
            print("RECO MODE", pmt)
            selected = formdata[pmt[0]]
            cutter: Cutter = selected["cutter"]
            sampler_params = selected["sampler"]
            #print("SAMPLER PARAMS", sampler_params)
            model_wrapper: ReconstructionModelWrapper = selected["model"]
            #print("FORMDATA", formdata)

            identifier = f"{self._filelist[self.pointer]}_{pmt[0]}"

            observed = np.array(self._loaded_data0), np.array(self._loaded_ut0)
            start, end = cutter.cut(cut_workon(observed[0], pmt))
            print("CUTTER_WORK",cutter, start,end)
            broken = self.track_plotter.get_broken()
            sampler = model_wrapper.get_pymc_model(observed, start, end, broken, pmt, self)
            sampler.sample(**sampler_params)
            self._traces[identifier] = sampler
            self.redraw_traces()
            self.redraw_tracks()
        except ValueError as e:
            print("Skipped")
            print("Reason: ValueError:", e)
            print("Traceback:")
            print(traceback.format_exc())

        except KeyboardInterrupt:
            print("Skipped")


    def _show_h5(self, h5file, filename):
        self._loaded_data0 = h5file["data0"][:]
        self._loaded_ut0 = h5file["UT0"][:]
        flattened = np.max(self._loaded_data0, axis=0)
        alive_pixels = (flattened != 0)

        self._pmt_c.set(get_track_attr(h5file, "bottom_left"))  # h5file.attrs["bottom_left"]
        self._pmt_d.set(get_track_attr(h5file, "bottom_right"))  # h5file.attrs["bottom_right"]
        self._pmt_a.set(get_track_attr(h5file, "top_left"))  # h5file.attrs["top_left"]
        self._pmt_b.set(get_track_attr(h5file, "top_right"))  # h5file.attrs["top_right"]

        self.track_plotter.buffer_matrix = flattened
        self.track_plotter.alive_pixels_matrix = alive_pixels
        self.track_plotter.update_matrix_plot(True)
        # self.update_track_locations()
        resolution = (self._loaded_ut0[1] - self._loaded_ut0[0]) * 1000

        self.track_plotter.axes.set_title(filename + "\n" +
                                          datetime.utcfromtimestamp(self._loaded_ut0[0]).strftime(
                                              DATETIME_FORMAT) + "\n" +
                                          f"Δt={round(resolution, 3)}ms")

        self.track_plotter.set_origin(0, 0, 0)
        self.track_plotter.clear_added_patches()
        self.on_orientation_update()
        self.track_plotter.draw()

    def on_prev(self):
        self.pointer -= 1
        self.show_event()

    def on_next(self):
        self.pointer += 1
        self.show_event()

    def plot_overlay(self, axes_proxy):
        if self._loaded_data0 is not None:
            self.track_plotter.clear_added_patches()
            modes = self.get_pmt_modes()
            for mode in modes:
                identifier = f"{self._filelist[self.pointer]}_{mode[0]}"
                if identifier in self._traces.keys():
                    obj = self._traces[identifier]
                    obj.reconstruction_overlay(axes_proxy)

    def redraw_tracks(self):
        self.plot_overlay(self.track_plotter)

    def postprocess_auxgrid(self, axes):
        self.plot_overlay(PlotProxy(axes))

    def ask_trace(self):
        options_to_select = self._traces.keys()
        if options_to_select:
            selector = SelectionDialog(self, options_to_select)
            return selector.result
        return None

    def on_show_trace(self):
        label = self.ask_trace()
        if label:
            trace = self._traces[label].idata
            axs = arviz.plot_trace(trace).flatten()
            #print(axs)
            fig = axs[0].get_figure()
            fig.canvas.manager.set_window_title(f'Trace {label}')
            for ax in axs:
                fig = ax.get_figure()
                fig.tight_layout()
            fig = axs[0].get_figure()
            fig.show()

    def on_pair_plot(self):
        label = self.ask_trace()
        if label:
            reco_result:ModelWithParameters = self._traces[label]
            trace = reco_result.idata
            needed_keys = reco_result.get_posterior_variables()
            selected_vars = CheckListDialog(self, needed_keys)
            if selected_vars.result:
                print("KEYS:", needed_keys)
                axs = arviz.plot_pair(trace, group="posterior", var_names=selected_vars.result)
                print(axs)
                if isinstance(axs,np.ndarray):
                    fig = axs[0][0].get_figure()
                else:
                    fig = axs.get_figure()
                fig.canvas.manager.set_window_title(f'Dual {label}')
                # for ax in axs:
                #     fig = ax.get_figure()
                #     fig.tight_layout()
                fig.show()

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

        self.redraw_tracks()
        self.invalidate_popup_plot()


    def on_load(self):
        filename = TRACKS_WORKSPACE.askopenfilename(
            filetypes=[(get_locale("reconstruction.filetypes.data"), "*.zip *.h5")])
        if filename:
            if self.loaded_file:
                self.loaded_file.close()
            if filename.endswith(".zip"):
                self.loaded_file = zipfile.ZipFile(filename, "r")
                self._filelist = self.loaded_file.namelist()
                self._length = len(self._filelist)
                self._is_archive = True
                self.pointer = 0
                self.on_traces_clear()
            elif filename.endswith(".h5"):
                self.loaded_file = h5py.File(filename, "r")
                self._filelist = [os.path.basename(filename)]
                self._length = 1
                self._is_archive = False
                self.pointer = 0
                self.on_traces_clear()
            else:
                raise RuntimeError("File extension must be .zip or .h5")
            self.show_event()

    def get_plot_data(self):
        if self._loaded_data0 is None:
            return None
        if self._use_time.get():
            xs = self._loaded_ut0 - self._loaded_ut0[0]
        else:
            xs = np.arange(self._loaded_data0.shape[0])
        return xs, self._loaded_data0

    def postprocess_plot(self, axes):
        modes = self.get_pmt_modes()
        for mode in modes:
            identifier = f"{self._filelist[self.pointer]}_{mode[0]}"
            if identifier in self._traces.keys():
                obj = self._traces[identifier]
                obj.postprocess(axes)
