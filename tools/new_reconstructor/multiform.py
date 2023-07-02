import tkinter as tk
from tkinter import ttk

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.common_GUI.tk_forms_assist import FormNode
from .form_diff import diff, patch

TABS = "MABCD"

class FormUnit(tk.Frame):
    def __init__(self, master, form_conf_cls):
        super().__init__(master)
        self._parser = form_conf_cls()
        self._parser: FormNode
        self._form = TkDictForm(self, self._parser.get_configuration_root())
        self._form.pack(fill="both", expand=True)

    def set_listener(self, on_commit):
        self._form.on_commit = on_commit

    def get_values(self):
        return self._form.get_values()

    def get_data(self):
        v = self.get_values()
        self._parser.parse_formdata(v)
        return self._parser.get_data()

    def set_values(self, v):
        self._form.set_values(v)




class Multiform(tk.Frame):
    def __init__(self, master, form_conf_cls):
        super().__init__(master)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)
        self._tabs = []

        for tab in TABS:
            tab_frame = FormUnit(self, form_conf_cls)
            self.notebook.add(tab_frame, text=tab)
            self._tabs.append(tab_frame)

        self._tabs[0].set_listener(self.on_master_change)
        self._last_mdata = self._tabs[0].get_values()
        self.propagate_master_change = True

    def on_master_change(self):
        if self.propagate_master_change:
            mdata = self._tabs[0].get_values()
            differ = diff(self._last_mdata, mdata)
            for i in range(1, 5):
                tdata = self._tabs[i].get_values()
                patched = patch(tdata, differ)
                self._tabs[i].set_values(patched)
            self._last_mdata = mdata

    def reset_individuals(self):
        mdata = self._tabs[0].get_values()
        for i in range(1, 5):
            self._tabs[i].set_values(mdata)

    def get_values(self):
        res = dict()
        for i, k in enumerate(TABS):
            res[k] = self._tabs[i].get_values()
        return res

    def set_values(self, v):
        for i, k in enumerate(TABS):
            self._tabs[i].set_values(v[k])

    def get_data(self):
        res = dict()
        for i, k in enumerate(TABS):
            res[k] = self._tabs[i].get_data()
        return res
