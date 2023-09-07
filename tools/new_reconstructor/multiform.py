import tkinter as tk
from tkinter import ttk

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.common_GUI.tk_forms_assist import FormNode
from .form_diff import diff, patch

TABS = "MABCD"

class FormUnit(tk.Frame):
    def __init__(self, master, form_conf_cls, protection):
        super().__init__(master)
        self._parser = form_conf_cls()
        self._parser: FormNode
        self._form = TkDictForm(self, self._parser.get_configuration_root(), protection=protection)
        self._form.pack(fill="both", expand=True)

    def set_listener(self, on_commit):
        self._form.on_commit = on_commit

    def get_values(self):
        return self._form.get_values()

    def get_data(self):
        v = self.get_values()
        self._parser.parse_formdata(v)
        return self._parser.get_data()

    def set_values(self, v, force=False):
        self._form.set_values(v, force=force)




class Multiform(tk.Frame):
    def __init__(self, master, form_conf_cls_creator):
        super().__init__(master)
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)
        self._tabs = []

        for tab in TABS:
            tab_frame = FormUnit(self, form_conf_cls_creator(), protection=(tab != "M"))
            self.notebook.add(tab_frame, text=tab)
            self._tabs.append(tab_frame)

        self._tabs[0].set_listener(self.on_master_change)
        self._last_mdata = self._tabs[0].get_values()
        self.propagate_master_change = True

    def on_master_change(self):
        if self.propagate_master_change:
            pass
            # self.reset_individuals()

    def reset_individuals(self):
        mdata = self._tabs[0].get_values()
        for i in range(1, 5):
            self._tabs[i].set_values(mdata)

    def get_values(self):
        res = dict()
        for i, k in enumerate(TABS):
            res[k] = self._tabs[i].get_values()
        return res

    def set_values(self, v, force=False):
        for i, k in enumerate(TABS):
            self._tabs[i].set_values(v[k], force=force)

    def get_data(self):
        res = dict()
        for i, k in enumerate(TABS):
            res[k] = self._tabs[i].get_data()
        return res

    def activate_tabs_tabs(self,*args):
        current_i = self.notebook.index(self.notebook.select())
        candidate = None
        needs_switch = False
        for i,v in enumerate(args):
            if v:
                self.notebook.tab(i,state="normal")
                if candidate is None:
                    candidate = i
            else:
                self.notebook.tab(i,state="disabled")
                if i==current_i:
                    needs_switch = True

        if needs_switch and candidate is not None:
            self.notebook.select(self._tabs[candidate])
