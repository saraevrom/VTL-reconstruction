import tkinter as tk
import tkinter.simpledialog as simpledialog
from vtl_common.common_GUI.tk_forms import ScrollView

class SelectionDialog(simpledialog.Dialog):
    def __init__(self, master, options):
        self.selection_options = options
        self.result = None
        self.selector = None
        super().__init__(master)

    def body(self, master: tk.Frame):
        self.selector = tk.Listbox(master)
        for opt in self.selection_options:
            self.selector.insert(tk.END, opt)
        self.selector.pack(fill=tk.BOTH, expand=True)
        self.selector.bind("<Double-Button-1>", self.ok)  # Qol tweak: double-clicking acts as OK button

    def apply(self):
        selection = self.selector.curselection()
        if selection:
            self.result = self.selector.get(selection[0])

class CheckListDialog(simpledialog.Dialog):
    def __init__(self, master, keys):
        self.check_keys = keys
        self.result = None
        self.variables = dict()
        super().__init__(master)

    def body(self, master):
        scrollview = ScrollView(master)
        scrollview.pack(expand=True, fill="both")
        index = 0
        for check in self.check_keys:
            var = tk.IntVar(master)
            self.variables[check] = var
            var.set(0)
            checkbox = tk.Checkbutton(scrollview.contents, text=check, variable=var)
            checkbox.grid(row=index, column=0, sticky="nsw")
            index += 1

    def apply(self):
        self.result = [k for k in self.variables.keys() if self.variables[k].get()]
