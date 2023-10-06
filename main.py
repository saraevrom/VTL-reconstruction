#!/usr/bin/env python3

from fastprogress.fastprogress import force_console_behavior
master_bar, progress_bar = force_console_behavior()

from add_parameters import ADD_PARAMS
import vtl_common.parameters as parameters
parameters.add_parameters(*ADD_PARAMS)
parameters.load_settings()

from vtl_common import localization
import os.path as ospath
localization.set_locale(parameters.LOCALE)
localization.SEARCH_DIRS.append(ospath.join(ospath.dirname(ospath.abspath(__file__)),"localization"))
parameters.localize_parameters_fields()


import tkinter as tk
from tkinter import ttk
from vtl_common.workspace_manager import Workspace
from vtl_common.parameters import add_parameters_menu
from vtl_common.localization import get_locale
from tools import add_tools



class App(tk.Tk):

    def __init__(self):
        super(App, self).__init__()

        self.topmenu = tk.Menu(self)
        self.title(get_locale("app.title"))
        add_parameters_menu(self.topmenu)
        self.config(menu=self.topmenu)

        self.main_notebook = ttk.Notebook(self)
        self.main_notebook.pack(side="top", fill="both", expand=True)
        self.tools = add_tools(self.main_notebook,self)

        Workspace.initialize_workspace()

    def ask_from_tool(self, index, what):
        if hasattr(self.tools[index],"ask_"+what):
            asker = getattr(self.tools[index],"ask_"+what)
            return asker()


if __name__ == "__main__":
    root = App()
    root.mainloop()
