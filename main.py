import vtl_common.parameters as parameters

parameters.load_settings()


import tkinter as tk
from tkinter import ttk
import os.path as ospath
from vtl_common import localization
from vtl_common.workspace_manager import Workspace
from vtl_common.parameters import add_parameters_menu
from vtl_common.localization import get_locale
from tools import add_tools

localization.set_locale(parameters.LOCALE)
localization.SEARCH_DIRS.append(ospath.join(ospath.dirname(ospath.abspath(__file__)),"localization"))
parameters.localize_parameters_fields()


class App(tk.Tk):

    def __init__(self):
        super(App, self).__init__()

        self.topmenu = tk.Menu(self)
        self.title(get_locale("app.title"))
        add_parameters_menu(self.topmenu)
        self.config(menu=self.topmenu)

        self.main_notebook = ttk.Notebook(self)
        self.main_notebook.pack(side="top", fill="both", expand=True)
        self.tools = add_tools(self.main_notebook)

        Workspace.initialize_workspace()


if __name__ == "__main__":
    root = App()
    root.mainloop()
