import tkinter as tk

from ..orientation.parameters import ParametersForm
from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.workspace_manager import Workspace
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE

ORIENTATION_WORKSPACE = Workspace("orientation")

class Orientation(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._form = SaveableTkDictForm(self, ParametersForm().get_configuration_root(), use_scrollview=False,
                                        load_label="form.orientation.load",
                                        save_label="form.orientation.save",
                                        file_asker=ORIENTATION_WORKSPACE)
        self._form.pack(fill="both", expand=True)

    def set_commit(self, callback):
        self._form.on_commit = callback

    def get_data(self):
        return self._form.get_values()

class LocationForm(FormNode):
    FIELD__lat = create_value_field(FloatNode, "LAT [°]", MAIN_LATITUDE)
    FIELD__lon = create_value_field(FloatNode, "LON [°]", MAIN_LONGITUDE)

class Location(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._form = TkDictForm(self, LocationForm().get_configuration_root(), use_scrollview=False)
        self._form.pack(fill="both", expand=True)

    def set_commit(self, callback):
        self._form.on_commit = callback

    def get_data(self):
        return self._form.get_values()