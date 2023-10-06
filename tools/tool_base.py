import tkinter as tk
from vtl_common.localization import get_locale

class ToolBase(tk.Frame):
    TOOL_KEY="fixme"

    def __init__(self, master):
        super(ToolBase, self).__init__(master)
        self.controller = None

    def connect_controller(self,controller):
        self.controller = controller

    @classmethod
    def get_name(cls):
        return get_locale(cls.TOOL_KEY)

    @classmethod
    def add_to_notebook(cls, notebook,controller):
        name = cls.get_name()
        tool_frame = cls(notebook)
        tool_frame.connect_controller(controller)
        tool_frame.pack(fill="both", expand=True)
        notebook.add(tool_frame, text=name)
        return tool_frame

    def ask_from_tool(self, index, what):
        return self.controller.ask_from_tool(index, what)