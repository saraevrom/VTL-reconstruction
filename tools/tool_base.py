import tkinter as tk
from vtl_common.localization import get_locale

class ToolBase(tk.Frame):
    TOOL_KEY="fixme"

    def __init__(self, master):
        super(ToolBase, self).__init__(master)

    @classmethod
    def get_name(cls):
        return get_locale(cls.TOOL_KEY)

    @classmethod
    def add_to_notebook(cls, notebook):
        name = cls.get_name()
        tool_frame = cls(notebook)
        tool_frame.pack(fill="both", expand=True)
        notebook.add(tool_frame, text=name)
        return tool_frame