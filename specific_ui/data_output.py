import tkinter as tk
from vtl_common.common_GUI.tk_forms import ScrollView

class DataOutput(tk.Frame):
    def __init__(self,master):
        super().__init__(master)
        self.scrollview = ScrollView(self)
        self.scrollview.pack(fill="both",expand=True)
        self.entries = []
        self.columnconfigure(1,weight=1)

    def add_entry(self, label_text, value):
        label = tk.Label(self.scrollview.contents, text=label_text,padx=10,anchor="w")
        row = len(self.entries)
        label.grid(row=row,column=0,sticky="nsew")
        display = tk.Label(self.scrollview.contents, text=str(value), bg="#FFFFFF", width=30, anchor="e")
        display.grid(row=row,column=1,sticky="nsew")
        self.entries.append((label, display))

    def add_separator(self, label_text):
        label = tk.Label(self.scrollview.contents, text=label_text, padx=10, anchor="w", font="TkDefaultFont 10 bold")
        row = len(self.entries)
        label.grid(row=row, column=0,columnspan=2,sticky="nsew")
        self.entries.append((label, None))

    def clear(self):
        for left,right in self.entries:
            left.destroy()
            if right is not None:
                right.destroy()
        self.entries.clear()

    def add_dict(self, entries:dict):
        for (k,v) in enumerate(entries):
            self.add_entry(k,v)
