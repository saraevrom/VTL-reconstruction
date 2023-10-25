import tkinter as tk
from vtl_common.common_GUI.tk_forms import ScrollView


H_STYLE_FMT = "TkDefaultFont {size} bold"



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
        display = tk.Entry(self.scrollview.contents, bg="#FFFFFF", width=30)
        display.insert(0,str(value))
        display.configure(state="readonly")
        display.grid(row=row,column=1,sticky="nsew")
        self.entries.append((label, display))

    def add_separator(self, label_text, size=10):
        font = H_STYLE_FMT.format(size=size)
        label = tk.Label(self.scrollview.contents, text=label_text, padx=10, anchor="w", font=font)
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
