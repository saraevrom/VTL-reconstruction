import tkinter as tk
from tkinter.simpledialog import Dialog

class OptionDialog(Dialog):
    def __init__(self, parent, options):
        self.result = None
        self.options = options
        self.index_var = None
        super().__init__(parent)

    def body(self, parent:tk.Frame):
        self.index_var = tk.IntVar(parent)
        parent.columnconfigure(0,weight=1)
        for i, opt in enumerate(self.options):
            radio = tk.Radiobutton(parent,variable=self.index_var,value=i, text=opt)
            radio.grid(row=i,column=0)
            parent.rowconfigure(i,weight=1)

    def apply(self):
        self.result = self.options[self.index_var.get()]