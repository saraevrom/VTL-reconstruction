import tkinter as tk
import tkinter.simpledialog as simpledialog


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
