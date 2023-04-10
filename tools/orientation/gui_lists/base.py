import tkinter as tk
from vtl_common.localization import get_locale


class GUIList(tk.Frame):
    def __init__(self, master, title_key):
        super().__init__(master)
        title_label = tk.Label(self, text=get_locale(title_key))
        title_label.pack(side=tk.TOP, fill=tk.X)
        self._container = []
        self._list_values = tk.Variable(self)
        self._display_list = tk.Listbox(self, selectmode=tk.EXTENDED, listvariable=self._list_values)
        self._display_list.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        button_panel = tk.Frame(self)
        button_panel.pack(side=tk.BOTTOM, fill=tk.X)
        add_btn = tk.Button(button_panel, text="+", command=self.on_add)
        add_btn.pack(side="left")

        remove_btn = tk.Button(button_panel, text="-", command=self.on_remove)
        remove_btn.pack(side="right")
        remove_btn = tk.Button(button_panel, text="E", command=self.clear)
        remove_btn.pack(side="right")

    def obtain_items(self):
        '''
        Obtain new item for list
        :return: item list. Or None if ypu need to cancel addition
        '''
        raise NotImplementedError("Cannot obtain new item")

    def get_items(self):
        return self._container

    @staticmethod
    def represent_item(item):
        raise NotImplementedError("Cannot show items")

    @staticmethod
    def modify_contents(container, new_item):
        raise NotImplementedError("Cannot modify container")

    def clear(self):
        self._display_list.delete(0, tk.END)
        self._container.clear()

    def redraw(self):
        self._display_list.delete(0, tk.END)
        for item in self._container:
            self._display_list.insert(tk.END, self.represent_item(item))

    def on_add(self):
        got_item = self.obtain_items()
        if got_item is None:
            return
        for item in got_item:
            self.modify_contents(self._container, item)
        self.redraw()

    def on_remove(self):
        cur = self._display_list.curselection()
        if cur:
            self._display_list.selection_clear(0, tk.END)
            for i in reversed(cur): #Reverse needed for extended deletion
                self._container.pop(i)
                self._display_list.delete(i)