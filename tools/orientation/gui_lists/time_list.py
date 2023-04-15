from .base import GUIList
from tkinter.simpledialog import Dialog
import tkinter as tk
from vtl_common.datetime_parser import parse_datetimes_dt
from vtl_common.localization import get_locale
from vtl_common.common_GUI.modified_base import EntryWithEnterKey
from orientation.time_interval import TimeRange



class AskTimeInterval(Dialog):
    
    def __init__(self, master, reference_time):
        self.reference_time = reference_time
        self.result = None
        self._stride = 1
        super().__init__(master)
    
    def body(self, master: tk.Frame):
        self.start = tk.StringVar(self)
        self.end = tk.StringVar(self)
        self.stride = tk.StringVar(self)
        self.stride.set("1")
        self.stride.trace("w",self.validate_integer)

        self._add_label(master, text_key="orientation.selection.time.start", row=0)
        self._add_entry(master, self.start, 1)
        self._add_label(master, text_key="orientation.selection.time.end", row=2)
        self._add_entry(master, self.end, 3)
        self._add_label(master, text_key="orientation.selection.time.stride", row=4)
        self._add_entry(master, self.stride, 5)
        master.columnconfigure(0, weight=1)

    def validate_integer(self,*args):
        v = self.stride.get()
        if not v or v=="-":
            self._stride = 0
            return
        try:
            actual_integer = int(v)
            self._stride = actual_integer
        except ValueError:
            self.stride.set(str(self._stride))


    def apply(self):
        start_dt = parse_datetimes_dt(self.start.get(), self.reference_time)
        end_dt = parse_datetimes_dt(self.end.get(), start_dt)
        stride = self._stride
        if stride<=1:
            stride = 1
        self.result = TimeRange(start_dt, end_dt, stride=stride)
        print(self.result.name())

    def _add_label(self, body, text_key, row):
        label = tk.Label(body, text=get_locale(text_key))
        label.grid(row=row, column=0, sticky="ew")

    def _add_entry(self, body, var, row):
        entry = tk.Entry(body, textvariable=var)
        entry.grid(row=row, column=0, sticky="ew")


class TimeList(GUIList):
    def __init__(self, master):
        super().__init__(master, "orientation.selection.time")
        self.reference_time = None
        self.limits = None

    def obtain_items(self):
        if self.reference_time is None:
            return None
        range_res = AskTimeInterval(self,self.reference_time)
        if range_res.result is None:
            return None
        print("Getting range")
        return [range_res.result]

    @staticmethod
    def modify_contents(container, new_item):
        container.append(new_item)

    @staticmethod
    def represent_item(item):
        return item.name()

    def postprocess(self, container: list):
        container.sort()
        pointer = 0
        while (pointer<len(container)-1) and (len(container)>1):
            first = container[pointer]
            print("PTR", pointer, "LEN", len(container))
            second = container[pointer+1]
            if first.intersects(second):
                second = container.pop(pointer+1)
                container[pointer] = first.unify(second)
                print("Unified intervals")
            else:
                pointer+= 1

        if self.limits:
            for i in range(len(container)-1,-1,-1):
                inters = self.limits.intersection(container[i])
                if inters:
                    container[i] = inters
                else:
                    container.pop(i)