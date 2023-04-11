from .base import GUIList
from tkinter.simpledialog import Dialog
import tkinter as tk
from vtl_common.datetime_parser import parse_datetimes_dt
from vtl_common.parameters import DATETIME_FORMAT
from vtl_common.localization import get_locale
import datetime

class TimeRange(object):
    def __init__(self, start:datetime.datetime, end:datetime.datetime):
        self.start = start
        self.end = end

    @staticmethod
    def from_unixtime(ut_start, ut_end):
        start = datetime.datetime.utcfromtimestamp(ut_start)
        end = datetime.datetime.utcfromtimestamp(ut_end)
        return TimeRange(start,end)

    def intersects(self, other):
        '''
        Checks if given interval intersects with another one
        :param other:
        :return:
        '''
        if other.start<=self.start<=other.end:
            return True
        if other.start <= self.end <= other.end:
            return True
        if self.start<=other.start<=self.end:
            return True
        if self.start <= other.end <= self.end:
            return True
        return False

    def name(self):
        if self.start == self.end:
            return self.start.strftime(DATETIME_FORMAT)
        else:
            return self.start.strftime(DATETIME_FORMAT)+" — "+self.end.strftime(DATETIME_FORMAT)

    def unify(self,other):
        return TimeRange(min(self.start,other.start), max(self.end, other.end))

    def intersection(self, other):
        '''
        Returns intersection interval in case of intersection with another one
        :param other:
        :return:
        '''
        if self.intersects(other):
            return TimeRange(max(self.start, other.start), min(self.end, other.end))

        return None

    def __lt__(self, other):
        return self.start < other.start

    def __le__(self, other):
        return self.start <= other.start

    def __gt__(self, other):
        return self.start > other.start

    def __ge__(self, other):
        return self.start >= other.start

    def __eq__(self, other):
        return self.start == other.start



class AskTimeInterval(Dialog):
    
    def __init__(self, master, reference_time):
        self.reference_time = reference_time
        self.result = None
        super().__init__(master)
    
    def body(self, master: tk.Frame):
        self.start = tk.StringVar(self)
        self.end = tk.StringVar(self)

        self._add_label(master, text_key="orientation.selection.time.start", row=0)
        self._add_entry(master, self.start, 1)
        self._add_label(master, text_key="orientation.selection.time.end", row=2)
        self._add_entry(master, self.end, 3)
        master.columnconfigure(0, weight=1)

    def apply(self):
        start_dt = parse_datetimes_dt(self.start.get(), self.reference_time)
        end_dt = parse_datetimes_dt(self.end.get(), start_dt)
        self.result = TimeRange(start_dt, end_dt)
        print(self.result.name())


    def _add_label(self, body, text_key, row):
        label = tk.Label(body, text=get_locale(text_key))
        label.grid(row=row, column=0, sticky="ew")

    def _add_entry(self,body,var,row):
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