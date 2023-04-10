import tkinter as tk
from tkinter import ttk
import datetime
from vtl_common.common_GUI.modified_base import EntryWithEnterKey
from vtl_common.localization import get_locale

class DatetimeEntry(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.year = tk.StringVar(self)
        self.month = tk.StringVar(self)
        self.day = tk.StringVar(self)
        self.hour = tk.StringVar(self)
        self.minute = tk.StringVar(self)
        self.second = tk.StringVar(self)
        self._index_coonstruct = 0
        self._add_field(self.year)
        self._add_label("/")
        self._add_field(self.month)
        self._add_label("/")
        self._add_field(self.day)
        self._add_label("", 1)
        self._add_field(self.hour)
        self._add_label(":")
        self._add_field(self.minute)
        self._add_label(":")
        self._add_field(self.second)
        self._add_label("", 10)
        self._datetime = None
        self.set_datetime(datetime.datetime.utcnow())
        self.on_commit = None

    def entry_on_commit(self):
        new_d = self.parse_datetime()
        if new_d is not None:
            self._datetime = new_d
        self._sync_display()
        if self.on_commit:
            self.on_commit()

    def _add_field(self, var):
        entry_field = EntryWithEnterKey(self, textvariable=var, width=5)
        entry_field.grid(row=0, column=self._index_coonstruct, sticky="nsew")
        entry_field.on_commit = self.entry_on_commit
        #self.columnconfigure(self._index, weight=1)
        self._index_coonstruct += 1

    def _add_label(self, label, stretchy=0):
        field = tk.Label(self, text=label)
        field.grid(row=0, column=self._index_coonstruct, sticky="nsew")
        if stretchy:
            self.columnconfigure(self._index_coonstruct, weight=stretchy)
        self._index_coonstruct += 1

    def _sync_display(self):
        dt = self._datetime
        self.year.set(str(dt.year))
        self.month.set(str(dt.month))
        self.day.set(str(dt.day))

        self.hour.set(str(dt.hour))
        self.minute.set(str(dt.minute))
        self.second.set(str(dt.second))

    def set_datetime(self, dt:datetime.datetime):
        self._datetime = dt
        self._sync_display()

    def parse_datetime(self):
        try:
            year = int(self.year.get())
            month = int(self.month.get())
            day = int(self.day.get())
            hour = int(self.hour.get())
            minute = int(self.minute.get())
            second = int(self.second.get())
            dt = datetime.datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
            print(dt)
            return dt
        except ValueError:
            return None

    def get_datetime(self):
        return self._datetime


UNITS = [
    ["second", datetime.timedelta(seconds=1)],
    ["minute", datetime.timedelta(minutes=1)],
    ["hour", datetime.timedelta(hours=1)],
    ["day", datetime.timedelta(days=1)],
    #["month",datetime.timedelta(months=1)],   #timedelta does not work with such arguments
    #["year",datetime.timedelta(years=1)],
]

class DatetimeStepper(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        left_button = tk.Button(self, text="<", command=self.step_left)
        left_button.grid(row=0, column=0, sticky="w")

        lab = tk.Label(self, text=get_locale("orientation.step"))
        lab.grid(row=0, column=1, sticky="w")

        values = [get_locale("orientation.step."+x[0]) for x in UNITS]
        self.units = ttk.Combobox(self, values=values, state="readonly")
        self.units.grid(row=0, column=2, sticky="ew")
        self.units.current(0)

        right_button = tk.Button(self, text=">", command=self.step_right)
        right_button.grid(row=0, column=3, sticky="e")
        self.on_step = None

    def step_left(self):
        index = self.units.current()
        self.on_step(-UNITS[index][1])

    def step_right(self):
        index = self.units.current()
        self.on_step(UNITS[index][1])

class DatetimePicker(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.datetime_entry = DatetimeEntry(self)
        self.datetime_entry.pack(side=tk.TOP, fill=tk.X)
        self.datetime_entry.on_commit = self.on_datetime_commit

        self.stepper = DatetimeStepper(self)
        self.stepper.pack(side=tk.BOTTOM, fill=tk.X)
        self.stepper.on_step = self.on_step

        self.on_commit = None
        self.unixtime_limit = None

    def on_datetime_commit(self):
        if self.on_commit:
            self.on_commit()

    def set_limits(self,start,end):
        self.unixtime_limit = start,end

    def set_datetime(self, dt):
        self.datetime_entry.set_datetime(dt)

    def get_datetime(self):
        return self.datetime_entry.get_datetime()

    def on_step(self, step):
        self.set_datetime(self.get_datetime()+step)
        self.on_datetime_commit()

    def get_unixtime(self):
        dt = self.datetime_entry.get_datetime()
        ref_dt = datetime.datetime(1970, 1, 1)
        unixtime = (dt - ref_dt).total_seconds()
        if self.unixtime_limit is not None:
            if unixtime<self.unixtime_limit[0]:
                unixtime = self.unixtime_limit[0]
                self.datetime_entry.set_datetime(ref_dt+datetime.timedelta(seconds=unixtime))
            if unixtime>self.unixtime_limit[1]:
                unixtime = self.unixtime_limit[1]
                self.datetime_entry.set_datetime(ref_dt + datetime.timedelta(seconds=unixtime))
        return unixtime
