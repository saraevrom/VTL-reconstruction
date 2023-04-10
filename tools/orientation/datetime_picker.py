import tkinter as tk
import datetime
from vtl_common.common_GUI.modified_base import EntryWithEnterKey

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
        self._add_label("", True)
        self._add_field(self.hour)
        self._add_label(":")
        self._add_field(self.minute)
        self._add_label(":")
        self._add_field(self.second)
        self._datetime = None
        self.set_datetime(datetime.datetime.utcnow())

    def entry_on_commit(self):
        new_d = self.parse_datetime()
        if new_d is not None:
            self._datetime = new_d
        self._sync_display()

    def _add_field(self, var):
        entry_field = EntryWithEnterKey(self, textvariable=var, width=5)
        entry_field.grid(row=0, column=self._index_coonstruct, sticky="nsew")
        entry_field.on_commit = self.entry_on_commit
        #self.columnconfigure(self._index, weight=1)
        self._index_coonstruct += 1

    def _add_label(self, label, stretchy=False):
        field = tk.Label(self, text=label)
        field.grid(row=0, column=self._index_coonstruct, sticky="nsew")
        if stretchy:
            self.columnconfigure(self._index_coonstruct, weight=1)
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


class DatetimePicker(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self.datetime_entry = DatetimeEntry(self)
        self.datetime_entry.pack(side=tk.TOP, fill=tk.X)


    def set_datetime(self, dt):
        self.datetime_entry.set_datetime(dt)

    def get_unixtime(self):
        dt = self.datetime_entry.get_datetime()
        unixtime = (dt - datetime.datetime(1970, 1, 1)).total_seconds()
        return unixtime
