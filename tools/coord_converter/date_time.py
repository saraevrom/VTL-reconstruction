import tkinter as tk
from datetime import datetime

from vtl_common.common_GUI.modified_base import EntryWithEnterKey
from vtl_common.datetime_parser import parse_datetimes_dt
from common_functions import datetime_to_copystr

class DatetimeEntry(tk.Frame):
    def __init__(self, master):
        super().__init__(master)
        self._dtstr = tk.StringVar(self)
        self._time = datetime.utcnow()
        self._dtstr.set(datetime_to_copystr(self._time))
        self._dtstr.trace("w", self.on_dt_change)
        entry = EntryWithEnterKey(self, textvariable=self._dtstr)
        entry.pack(fill="both", expand=True)
        self.commit = None

    def on_dt_change(self,*args):
        dts = self._dtstr.get()
        now = datetime.utcnow()
        self._time = parse_datetimes_dt(dts, now)
        if self.commit is not None:
            self.commit()

    def get_time(self):
        return self._time