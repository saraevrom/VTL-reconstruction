from .base import GUIList


class TimeList(GUIList):
    def __init__(self, master):
        super().__init__(master, "orientation.selection.time")