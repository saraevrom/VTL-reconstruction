from ..tool_base import ToolBase
from vtl_common.localized_GUI import GridPlotter

class ReconstructorTool(ToolBase):
    TOOL_KEY = "tools.reconstruction"

    def __init__(self, master):
        super().__init__(master)
        self.track_plotter = GridPlotter(self)
        self.track_plotter.pack(fill="both", expand=True)