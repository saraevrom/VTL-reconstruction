
from ..tool_base import ToolBase

class OrientationTool(ToolBase):
    TOOL_KEY = "tools.orientation"

    def __init__(self, master):
        super().__init__(master)