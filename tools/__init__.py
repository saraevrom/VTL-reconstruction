from .reconstructor import ReconstructorTool
from .orientation import OrientationTool


def add_tools(notebook):
    return [
        ReconstructorTool.add_to_notebook(notebook),
        OrientationTool.add_to_notebook(notebook)
    ]