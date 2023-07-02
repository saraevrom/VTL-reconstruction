from .reconstructor import ReconstructorTool
from .orientation import OrientationTool
from .new_reconstructor import NewReconstructorTool


def add_tools(notebook):
    return [
        NewReconstructorTool.add_to_notebook(notebook),
        #ReconstructorTool.add_to_notebook(notebook),
        OrientationTool.add_to_notebook(notebook)
    ]
