from .orientation import OrientationTool
from .new_reconstructor import NewReconstructorTool
from .coord_converter import CoordConverter
from .linear_trajectory_estimator import LinearEstimator

def add_tools(notebook,controller):
    return [
        NewReconstructorTool.add_to_notebook(notebook,controller),
        OrientationTool.add_to_notebook(notebook,controller),
        CoordConverter.add_to_notebook(notebook,controller),
        LinearEstimator.add_to_notebook(notebook,controller)
    ]
