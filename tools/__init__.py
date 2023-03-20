from .reconstructor import ReconstructorTool


def add_tools(notebook):
    return [
        ReconstructorTool.add_to_notebook(notebook)
    ]