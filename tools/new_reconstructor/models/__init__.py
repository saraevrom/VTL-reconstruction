from vtl_common.common_GUI.tk_forms_assist import AlternatingNode
from vtl_common.localization import get_locale
from .linear_track_model import LinearTrackModel
from .track_3d_model import SpatialTrackModel

#
# class ModelSelect(AlternatingNode):
#     DISPLAY_NAME = get_locale("reconstruction_model")
#     SEL__linear = LINEAR_TRACK_FORM.generate_subform()
#     SEL__linear_3d = TRACK_3d_FORM.generate_subform()

def create_selector():
    class ModelSelect(AlternatingNode):
        DISPLAY_NAME = get_locale("reconstruction_model")
        SEL__linear = LinearTrackModel().generate_subform()
        SEL__linear_3d = SpatialTrackModel().generate_subform()
    return ModelSelect
