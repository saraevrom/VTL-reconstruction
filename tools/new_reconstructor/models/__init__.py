from vtl_common.common_GUI.tk_forms_assist import AlternatingNode
from vtl_common.localization import get_locale
from .linear_track_model import LINEAR_TRACK_FORM
from .track_3d_model import TRACK_3d_FORM


class ModelSelect(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction_model")
    SEL__linear = LINEAR_TRACK_FORM.generate_subform()
    SEL__linear_3d = TRACK_3d_FORM.generate_subform()
