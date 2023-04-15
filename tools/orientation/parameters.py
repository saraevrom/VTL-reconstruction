from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode
from vtl_common.parameters import MAIN_LATITUDE, MAIN_LONGITUDE
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
import pymc as pm


PARAM_LIST = [
    ["VIEW_LATITUDE", FloatNode, MAIN_LATITUDE],
    ["VIEW_LONGITUDE", FloatNode, MAIN_LONGITUDE],
    ["SELF_ROTATION", FloatNode, 0.0],
    ["FOCAL_DISTANCE", FloatNode, 165.0],
    ["PSF", FloatNode, 2.85],
    ["MULTIPLIER", FloatNode, 1.0],
    ["OFFSET", FloatNode, 0.0],
]


class ParametersForm(FormNode):
    pass


for param, ptype, default in PARAM_LIST:
    setattr(ParametersForm,"FIELD__"+param, create_value_field(ptype,get_locale("orientation.parameter."+param), default))