import numpy as np

from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode, BoolNode
from vtl_common.localization import get_locale
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from ..orientation.parameters import ParametersForm as OrientationParametersForm


class RightAnscension(FloatNode):
    DISPLAY_NAME = "RA [h]"
    DEFAULT_VALUE = 0.0

    def get_data(self):
        d = super().get_data()
        return d*np.pi/12


class Declination(FloatNode):
    DISPLAY_NAME = "Dec [°]"
    DEFAULT_VALUE = 0.0

    def get_data(self):
        d = super().get_data()
        return d*np.pi/180


class OrientedPoint(FormNode):
    DISPLAY_NAME = ""
    FIELD__show = create_value_field(BoolNode, get_locale("reconstruction.direction.show"), False)
    #FIELD__orientation = OrientationParametersForm
    FIELD__ra = RightAnscension
    FIELD__dec = Declination

