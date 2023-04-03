from reconstruction import ReconstructionForm
from vtl_common.common_GUI.tk_forms_assist import FormNode, OptionNode, IntNode, BoolNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
from .cutters import RangeCutter, WholeCutter

class RangeSelect(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.range")
    FIELD__start = create_value_field(IntNode, get_locale("reconstruction.form.cutter.start"), 0)
    FIELD__end = create_value_field(IntNode, get_locale("reconstruction.form.cutter.end"), -1)

    def get_data(self):
        data = super().get_data()
        return RangeCutter(**data)

class RangeOption(OptionNode):
    ITEM_TYPE = RangeSelect

    def get_data(self):
        data = super().get_data()
        if data is None:
            return WholeCutter()
        else:
            return data

class Cutter(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter")
    FIELD__bl = create_value_field(RangeOption, "BL", None)
    FIELD__br = create_value_field(RangeOption, "BR", None)
    FIELD__tl = create_value_field(RangeOption, "TL", None)
    FIELD__tr = create_value_field(RangeOption, "TR", None)

class ControlForm(ReconstructionForm):
    USE_SCROLLVIEW = True
    FIELD__cutter = Cutter
    FIELD__split_chains = create_value_field(BoolNode, get_locale("reconstruction.form.split_chains"), True)
    FIELD__overwrite = create_value_field(BoolNode, get_locale("reconstruction.form.overwrite"), False)