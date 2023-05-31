from reconstruction import ReconstructionForm
from vtl_common.common_GUI.tk_forms_assist import FormNode, IntNode, BoolNode, AlternatingNode, LabelNode, FloatNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field, kwarg_builder
from vtl_common.localization import get_locale
from .cutters import RangeCutter, WholeCutter, ThresholdCutter


@kwarg_builder(RangeCutter)
class RangeSelect(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.range")
    FIELD__start = create_value_field(IntNode, get_locale("reconstruction.form.cutter.start"), 0)
    FIELD__end = create_value_field(IntNode, get_locale("reconstruction.form.cutter.end"), -1)




class WholeSelect(LabelNode):
    DISPLAY_NAME = "---"

    def get_data(self):
        return WholeCutter()


@kwarg_builder(ThresholdCutter)
class ThreshSelect(FormNode):
    FIELD__threshold = create_value_field(FloatNode, get_locale("reconstruction.form.cutter.threshold"), 3.5)



class CutterSelection(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.selection")
    SEL__entire = WholeSelect
    SEL__range = RangeSelect
    SEL__threshold = ThreshSelect


class Cutter(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter")
    FIELD__bl = create_value_field(CutterSelection, "BL", None)
    FIELD__br = create_value_field(CutterSelection, "BR", None)
    FIELD__tl = create_value_field(CutterSelection, "TL", None)
    FIELD__tr = create_value_field(CutterSelection, "TR", None)


class ControlForm(ReconstructionForm):
    USE_SCROLLVIEW = True
    FIELD__cutter = Cutter
    FIELD__split_chains = create_value_field(BoolNode, get_locale("reconstruction.form.split_chains"), True)
    FIELD__overwrite = create_value_field(BoolNode, get_locale("reconstruction.form.overwrite"), False)