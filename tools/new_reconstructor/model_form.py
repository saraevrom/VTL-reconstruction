from vtl_common.common_GUI.tk_forms_assist import FormNode, IntNode, LabelNode, AlternatingNode, FloatNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.common_GUI.tk_forms_assist.factory import kwarg_builder
from vtl_common.localization import get_locale
from common_forms import Sampler
from .cutters import RangeCutter, WholeCutter, ThresholdCutter
from .models import create_selector


@kwarg_builder(RangeCutter)
class RangeSelect(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.range")
    LOCKED = True
    FIELD__start = create_value_field(IntNode, get_locale("reconstruction.form.cutter.start"), 0)
    FIELD__end = create_value_field(IntNode, get_locale("reconstruction.form.cutter.end"), -1)




class WholeSelect(LabelNode):
    DISPLAY_NAME = "---"
    LOCKED = True

    def get_data(self):
        return WholeCutter()


@kwarg_builder(ThresholdCutter)
class ThreshSelect(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.threshold")
    LOCKED = True
    FIELD__threshold = create_value_field(FloatNode, get_locale("reconstruction.form.cutter.threshold"), 3.5)



class CutterSelection(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction.form.cutter.selection")
    LOCKED = True
    SEL__entire = WholeSelect
    SEL__range = RangeSelect
    SEL__threshold = ThreshSelect



def create_reco_params():
    class ReconstructionParameters(FormNode):
        FIELD__cutter = CutterSelection
        FIELD__sampler = Sampler
        FIELD__model = create_selector()

    return ReconstructionParameters