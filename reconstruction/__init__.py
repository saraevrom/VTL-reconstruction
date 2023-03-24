from .linear_track import LinearTrackModelForm, LinearTrackAltModelForm
from vtl_common.common_GUI.tk_forms_assist import AlternatingNode, FormNode, IntNode, FloatNode, OptionNode, ComboNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale

class ModelSelector(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction_model")
    SEL__linear = LinearTrackModelForm
    SEL__linear_alt = LinearTrackAltModelForm


class OptInt(OptionNode):
    DEFAULT_VALUE = None
    ITEM_TYPE = create_value_field(IntNode, "", 1)


class NUTSSampler(ComboNode):
    SELECTION_READONLY = True
    VALUES = ["pymc", "nutpie", "blackjax", "numpyro"]

SAMPLER_PARAMETERS = [
    (IntNode, "draws", 2000),
    (IntNode, "tune", 2000),
    (IntNode, "chains", 4),
    (IntNode, "random_seed", 5),
    (FloatNode, "target_accept", 0.95),
    (NUTSSampler, "nuts_sampler", "pymc")
]


class Sampler(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.sample")


for par_type, par_name, default_value in SAMPLER_PARAMETERS:
    setattr(Sampler, "FIELD__"+par_name,
            create_value_field(par_type, get_locale(f"reconstruction.sample.{par_name}"), default_value))


class ReconstructionForm(FormNode):
    USE_SCROLLVIEW = True
    FIELD__model = ModelSelector
    FIELD__sampler = Sampler