from .linear_constant_track import LINEAR_TRACK_CONSTANT_BRIGHTNESS
from .linear_track_hard import LINEAR_HARD_MODEL
from vtl_common.common_GUI.tk_forms_assist import AlternatingNode, FormNode, IntNode, FloatNode, OptionNode, ComboNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale

class ModelSelector(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction_model")
    SEL__linear_constant = LINEAR_TRACK_CONSTANT_BRIGHTNESS.generate_subform()
    #SEL__linear_hard = LINEAR_HARD_MODEL.generate_subform()

    def get_data(self):
        data = super().get_data()
        return {
            "A": data,
            "B": data,
            "C": data,
            "D": data,
        }

class IndependentSelector(FormNode):
    DISPLAY_NAME = get_locale("reconstruction_model")
    FIELD__A = create_value_field(ModelSelector,"A")
    FIELD__B = create_value_field(ModelSelector,"B")
    FIELD__C = create_value_field(ModelSelector,"C")
    FIELD__D = create_value_field(ModelSelector,"D")


class AlternativeModelSelector(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction")
    SEL__common = ModelSelector
    SEL__independent = IndependentSelector

class OptInt(OptionNode):
    DEFAULT_VALUE = None
    ITEM_TYPE = create_value_field(IntNode, "", 1)


class NUTSSampler(ComboNode):
    SELECTION_READONLY = True
    VALUES = ["pymc", "nutpie", "blackjax", "numpyro"]

class NUTSInit(ComboNode):
    SELECTION_READONLY = True
    VALUES = [
              "auto",
              "adapt_diag",
              "jitter+adapt_diag",
              "jitter+adapt_diag_grad",
              "advi+adapt_diag",
              "advi",
              "advi_map",
              "map",
              "adapt_full",
              "jitter+adapt_full",
             ]


SAMPLER_PARAMETERS = [
    (IntNode, "draws", 2000),
    (IntNode, "tune", 2000),
    (IntNode, "chains", 4),
    (IntNode, "random_seed", 5),
    (FloatNode, "target_accept", 0.95),
    (NUTSSampler, "nuts_sampler", "pymc"),
    (NUTSInit, "init", "auto"),
]


class Sampler(FormNode):
    DISPLAY_NAME = get_locale("reconstruction.sample")


for par_type, par_name, default_value in SAMPLER_PARAMETERS:
    setattr(Sampler, "FIELD__"+par_name,
            create_value_field(par_type, get_locale(f"reconstruction.sample.{par_name}"), default_value))


class ReconstructionForm(FormNode):
    USE_SCROLLVIEW = True
    FIELD__model = AlternativeModelSelector
    FIELD__sampler = Sampler