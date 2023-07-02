from vtl_common.common_GUI.tk_forms_assist import AlternatingNode, FormNode, IntNode, FloatNode, OptionNode, ComboNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale



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
