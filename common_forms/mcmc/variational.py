from vtl_common.common_GUI.tk_forms_assist import AlternatingNode, FormNode, IntNode, FloatNode, OptionNode, ComboNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
from .pymc_sample_abstraction import VariationalSample

class VariationalMethod(ComboNode):
    DISPLAY_NAME = "method"
    SELECTION_READONLY = True
    DEFAULT_VALUE = "advi"
    VALUES = [
                "advi",
                "fullrank_advi",
                "svgd",
                "asvgd"
             ]


class SampleForm(FormNode):
    DISPLAY_NAME = "sample"
    FIELD__draws = create_value_field(IntNode, "draws", 1000)
    FIELD__random_seed = create_value_field(IntNode, "random_seed", 5)


class FitForm(FormNode):
    DISPLAY_NAME = "fit"
    FIELD__random_seed = create_value_field(IntNode, "random_seed", 5)
    FIELD__method = VariationalMethod
    FIELD__n = create_value_field(IntNode, "n", 10000)


class VariationalSampler(FormNode):
    DISPLAY_NAME = ""
    FIELD__fit_params = FitForm
    FIELD__sample_params = SampleForm

    def get_data(self):
        data = super().get_data()
        sample_kwargs = data["sample_params"]
        fit_kwargs = data["fit_params"]
        return VariationalSample(fit_kwargs=fit_kwargs, sample_kwargs=sample_kwargs)
