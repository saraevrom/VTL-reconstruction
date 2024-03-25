from vtl_common.common_GUI.tk_forms_assist import FormNode, BoolNode, FloatNode, OptionNode, AlternatingNode, ComboNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale
from common_forms import Sampler
import pymc as pm
from .projectors import ProjectionSelector

def create_range(name_key, min_value, max_value):
    class Range(FormNode):
        DISPLAY_NAME = get_locale(name_key)
        FIELD__min = create_value_field(FloatNode, get_locale("orientation.form.tune.range.min"), min_value)
        FIELD__max = create_value_field(FloatNode, get_locale("orientation.form.tune.range.max"), max_value)

        def get_data(self):
            data = super().get_data()
            return data["min"], data["max"]

    return Range



def wrap_in_option(inner, name_key):
    class TunableParameter(OptionNode):
        DISPLAY_NAME = get_locale(name_key)
        ITEM_TYPE = inner
        def get_data(self):
            data = super().get_data()
            if data is None:
                return lambda x, y: y
            else:
                return data
    return TunableParameter

def create_gauss(default_value, min_value=None, max_value=None):
    if min_value is not None and max_value is not None:
        class TuneForm(FormNode):
            DISPLAY_NAME = ""
            FIELD__sigma = create_value_field(FloatNode, "σ", default_value)
            FIELD__range = create_range("orientation.form.tune.range",min_value, max_value)

            def get_data(self):
                data = super().get_data()
                min_, max_ = data["range"]
                return lambda x, y: pm.TruncatedNormal(x, mu=y, sigma=data["sigma"], lower=min_, upper=max_)
    else:
        class TuneForm(FormNode):
            DISPLAY_NAME = ""
            FIELD__sigma = create_value_field(FloatNode, "σ", default_value)

            def get_data(self):
                data = super().get_data()
                return lambda x, y: pm.Normal(x, mu=y, sigma=data["sigma"])
    return TuneForm


def create_tunable_parameter(name_key, default_value, min_value=None, max_value=None):
    TuneForm = create_gauss(default_value, min_value, max_value)
    return wrap_in_option(TuneForm, name_key)


class MasterCoeffTunerUniform(FormNode):
    #DISPLAY_NAME = get_locale("orientation.form.tune.a")
    DISPLAY_NAME = ""
    FIELD__higher = create_value_field(FloatNode, get_locale("orientation.form.tune.range.max"), 100.0)

    def get_data(self):
        data = super().get_data()
        return lambda x, y: pm.Uniform(x, upper=data["higher"])


class MasterCoeffTuner(AlternatingNode):
    DISPLAY_NAME = get_locale("orientation.form.tune.dist")
    SEL__uniform = MasterCoeffTunerUniform
    SEL__norm = create_gauss(10)


class FinalDistSelect(ComboNode):
    DISPLAY_NAME = get_locale("orientation.form.tune.final_dist")
    VALUES = ["gauss", "laplace", "student"]
    DEFAULT_VALUE = "gauss"
    SELECTION_READONLY = True

class Tuner(FormNode):
    DISPLAY_NAME = get_locale("orientation.form.tune")
    FIELD__tune_lat = create_tunable_parameter("orientation.form.tune.lat", 10.0, -90.0, 90.0)
    FIELD__tune_lon = create_tunable_parameter("orientation.form.tune.lon", 10.0, -180.0, 180.0)
    FIELD__tune_rot = create_tunable_parameter("orientation.form.tune.rot", 10.0, -180.0, 180.0)
    FIELD__tune_f = create_tunable_parameter("orientation.form.tune.f", 10.0, 140.0, 180.0)
    FIELD__tune_psf = create_tunable_parameter("orientation.form.tune.psf", 0.25, 0.1, 3.0)
    FIELD__tune_kappa = wrap_in_option(MasterCoeffTuner, "orientation.form.tune.a")
    FIELD__tune_b = create_tunable_parameter("orientation.form.tune.b", 0.25, 0.1, 3.0)
    FIELD__tune_b_auto_assume = create_value_field(BoolNode, get_locale("orientation.form.tune.b_auto_assume"),True)
    FIELD__uniform_error = create_value_field(BoolNode, get_locale("orientation.form.tune.uniform_error"), True)
    #FIELD__use_laplace = create_value_field(BoolNode, get_locale("orientation.form.tune.use_laplace"), False)
    FIELD__final_dist = FinalDistSelect


class OrientationForm(FormNode):
    FIELD__use_ff = create_value_field(BoolNode, get_locale("orientation.form.use_ff"), True)
    FIELD__tuner = Tuner
    FIELD__sampler = Sampler
    FIELD__projection = ProjectionSelector
    FIELD__scale = create_value_field(FloatNode, "Mag multiplier", 1.0)
