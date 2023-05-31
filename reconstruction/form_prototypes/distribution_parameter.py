from vtl_common.common_GUI.tk_forms_assist import FormNode, IntNode, FloatNode, BoolNode, OptionNode, AlternatingNode
from vtl_common.common_GUI.tk_forms_assist.factory import kwarg_builder, create_value_field
import pymc as pm

MU = create_value_field(FloatNode,"mu",0.0)
SIGMA = create_value_field(FloatNode,"sigma",1.0)
LOWER = create_value_field(FloatNode,"lower",0.0)
UPPER = create_value_field(FloatNode,"upper",1.0)
VALUE = create_value_field(FloatNode,"value",1.0)

class LowerOpt(OptionNode):
    DISPLAY_NAME = "lower"
    ITEM_TYPE = VALUE

class UpperOpt(OptionNode):
    DISPLAY_NAME = "upper"
    ITEM_TYPE = VALUE

class DistFactory(object):

    def __init__(self, dist_cls):
        self.dist_cls = dist_cls

    def __call__(self, **kwargs):
        def wrapper(name):
            return self.dist_cls(name, **kwargs)

        return wrapper


@kwarg_builder(DistFactory(pm.Normal))
class NormalBuilder(FormNode):
    DISPLAY_NAME = "Normal"
    FIELD__mu = MU
    FIELD__sigma = SIGMA

@kwarg_builder(DistFactory(pm.Cauchy))
class CauchyBuilder(FormNode):
    DISPLAY_NAME = "Normal"
    FIELD__mu = MU
    FIELD__sigma = SIGMA


@kwarg_builder(DistFactory(pm.HalfNormal))
class HalfNormalBuilder(FormNode):
    DISPLAY_NAME = "HalfNormal"
    FIELD__sigma = SIGMA


@kwarg_builder(DistFactory(pm.Uniform))
class UniformBuilder(FormNode):
    DISPLAY_NAME = "Uniform"
    FIELD__upper = UPPER
    FIELD__lower = LOWER

class TruncatedBuilder(FormNode):
    DISPLAY_NAME = "Truncated"
    FIELD__lower = LowerOpt
    FIELD__upper = UpperOpt


class DistBuilder(AlternatingNode):
    DISPLAY_NAME = "dist"
    SEL__normal = NormalBuilder
    SEL__halfnormal = HalfNormalBuilder
    SEL__uniform = UniformBuilder
    SEL__truncated = TruncatedBuilder
    SEL__cauchy = CauchyBuilder


TruncatedBuilder.FIELD__dist = DistBuilder
