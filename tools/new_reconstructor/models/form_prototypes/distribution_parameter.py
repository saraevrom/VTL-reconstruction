import pymc

from vtl_common.common_GUI.tk_forms_assist import FormNode, IntNode, FloatNode, BoolNode, OptionNode, AlternatingNode
from vtl_common.common_GUI.tk_forms_assist.factory import kwarg_builder, create_value_field
import pymc as pm
import pytensor.tensor as pt

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


class DistParameterCreator(object):
    def __call__(self, name, const_storage=None):
        raise NotImplementedError()

    def get_dist(self):
        raise NotImplementedError()

class DistCreator(DistParameterCreator):
    def __init__(self, dist_cls,negative=False, **kwargs):
        self.kwargs = kwargs
        self.dist_cls = dist_cls
        self.apply_negate = negative

    def __call__(self, name, const_storage=None):
        if self.apply_negate:
            neg = self.dist_cls(name+"_neg_", **self.kwargs)
            res = pm.Deterministic(name,-neg)
            return res
        return self.dist_cls(name, **self.kwargs)

    def get_dist(self):
        return self.dist_cls.dist(**self.kwargs)


class ConstantDistCreator(DistParameterCreator):
    def __init__(self, value):
        self.value = value

    def __call__(self, name, const_storage=None):
        # if const_storage is not None:
        #     const_storage[name] = self.value
        tensor = pt.constant(self.value)
        return pymc.Deterministic(name,tensor)
        #return self.value

    def get_dist(self):
        raise RuntimeError("This is constant")

class DistFactory(object):

    def __init__(self, dist_cls):
        self.dist_cls = dist_cls

    def __call__(self, **kwargs):
        return DistCreator(self.dist_cls, **kwargs)




@kwarg_builder(DistFactory(pm.Normal))
class NormalBuilder(FormNode):
    DISPLAY_NAME = "Normal"
    FIELD__mu = MU
    FIELD__sigma = SIGMA

@kwarg_builder(DistFactory(pm.Cauchy))
class CauchyBuilder(FormNode):
    DISPLAY_NAME = "Normal"
    FIELD__alpha = create_value_field(FloatNode, "alpha", 0.0)
    FIELD__beta = create_value_field(FloatNode, "beta", 0.0)


@kwarg_builder(DistFactory(pm.HalfNormal))
class HalfNormalBuilder(FormNode):
    DISPLAY_NAME = "HalfNormal"
    FIELD__sigma = SIGMA
    FIELD__negative = create_value_field(BoolNode,"Negative",False)


@kwarg_builder(DistFactory(pm.Uniform))
class UniformBuilder(FormNode):
    DISPLAY_NAME = "Uniform"
    FIELD__lower = LOWER
    FIELD__upper = UPPER


class ExponentialReverse(FloatNode):
    DISPLAY_NAME = "mu"
    DEFAULT_VALUE = 1.0

    def get_data(self):
        data = super().get_data()
        return 1/data

class ExponentialAlt(AlternatingNode):
    DISPLAY_NAME = "parameter"
    SEL__mu = ExponentialReverse
    SEL__lambda = create_value_field(FloatNode,"lambda",1.0)

@kwarg_builder(DistFactory(pm.Exponential))
class ExponentialBuilder(FormNode):
    DISPLAY_NAME = "Exponent"
    FIELD__lam = ExponentialAlt
    FIELD__negative = create_value_field(BoolNode,"Negative",False)

class TruncatedBuilder(FormNode):
    DISPLAY_NAME = "Truncated"
    FIELD__lower = LowerOpt
    FIELD__upper = UpperOpt

    def get_data(self):
        data = super().get_data()
        if isinstance(data["dist"], ConstantDistCreator):
            return data["dist"]  # Ignore truncation
        else:
            data["dist"] = data["dist"].get_dist()
        return DistCreator(pm.Truncated, **data)


class ConstantBuilder(FloatNode):
    DISPLAY_NAME = "Const"
    DEFAULT_VALUE = 0.0

    def get_data(self):
        data = super().get_data()
        return ConstantDistCreator(data)


class DistBuilder(AlternatingNode):
    DISPLAY_NAME = ""
    SEL__const = ConstantBuilder
    SEL__normal = NormalBuilder
    SEL__halfnormal = HalfNormalBuilder
    SEL__uniform = UniformBuilder
    SEL__cauchy = CauchyBuilder
    SEL__exponential = ExponentialBuilder
    SEL__truncated = TruncatedBuilder

TruncatedBuilder.FIELD__dist = DistBuilder
