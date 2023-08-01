from vtl_common.common_GUI.tk_forms_assist import IntNode, FloatNode, BoolNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from .base import FieldPrototype
from .distribution_parameter import DistBuilder


class PassthroughField(FieldPrototype):
    def __init__(self, required_class_creator):
        self.req_class_creator = required_class_creator

    def generate(self, display_name):
        return self.req_class_creator()

class BasicField(FieldPrototype):
    def __init__(self, req_type, default_value):
        self.req_type = req_type
        self.default_value = default_value

    def generate(self, display_name):
        return create_value_field(self.req_type, display_name, self.default_value)


class FloatField(BasicField):
    def __init__(self, default_value):
        super().__init__(FloatNode, default_value)


class IntField(BasicField):
    def __init__(self, default_value):
        super().__init__(IntNode, default_value)


class BoolField(BasicField):
    def __init__(self, default_value):
        super().__init__(BoolNode, default_value)


class DistributionField(BasicField):
    def __init__(self, distribution, *args, **params):
        if args:
            default_value = {
                "selection_type": distribution,
                "value": args[0]
            }
        else:
            default_value = {
                "selection_type": distribution,
                "value": params
            }
        super().__init__(DistBuilder, default_value)