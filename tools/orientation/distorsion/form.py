from vtl_common.common_GUI.tk_forms_assist import FormNode, ArrayNode, OptionNode
from tools.new_reconstructor.models.form_prototypes.distribution_parameter import DistBuilder, DistParameterCreator

class CoefficientArray(ArrayNode):
    ITEM_TYPE = DistBuilder

class RadialCoefficientArray(CoefficientArray):
    DISPLAY_NAME = "R"
