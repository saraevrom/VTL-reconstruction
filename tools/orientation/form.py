from vtl_common.common_GUI.tk_forms_assist import FormNode, BoolNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from vtl_common.localization import get_locale

class OrientationForm(FormNode):
    FIELD__use_ff = create_value_field(BoolNode, get_locale("orientation.form.use_ff"), True)