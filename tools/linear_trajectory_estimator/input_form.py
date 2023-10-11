import datetime

from vtl_common.common_GUI.tk_forms import TkDictForm
from vtl_common.localized_GUI.tk_forms import SaveableTkDictForm
from vtl_common.common_GUI.tk_forms_assist import FormNode, FloatNode, StringNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field
from tools.new_reconstructor.orientation_pointing_form import OrientedPoint
from vtl_common.localization import get_locale
from vtl_common.workspace_manager import Workspace

TRACK2D_PARAMETERS = Workspace("track_2d_parameters")

DT_DEFAULT = "%Y-%m-%d %H:%M:%S:0"
DEFAULT_DATETIME = datetime.datetime.now().replace(microsecond=0).strftime(DT_DEFAULT)

M_MM = get_locale("measurements.mm")
M_KM = get_locale("measurements.km")
M_S = get_locale("measurements.s")
M_MRAD = get_locale("measurements.mrad")
M_PIX = get_locale("measurements.pix")
M_FR = get_locale("measurements.fr")

class DataForm(FormNode):
    FIELD__tres = create_value_field(FloatNode,get_locale("tools.trajectory_calculator.temporal_resolution"),0.2)
    FIELD__a = create_value_field(FloatNode,f"A [{M_PIX}/{M_FR}^2]",0.0)
    FIELD__u_z = create_value_field(FloatNode,f"u_z [1/{M_FR}]",0.0)
    FIELD__u0=create_value_field(FloatNode,f"U0 [{M_PIX}/{M_FR}]",0.1)
    FIELD__phi0=create_value_field(FloatNode,"phi_0 [°]",0.0)
    FIELD__x0=create_value_field(FloatNode,f"x0 [{M_MM}]",0.0)
    FIELD__y0=create_value_field(FloatNode,f"y0 [{M_MM}]",0.0)
    FIELD__k0 = create_value_field(FloatNode,f"k0 [{M_FR}]",0.0)
    FIELD__k = create_value_field(FloatNode,f"k [{M_FR}]",0.0)
    FIELD__v = create_value_field(FloatNode,f"v [{M_KM}/{M_S}]",35.0)
    FIELD__radiant = create_value_field(OrientedPoint, get_locale("tools.trajectory_calculator.radiant"))
    FIELD__k0_datetime = create_value_field(StringNode, get_locale("tools.trajectory_calculator.datetime"),
                                            DEFAULT_DATETIME)


def make_form(master):
    parser = DataForm()
    conf = parser.get_configuration_root()
    return SaveableTkDictForm(master,conf, file_asker=TRACK2D_PARAMETERS), parser
