from fixed_rotator.astro_math_z_aligned import Vector3, Quaternion, unixtime_to_era, eci_to_ocef, ocef_to_altaz, detector_plane_to_ocef
from fixed_rotator.astro_math_z_aligned import Vector2, ocef_to_eci, eci_to_ocef
import numpy as np

from vtl_common.common_GUI.tk_forms_assist import FormNode, LabelNode, FloatNode, AlternatingNode
from vtl_common.common_GUI.tk_forms_assist.factory import create_value_field, kwarg_builder

class Projection(object):
    def polar_projection(self, ocef:Vector3):
        raise NotImplementedError

    def radial_range(self):
        return NotImplementedError

class Equidistant(Projection):
    def polar_projection(self, ocef:Vector3):
        alt, az = ocef_to_altaz(ocef)
        rs = 90 - alt * 180 / np.pi
        return az, rs

    def radial_range(self):
        return 0.0, 90.0


class Gnomonic(Projection):
    def __init__(self, cutoff):
        self.cutoff = cutoff

    def polar_projection(self, ocef:Vector3):
        alt, az = ocef_to_altaz(ocef)
        theta = np.pi/2 - alt
        rs = np.tan(theta)
        return az, rs

    def radial_range(self):
        return 0., self.cutoff


class Ortographic(Projection):
    def radial_range(self):
        return 0., 1.

    def polar_projection(self, ocef:Vector3):
        alt, az = ocef_to_altaz(ocef)
        theta = np.pi/2 - alt
        rs = np.sin(theta)
        rs[np.where(theta>np.pi/2)] = 2
        return az, rs



@kwarg_builder(Equidistant)
class EquidistantForm(FormNode):
    DISPLAY_NAME = "Equidistant projection"

@kwarg_builder(Gnomonic)
class GnomonicForm(FormNode):
    DISPLAY_NAME = "Gnomonic projection"
    FIELD__cutoff = create_value_field(FloatNode, "Radial cutoff", default_value=1.0)


@kwarg_builder(Ortographic)
class OrtographicForm(FormNode):
    DISPLAY_NAME = "Ortographic projection"


class ProjectionSelector(AlternatingNode):
    DISPLAY_NAME = "Projection"
    SEL__equidistant = EquidistantForm
    SEL__gnomonic = GnomonicForm
    SEL__ortographic = OrtographicForm

