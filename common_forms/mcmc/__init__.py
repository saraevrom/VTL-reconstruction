from vtl_common.common_GUI.tk_forms_assist import AlternatingNode
from vtl_common.localization import get_locale
from .just_sample import JustSampler
from .variational import VariationalSampler


class Sampler(AlternatingNode):
    DISPLAY_NAME = get_locale("reconstruction.sample")
    SEL__normal = JustSampler
    SEL__variational = VariationalSampler
