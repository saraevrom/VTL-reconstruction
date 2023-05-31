import pymc as pm
from vtl_common.common_GUI.tk_forms_assist import ComboNode
from vtl_common.localization import get_locale
from .base import FieldPrototype



def add_kwarg(tgt,key,value):
    if key not in tgt.keys():
        tgt[key] = value

def wrap_normal(*args, **kwargs):
    '''
    accepts same args as pymc.Normal
    :param args:
    :param kwargs:
    :return:
    '''
    sigma0 = pm.HalfNormal('Sigma0', 1.)
    add_kwarg(kwargs, "sigma", sigma0)
    return pm.Normal(*args, **kwargs)


def wrap_student(*args, **kwargs):
    '''
        accepts same args as pymc.StudentT
        :param args:
        :param kwargs:
        :return:
        '''
    sigma0 = pm.HalfNormal('Sigma0', 1.)
    nu = pm.Exponential('nu', 0.2)
    add_kwarg(kwargs, "sigma", sigma0)
    add_kwarg(kwargs, "nu", nu)
    return pm.StudentT(*args, **kwargs)

def wrap_cauchy(*args, **kwargs):
    '''
    accepts same args as pymc.Cauchy
    :param args:
    :param kwargs:
    :return:
    '''
    mu = kwargs["mu"]  # Since mean of cauchy is alpha, we will replace mu with alpha
    beta = pm.HalfNormal("beta", 1.0)
    add_kwarg(kwargs, "alpha", mu)
    del kwargs["mu"]
    add_kwarg(kwargs, "beta", beta)
    pm.Cauchy(*args, **kwargs)

WRAPPERS = {
    "normal": wrap_normal,
    "student": wrap_student,
    "cauchy": wrap_cauchy
}

class FinalDistributionNode(ComboNode):
    DISPLAY_NAME = get_locale("reconstruction.form.distribution")
    VALUES = list(WRAPPERS.keys())
    DEFAULT_VALUE = list(WRAPPERS.keys())[0]
    SELECTION_READONLY = True

    def get_data(self):
        val = super().get_data()
        return WRAPPERS[val]

class FinalDistributionField(FieldPrototype):


    def generate(self, display_name):
        return FinalDistributionNode