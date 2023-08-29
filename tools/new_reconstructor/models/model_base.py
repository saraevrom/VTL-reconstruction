import inspect
from pprint import pprint

import arviz as az
import numpy as np
import pymc as pm
from vtl_common.common_GUI.tk_forms_assist import FormNode
from .form_prototypes import FinalDistributionField


def create_deterministic(varname, value, consts):
    if isinstance(value, float):
        consts[varname] = value
        return value
    else:
        return pm.Deterministic(varname, value)

def pmt_break_mask(pmt: str):
    res = np.full(shape=(16, 16), fill_value=True)
    if "A" in pmt:
        res[:8, 8:] = False
    elif "B" in pmt:
        res[8:, 8:] = False
    elif "C" in pmt:
        res[:8, :8] = False
    elif "D" in pmt:
        res[8:, :8] = False
    return res

class FormPrototype(object):
    @classmethod
    def get_static(cls):
        res = dict(cls.__dict__)
        for base in inspect.getmro(cls):
            if base != cls:
                if hasattr(base, "get_static"):
                    res.update(base.get_static())
        return res

    def generate_subform(self):
        statics = type(self).get_static()

        class ParametersForm(FormNode):
            DISPLAY_NAME = ""

            def get_data(self_internal):
                data = super().get_data()
                #print("DATA", data)
                self.set_params(data)
                return self.get_form_result()

        for key in statics.keys():
            if hasattr(statics[key], "generate"):
                setattr(ParametersForm, "FIELD__" + key, statics[key].generate(key))

        return ParametersForm

    def get_form_result(self):
        return self

    def set_params(self, params_result:dict):
        for k in params_result.keys():
            setattr(self, k, params_result[k])


class ModelWithParameters(object):
    def __init__(self, model, parameters, consts):
        self.model = model
        self.parameters = parameters
        self.pmt = None
        self.idata = None
        self.parent = None
        self.cut_range = None
        self.consts = consts
        print("PARAMETERS SET")
        pprint(parameters)
        print("CONSTANTS")
        pprint(consts)

    def sample(self, *args, **kwargs):
        with self.model:
            self.idata = pm.sample(*args, **kwargs)

    def reconstruction_overlay(self, plotter):
        self.parent.reconstruction_overlay(plotter, self)

    def postprocess(self, axes, actual_x):
        return self.parent.postprocess(axes, self, actual_x)

    def get_estimation(self, key, use_float=True):
        if key in self.consts:
            return self.consts[key]
        elif use_float:
            return float(self.idata.posterior[key].median())
        else:
            return self.idata.posterior[key].median()

    def get_posterior_variables(self):
        return list(self.idata.posterior.keys())

    def ask_parameter(self, parameter_name):
        callable_name = f"ask_{parameter_name}"
        if hasattr(self.parent, callable_name):
            return getattr(self.parent,callable_name)(self)
        return None

    def get_summary(self):
        whole_summary = az.summary(self.idata.posterior, stat_focus="median")
        whole_summary.insert(0, 'parameter', whole_summary.index)
        whole_summary = whole_summary.reset_index(drop=True)
        return self.consts.copy(), whole_summary


class ReconstructionModelWrapper(FormPrototype):
    '''
    Wraps pymc model to create shared functionality.
    Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    '''

    final_distribution = FinalDistributionField() # One common field

    def get_form_result(self):
        return self

    def generate_pymc_model(self, observed, cut_start, cut_end, broken, pmt, reconstructor_main) -> ModelWithParameters:
        raise NotImplementedError("Cannot get model")

    @staticmethod
    def process_mask(broken, pmt):
        mask = pmt_break_mask(pmt)
        return np.logical_or(mask, broken)

    def get_pymc_model(self, observed, cut_start, cut_end, broken, pmt, reconstructor_main):
        b = self.process_mask(broken, pmt)
        res = self.generate_pymc_model(observed, cut_start, cut_end, b, pmt, reconstructor_main)
        res.pmt = pmt
        res.parent = self
        res.cut_range = cut_start, cut_end
        return res

    def reconstruction_overlay(self, plotter, model_params: ModelWithParameters):
        pass

    def postprocess(self, ax, model_params: ModelWithParameters, actual_x):
        pass