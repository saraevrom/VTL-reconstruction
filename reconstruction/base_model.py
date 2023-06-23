import inspect
from vtl_common.common_GUI.tk_forms_assist import FormNode
from .form_prototypes import FinalDistributionField

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


class ReconstructionModelWrapper(FormPrototype):
    '''
    Wraps pymc model to create shared functionality.
    Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    '''

    final_distribution = FinalDistributionField() # One common field

    def get_form_result(self):
        return self, self.get_pymc_model

    def get_pymc_model(self, observed, cut_start, cut_end, broken):
        raise NotImplementedError("Cannot get model")

    def reconstruction_overlay(self, plotter, i_trace, offset):
        x_off, y_off = offset

    def postprocess(self, ax, k_start, k_end, pmt, trace):
        pass