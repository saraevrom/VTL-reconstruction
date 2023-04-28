import inspect
from vtl_common.common_GUI.tk_forms_assist import FormNode
from .form_prototypes import FinalDistributionField


class ReconstructionModelWrapper(object):
    '''
    Wraps pymc model to create shared functionality.
    Static field inheriting FieldPrototype class will be passed to form and transformed automatically
    '''

    final_distribution = FinalDistributionField() # One common field

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
                return self, lambda x: self.get_model(x, data)

        for key in statics.keys():
            if hasattr(statics[key], "generate"):
                setattr(ParametersForm, "FIELD__"+key, statics[key].generate(key))

        return ParametersForm

    def get_model(self, observed, params_result:dict):
        for k in params_result.keys():
            setattr(self, k, params_result[k])
        return self.get_pymc_model(observed)

    def get_pymc_model(self, observed):
        raise NotImplementedError("Cannot get model")

    def reconstruction_overlay(self, plotter, i_trace, offset):
        x_off, y_off = offset