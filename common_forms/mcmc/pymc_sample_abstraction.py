import pymc
import pymc as pm
import pytensor.tensor as pt

class PyMCSampleAbstraction(object):
    def sample(self,model:pymc.Model):
        raise NotImplementedError

class JustSample(PyMCSampleAbstraction):
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs
        self.args = args

    def sample(self,model:pymc.Model):
        with model:
            idata = pm.sample(*self.args,**self.kwargs)
        return idata


class VariationalSample(PyMCSampleAbstraction):
    def __init__(self,fit_kwargs, sample_kwargs):
        self.fit_kwargs = fit_kwargs
        self.sample_kwargs = sample_kwargs

    def sample(self,model:pymc.Model):
        with model:
            approx = pm.fit(**self.fit_kwargs)
            idata = approx.sample(**self.sample_kwargs)
        return idata