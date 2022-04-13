#!/usr/bin/env python3

import json
import theano
import arviz as az
import pymc3 as pm
import numpy as np
import theano.tensor as tt


# inspo:# https://github.com/YuRiTan/prediction-uncertainty/blob/master/src/model/abstract_model.py
# https://github.com/parsing-science/pymc3_models
class AbstractModel:
    def fit(self, *args, **kwargs):
        raise NotImplementedError

    def sample_posterior_predictive(self, *args, **kwargs):
        raise NotImplementedError

    def save(self, *args, **kwargs):
        raise NotImplementedError

    def load(self, *args, **kwargs):
        raise NotImplementedError

class MultivariateInversion(AbstractModel):
    """
    https://arviz-devs.github.io/arviz/schema/PyMC3_schema_example.html
    """

    def __init__(self, forward):

        self.model  = None
        self.trace  = None # inference data
        self.priors = {}
        self.priors_dict = {}


        # forward model instance of PositiveDegreeDayModel
        self.forward = forward

        # Dummy shared variables
        shared_z   = theano.shared(np.zeros((1)))
        shared_RAM = theano.shared(np.zeros((1, 3)))
        self.shared_vars = {'z'  : shared_z,
                            'RAM': shared_RAM}

    def create_model(self, priors, fit_type):
        # save the dictionary of priors to attr. (usefull for plotting)
        self.priors_dict = priors

        # initialize the pymc3 model class
        self.model = pm.Model()
        # create priors for physical model params
        self.__physical_priors(priors)

        # Actual forward model prediction
        R, A, M = self.forward.forward(self.shared_vars['z'], **self.priors)
        mu_pred = tt.transpose(tt.stack([R, A, M]))


        if fit_type == 'MvNormal_correlated':
            # create priors of hyperparameters (observation variance)
            with self.model:
                sd_dist = pm.Exponential.dist(1.0, shape=3)

                chol, corr, stds = pm.LKJCholeskyCov('chol_cov',
                                                     n=3,
                                                     eta=2,
                                                     sd_dist=sd_dist,
                                                     compute_corr=True)

            # create the (multivariate) likelihood model
            with self.model:
                pm.MvNormal('obs', mu=mu_pred, chol=chol, observed=self.shared_vars['RAM'])

            return self.model

        elif fit_type == 'Normal_multidim':

            with self.model:
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # ----> Hyperparameters (likelihood related priors)
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                R_sigma = pm.HalfCauchy("R_sigma", 0.05 )
                A_sigma = pm.HalfCauchy("A_sigma", 0.75 )
                M_sigma = pm.HalfCauchy("M_sigma", 0.5 )

                sigmas  = tt.transpose(tt.stack([R_sigma, A_sigma, M_sigma]))

                # create the likelihood model
                pm.Normal("obs", mu=mu_pred, sigma=sigmas, observed=self.shared_vars['RAM'])

            return self.model


    def _set_shared_vars(self, values):
        for k, v in values.items():
            self.shared_vars[k].set_value(v)

    def __physical_priors(self, priors):
        """
        """
        # parameters that must always be passed.
        # Don't need to be priors, constants  are fine
        params = ['f_s', 'f_r', 'grad_A', 'A_mean']

        # Check that required params were passed
        not_passed = []
        for param in params:
            if param not in priors:
                not_passed.append(param)

        if not_passed:
            raise NameError(f"{not_passed}: Not passed as kwarg.")

        # Now check about "C" OR "f_i"
        if any(param in priors for param in ["C", "f_i"]):
            if all(param in priors for param in ["C", "f_i"]):
                raise NameError("Either C or f_i should be passed, not both.")


        with self.model:
            for key, dict in priors.items():
                self.priors[key] = getattr(pm, dict["type"])(key, **dict["kwargs"])

    # def __hyperpriors(self):
    #     """
    #     """
    #     with self.model:
    #         sd_dist = pm.Exponential.dist(1.0,
    #                                       shape=self.shared_vars['RAM'].shape[1])
    #
    #         chol, corr, stds = pm.LKJCholeskyCov('chol_cov',
    #                                              n=3,
    #                                              eta=2,
    #                                              sd_dist=sd_dist,
    #                                              compute_corr=True)
    #
    # def __likelihood(self):
    #     """
    #     """
    #     pm.MvNormal('vals', mu=mu_pred, chol=chol, observed=mu_obs)


    def fit(self, z, RAM, priors, fit_type, **fit_kwargs):
        """ SKlearn `fit()`-like method to sample with input features `x`
        and observations `y` using `pm.sample`
        """

        self._set_shared_vars({'z': z, 'RAM': RAM})
        self.model = self.create_model(priors, fit_type)
        self.trace = pm.sample(model=self.model, **fit_kwargs)

    def sample_posterior_predictive(self, **inference_kwargs):
        """ SKlearn `predict()`-like method to sample from the posterior
        predictive using `pm.sample_ppc`.
        """
        if self.trace is None:
            raise AttributeError("Please fit the model before predicting")
        if 'samples' not in inference_kwargs:
            inference_kwargs['samples'] = len(self.trace) * self.trace.nchains

        # actually sample the posterior
        with self.model:
            posterior = pm.sample_posterior_predictive(trace=self.trace, **inference_kwargs)

        # return shape: batch_size x sample_size
        posterior = posterior['obs']

    def save(self, filepath, **save_kwargs):
        """ Saves trace of the PyMC3 model. """
        pm.save_trace(self.trace, directory=filepath, **save_kwargs)

    def load(self, filepath):
        """ recreates the model, and loads the trace of the PyMC3 model. """
        self.model = self._create_model()
        self.trace = pm.load_trace(filepath, model=self.model)
        return self

class HierarchicalInversion(AbstractModel):
    pass
