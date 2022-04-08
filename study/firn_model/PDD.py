#!/usr/bin/env python3

import numpy as np

SEED = 1234567

class PositiveDegreeDayModel:
    def __init__(self, α, ref_z, T_ma, ΔTΔz, T_p, T_m, T_rs, T_σ):
        """Theano implementation of the (P)ositive (D)egree (D)ay model.

        Inputs:
            α     (float)   --> anual air temp. amplitude    []
            ref_z (float)   --> reference surface elevation  [m a.s.l.]
            T_ma  (float)   --> Mean annual air temp @ ref_z [C]
            ΔTΔz  (float)   --> air temp lapse rate          [K m^-1]
            T_p   (int)     --> DOY of annual temp peak      [DOY]
            T_m   (float)   --> Melting temp. threshold      [C]
            T_rs  (float)   --> rain to snow threshold       [C]
            T_σ   (ndarray) --> std. dev. of air temp as     [C]
                                fuction of DOY

        """
        self.α      = α       # anual air temp. amplitude    []
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [C]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [C]
        self.T_rs   = T_rs    # rain to snow threshold       [C]
        self.T_σ    = T_σ     # std. dev. of air temp        [C]

    def calc_air_temp(self, z):
        """"Evaluate the surface air temperature at a given elev. for a given
            day of the year

        Inputs:
            z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]

        Outputs:
            T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
                                              day of the yeat      [C]
        """
        # Need to seed the Generator every function call to ensure
        # all MCMC steps have the same temperature forcing
        random = np.random.default_rng(SEED)

        # Use array broadcasting to calc temp as function of elevation and doy
        doy  = np.arange(1,366)[:, np.newaxis]
        Temp = self.α * np.cos( 2*np.pi*(doy-self.T_p)/365 ) + \
               self.ΔTΔz * (self.ref_z-z) + self.T_ma + \
               random.normal(0, self.T_σ, (365,1))

        return Temp


    def forward(self, z, **kwargs):
        # parameters that must always be passed.
        # Don't need to be priors, constants  are fine
        params = ['f_s', 'f_r', 'grad_A', 'A_mean']

        # Check that required params were passed
        not_passed = []
        for param in params:
            if param not in kwargs:
                not_passed.append(param)
        if not_passed:
            raise NameError(f"{not_passed}: Not passed as kwarg.")

        # Unpack the kwargs (not pretty but works)
        f_s    = kwargs['f_s']
        f_r    = kwargs['f_r']
        grad_A = kwargs['grad_A']
        A_mean = kwargs['A_mean']

        # Now check about "C" OR "f_i"
        if any(param in kwargs for param in ["C", "f_i"]):
            if all(param in kwargs for param in ["C", "f_i"]):
                raise NameError("Either C or f_i should be passed, not both.")

            # unpack or create f_i depending on what was passed
            if "C" in kwargs:
                f_i = kwargs["C"]*f_s
            elif "f_i" in kwargs:
                f_i = kwargs["f_i"]

        # Temp. and PDDs calculation
        Temp = self.calc_air_temp(z)

        PDDs       = np.where(Temp>self.T_m, Temp-self.T_m, 0).sum(axis=0)
        accum_days = np.where(Temp<self.T_rs,     1/365., 0).sum(axis=0)

        # calculate snow accumulation
        A_snow=np.maximum((accum_days*A_mean)*
                          (1+(z-self.ref_z)*grad_A), 0.0)

        # calculate local surface melt assuming f_m = f_snow
        melt_local = PDDs * f_s
        # calculate refreezing
        R = np.minimum(f_r*A_snow, melt_local)
        # compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
        with np.errstate(divide='ignore', invalid='ignore'):
            r_s2m = np.where(melt_local==0.0, 1.0, A_snow / melt_local)

        # Compute nodal specific degree day factor
        f_m = np.where(r_s2m >= 1, f_s, f_i - (f_i - f_s) * r_s2m )
        # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
        M_melt = f_m*PDDs

        return A_snow, M_melt
