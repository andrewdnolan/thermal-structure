#!/usr/bin/env python3

import theano
import numpy as np
import theano.tensor as tt

class PDD_melt_model:
    """ Forward model used for hierarchical MCMC paramter inversion.

    __TO DO__:
        - implenent along "range" accumulation gradient
    """

    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, A_mean):
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.A_mean = A_mean  # Mean annual accum.   @ ref_z [kg m^-2 yr^-1]

        # Compile theano forward model
        self.__compiled = self.__compile_forward()

    def _air_temp(self, z):
        """"Evaluate the surface air temperature at a given evlevation for a given
            day of the year

        Inputs:
            z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]

        Outputs:
            T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
                                              day of the yeat      [C]
        """
        doy  = np.arange(1,366)[:, np.newaxis]
        Temp = self.α*np.cos( 2*np.pi*(doy-self.T_p)/365 ) + \
               self.ΔTΔz*(self.ref_z-z)+self.T_ma
        return Temp

    def __tt_component(self, z, f_snow, C, grad_a, f_r):
        """Theano implementation of the forward model which supports shared
           variables as input.

           This is an intermidiate function which returns a list of the three
           components of the PDD model (i.e. Accumulation, Refreezing and melt).

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow [-]
                                1 <= C < ?
            grad_a (float) --> alititudinal precip factor      [m^-1]
            f_r    (float) --> refreezing factor               [-]

        Outputs:
            [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
                List of theano tensors without numeric values for each component
                of the mass balance model [m i.e yr^-1]
        """
        f_ice = C*f_snow
        T     = self._air_temp(z) #+ T_bias

        PDDs = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)

        accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)

        # calculate snow accumulation
        A_snow = tt.maximum((accum_days*self.A_mean)*
                          (1+(z-self.ref_z)*grad_a), 0.0)

        # calculate local surface melt assuming f_m = f_snow
        melt_local = PDDs * f_snow

        # calculate refreezing
        R = tt.minimum(f_r*A_snow, melt_local)


        r_s2m = tt.switch(tt.eq(melt_local, 0.0), 1.0, A_snow/melt_local)

        f_m = tt.switch(tt.ge(r_s2m, 1.), f_snow,
                        f_ice - (f_ice - f_snow)*r_s2m)

        # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
        M_melt = f_m*PDDs

        # Return individual components of the mass balance model in [m i.e. / y]
        return [A_snow * (1/910), R * (1/910), M_melt * (1/910)]

    def tt_forward(self, z, f_snow, C, grad_a, f_r):
        """Wrapper for the individual component model. Computes the sum of the
           individual components, resulting in the net balance (i.e. target value
           in tunning procedure).

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow [-]
                                1 <= C < ?
            grad_a (float) --> alititudinal precip factor      [m^-1]
            f_r    (float) --> refreezing factor               [-]

        Outputs:
            MB (theano.tt) -->  theano tensor without numeric values [m i.e. / y]

        """
        A_snow_tt, R_tt, M_melt_tt = self.__tt_component(z, f_snow, C, grad_a, f_r)

        # Return the net balance
        return A_snow_tt + R_tt - M_melt_tt

    def __compile_forward(self):
        """Function to compile the theano implementation of the forward model.
        """
        z      = tt.vector('z')
        f_snow = tt.dscalar('T_m')
        C      = tt.dscalar('C')
        grad_a = tt.dscalar('grad_a')
        f_r    = tt.dscalar('f_r')

        compiled = theano.function([z, f_snow, C, grad_a, f_r],
                            self.__tt_component(z, f_snow, C, grad_a, f_r))

        return compiled

    def eval_forward(self, z, f_snow, C, grad_a, f_r, parts=True):
        """Compiled version of the forward model which takes numeric inputs as
           arguments and output numeric arguments.

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow
                                1 <= C < ?                     [-]
            grad_a (float) --> alititudinal precip factor      [m^-1]
            f_r    (float) --> refreezing factor               [-]
            parts  (bool)  --> whether to return individual
                               componetns of MB model.

        Outputs:
            MB     (ndarray) --> Nodal surface mass balance     [m i.e. / y]
        Optional Outputs (Set by `parts` kwarg):
            A_snow (ndarray) --> Snow accumulation (optional)   [m i.e. / y]
            R      (ndarray) --> Refreezing                     [m i.e. / y]
            M_melt (ndarray) --> Melting                        [m i.e. / y]
        """

        # Evaluate individual components of MB model
        A_snow, R, M_melt = self.__compiled(z, f_snow, C, grad_a, f_r)
        # Calculate the net balance from components
        MB = A_snow + R - M_melt

        if parts:
            return MB, A_snow, R, M_melt
        else:
            return MB
