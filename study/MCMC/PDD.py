#!/usr/bin/env python3

import theano
import warnings
import numpy as np
import theano.tensor as tt

warnings.filterwarnings(action='ignore')
theano.config.compute_test_value = 'warn'

class PDD_LA:
    """ Forward model used for hierarchical MCMC paramter inversion.

    __TO DO__:
        - implenent along "range" accumulation gradient
    """

    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, T_σ, A_mean):
        # These a arguments which all models will need. Any of the model
        # pamaeters which are dependent on the model formulation are passed
        # as kwargs.
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.T_σ    = T_σ
        self.A_mean = A_mean

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
               self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
               np.random.normal(0, self.T_σ, (365,1))

        return Temp

    def __tt_accumulation(self, z, grad_a, T):
        A_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)
        A_snow = tt.maximum((A_days*self.A_mean)*(1+(z-self.ref_z)*grad_a), 0.0)

        return A_snow

    def __tt_component(self, z, f_snow, C, f_r, grad_a):
        """Theano implementation of the forward model which supports shared
           variables as input.

           This is an intermidiate function which returns a list of the three
           components of the PDD model (i.e. Accumulation, Refreezing and melt).

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow [-]
                                1 <= C < ?
            f_r    (float) --> refreezing factor               [-]
            **kwargs
        Outputs:
            [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
                List of theano tensors without numeric values for each component
                of the mass balance model [m i.e yr^-1]
        """

        f_ice  = C*f_snow

        # temperature and PDDs calc
        T      = self._air_temp(z)
        PDDs   = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)

        # accumulation calc
        A_snow = self.__tt_accumulation(z, grad_a, T)

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
        return A_snow * (1/910), R * (1/910), M_melt * (1/910), f_m, r_s2m

    def tt_forward(self, z, f_snow, C, f_r, grad_a):
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
        A_snow_tt, R_tt, M_melt_tt, f_m, r_s2m = self.__tt_component(z, f_snow, C, f_r, grad_a)

        # Return the net balance
        return A_snow_tt + R_tt - M_melt_tt

    def __compile_forward(self):
        """Function to compile the theano implementation of the forward model.
        """
        z      = tt.vector('z')
        f_snow = tt.dscalar('f_s')
        C      = tt.dscalar('C')
        f_r    = tt.dscalar('f_r')
        grad_A = tt.dscalar('grad_A')

        compiled = theano.function([z, f_snow, C, f_r, grad_A],
                            self.__tt_component(z, f_snow, C, f_r, grad_A))

        return compiled

    def eval_forward(self, z, f_snow, C, f_r, grad_a, parts=True):
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
        A_snow, R, M_melt, f_m, r_s2m = self.__compiled(z, f_snow, C, f_r, grad_a)
        # Calculate the net balance from components
        MB = A_snow + R - M_melt

        if parts:
            return MB, A_snow, R, M_melt, f_m, r_s2m
        else:
            return MB

class PDD_PWA:
    """ Forward model used for hierarchical MCMC paramter inversion.

    __TO DO__:
        - implenent along "range" accumulation gradient
    """

    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, T_σ, z_ELA):
        # These a arguments which all models will need. Any of the model
        # pamaeters which are dependent on the model formulation are passed
        # as kwargs.
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.T_σ    = T_σ
        self.z_ELA  = z_ELA

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
               self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
               np.random.normal(0, self.T_σ, (365,1))

        return Temp

    def __tt_accumulation(self, z, T, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max):

        accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)

        A = tt.switch(tt.lt(self.z_ELA,z) & tt.lt(z, z_max), ΔPΔz_2*z + P0-ΔPΔz_2*self.z_ELA, z)
        A = tt.switch(tt.le(z, self.z_ELA), ΔPΔz_1*z + P0-ΔPΔz_1*self.z_ELA, A)
        A = tt.switch(tt.le(z_max, z), P_max, A)

        return A*accum_days*910

    def __tt_component(self, z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max):
        """Theano implementation of the forward model which supports shared
           variables as input.

           This is an intermidiate function which returns a list of the three
           components of the PDD model (i.e. Accumulation, Refreezing and melt).

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
            C      (float) --> factor relating f_ice to f_snow [-]
                                1 <= C < ?
            f_r    (float) --> refreezing factor               [-]
            **kwargs
        Outputs:
            [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
                List of theano tensors without numeric values for each component
                of the mass balance model [m i.e yr^-1]
        """

        f_ice  = C*f_snow

        # PDDs and accumulation calc
        T      = self._air_temp(z)
        PDDs   = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)

        A_snow = self.__tt_accumulation(z, T, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max)

        # calculate local surface melt assuming f_m = f_snow
        melt_local = PDDs * f_snow

        # calculate refreezing
        R = tt.minimum(f_r*A_snow, melt_local)

        #place folder variable
        v = tt.switch(tt.eq(melt_local, 0.0), A_snow, melt_local)
        r_s2m = A_snow/v

        f_m = tt.switch(tt.ge(r_s2m, 1.), f_snow,
                        f_ice - (f_ice - f_snow)*r_s2m)

        # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
        M_melt = f_m*PDDs

        # Return individual components of the mass balance model in [m i.e. / y]
        return A_snow * (1/910), R * (1/910), M_melt * (1/910), f_m, r_s2m

    def tt_forward(self, z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max):
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
        A_snow_tt, R_tt, M_melt_tt, f_m, r_s2m = self.__tt_component(
                                                 z, f_snow, C, f_r, z_max, P0,
                                                 ΔPΔz_1, ΔPΔz_2, P_max)

        # Return the net balance
        return A_snow_tt + R_tt - M_melt_tt

    def __compile_forward(self):
        """Function to compile the theano implementation of the forward model.
        """
        z      = tt.vector('z')
        f_snow = tt.dscalar('f_s')
        C      = tt.dscalar('C')
        f_r    = tt.dscalar('f_r')
        z_max  = tt.dscalar('z_max')
        P0     = tt.dscalar('P0')
        ΔPΔz_1 = tt.dscalar('ΔPΔz_1')
        ΔPΔz_2 = tt.dscalar('ΔPΔz_2')
        P_max  = tt.dscalar('P_max')

        compiled = theano.function([z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max],
                            self.__tt_component(z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max))

        return compiled

    def eval_forward(self, z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max, parts=True):
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
        A_snow, R, M_melt, f_m, r_s2m = self.__compiled(z, f_snow, C, f_r, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max)
        # Calculate the net balance from components
        MB = A_snow + R - M_melt

        if parts:
            return MB, A_snow, R, M_melt, f_m, r_s2m
        else:
            return MB

class Accumulation_default:

    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, T_σ, A_mean):
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.T_σ    = T_σ
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
               self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
               np.random.normal(0, self.T_σ, (365,1))

        return Temp

    def tt_forward(self, z, ΔPΔz):
        T     = self._air_temp(z) #+ T_bias

        PDDs = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)

        accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)

        # calculate snow accumulation
        A_snow = tt.maximum((accum_days*self.A_mean)* (1+(z-self.ref_z)*ΔPΔz), 0.0)
        return A_snow

    def __compile_forward(self):
        """Function to compile the theano implementation of the forward model.
        """
        z      = tt.vector('z')
        ΔPΔz   = tt.dscalar('ΔPΔz')
        # A_mean = tt.dscalar('A_mean')

        compiled = theano.function([z, ΔPΔz], self.tt_forward(z, ΔPΔz))

        return compiled

    def eval_forward(self, z, ΔPΔz):
        """Compiled version of the forward model which takes numeric inputs as
           arguments and output numeric arguments.

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            ΔPΔz   (float) --> alititudinal precip factor      [m^-1]

        Outputs:
            A_snow (ndarray) --> Snow accumulation (optional)   [m i.e. / y]

        """

        # Evaluate individual components of MB model
        A_snow = self.__compiled(z, ΔPΔz)

        return A_snow

class Accumulation_piecewise:


    def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, T_σ, z_ELA):
        self.α      = α       # anual air temp. amplitude    [K]
        self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
        self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
        self.T_p    = T_p     # DOY of annual temp peak      [DOY]
        self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
        self.T_m    = T_m     # Melting temp. threshold      [K]
        self.T_rs   = T_rs    # T_
        self.T_σ    = T_σ
        self.z_ELA = z_ELA   # ELA elevation (known a prior only in the test case)

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
               self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
               np.random.normal(0, self.T_σ, (365,1))

        return Temp

    def tt_forward(self, z, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max):

        T = self._air_temp(z) #+ T_bias
        accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)

        A = tt.switch(tt.lt(self.z_ELA,z) & tt.lt(z, z_max), ΔPΔz_2*z + P0-ΔPΔz_2*self.z_ELA, z)
        A = tt.switch(tt.le(z, self.z_ELA), ΔPΔz_1*z + P0-ΔPΔz_1*self.z_ELA, A)
        A = tt.switch(tt.le(z_max, z), P_max, A)

        return A*accum_days*910

    def __compile_forward(self):
        """Function to compile the theano implementation of the forward model.
        """
        z      = tt.vector('z')
        z_max  = tt.dscalar('z_max')
        P0     = tt.dscalar('P0')
        ΔPΔz_1 = tt.dscalar('ΔPΔz_1')
        ΔPΔz_2 = tt.dscalar('ΔPΔz_2')
        P_max  = tt.dscalar('P_max')

        compiled = theano.function([z, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max],
                   self.tt_forward(z, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max))

        return compiled

    def eval_forward(self, z, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max):
        """Compiled version of the forward model which takes numeric inputs as
           arguments and output numeric arguments.

        Inputs:
            z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
            ΔPΔz   (float) --> alititudinal precip factor      [m^-1]

        Outputs:
            A_snow (ndarray) --> Snow accumulation (optional)   [m i.e. / y]

        """

        # Evaluate individual components of MB model
        A_snow = self.__compiled(z, z_max, P0, ΔPΔz_1, ΔPΔz_2, P_max)

        return A_snow

# class PDD_piecewise(PDD_melt_model):
#     def __tt_component(self, z, f_snow, C, grad_a, f_r):
#         """Theano implementation of the forward model which supports shared
#            variables as input.
#
#            This is an intermidiate function which returns a list of the three
#            components of the PDD model (i.e. Accumulation, Refreezing and melt).
#
#         Inputs:
#             z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
#             f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
#             C      (float) --> factor relating f_ice to f_snow [-]
#                                 1 <= C < ?
#             grad_a (float) --> alititudinal precip factor      [m^-1]
#             f_r    (float) --> refreezing factor               [-]
#
#         Outputs:
#             [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
#                 List of theano tensors without numeric values for each component
#                 of the mass balance model [m i.e yr^-1]
#         """
#         f_ice = C*f_snow
#         T     = self._air_temp(z) #+ T_bias
#
#         PDDs = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)
#
#         accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)
#
#         # # Rounce et al. 2020 / Huss and Hock 2014
#         # tmp = tt.switch(tt.lt(self.T_rs-1, T) & tt.gt(T, self.T_rs+1),
#         #                 0.5+(T-self.T_rs)/2, T)
#         # tmp = tt.switch(tt.le(tmp, self.T_rs - 1), 1 , tmp)
#         # tmp = tt.switch(tt.ge(tmp, self.T_rs + 1), 0 , tmp)
#         # accum_days = (tmp / 365.).sum(axis=0)
#
#         # calculate snow accumulation
#         A_snow = tt.maximum((accum_days*self.A_mean)* (1+(z-self.ref_z)*grad_a),
#                             0.0)
#
#
#         # Rounce et al. 2020 / Huss and Hock 2014
#         tmp = tt.switch(tt.lt(self.T_rs-1, T) & tt.gt(T, self.T_rs+1),
#                         0.5+(T-self.T_rs)/2, T)
#         tmp = tt.switch(tt.le(tmp, self.T_rs - 1), 1 , tmp)
#         tmp = tt.switch(tt.ge(tmp, self.T_rs + 1), 0 , tmp)
#         accum_days = (tmp / 365.).sum(axis=0)
#
#         # calculate local surface melt assuming f_m = f_snow
#         melt_local = PDDs * f_snow
#
#         # calculate refreezing
#         R = tt.minimum(f_r*A_snow, melt_local)
#
#
#         r_s2m = tt.switch(tt.eq(melt_local, 0.0), 1.0, A_snow/melt_local)
#
#         f_m = tt.switch(tt.ge(r_s2m, 1.), f_snow,
#                         f_ice - (f_ice - f_snow)*r_s2m)
#
#         # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
#         M_melt = f_m*PDDs
#
#         # Return individual components of the mass balance model in [m i.e. / y]
#         return [A_snow * (1/910), R * (1/910), M_melt * (1/910)]
#
# class PDD_enhanced:
#     """ Forward model used for hierarchical MCMC paramter inversion.
#
#     __TO DO__:
#         - implenent along "range" accumulation gradient
#     """
#
#     def __init__(self, α, T_ma, ΔTΔz, T_p, ref_z, T_m, T_rs, A_mean):
#         self.α      = α       # anual air temp. amplitude    [K]
#         self.T_ma   = T_ma    # Mean annual air temp @ ref_z [K]
#         self.ΔTΔz   = ΔTΔz    # air temp lapse rate          [K m^-1]
#         self.T_p    = T_p     # DOY of annual temp peak      [DOY]
#         self.ref_z  = ref_z   # reference surface elevation  [m a.s.l.]
#         self.T_m    = T_m     # Melting temp. threshold      [K]
#         self.T_rs   = T_rs    # T_
#         self.A_mean = A_mean  # Mean annual accum.   @ ref_z [kg m^-2 yr^-1]
#
#         # Compile theano forward model
#         self.__compiled = self.__compile_forward()
#
#     def _air_temp(self, z):
#         """"Evaluate the surface air temperature at a given evlevation for a given
#             day of the year
#
#         Inputs:
#             z   (float or Nx1 ndarry) ---> Nodal surface elevation [m a.s.l.]
#
#         Outputs:
#             T   (floar or Nx365 ndarray) ---> Nodal surface elevation for each
#                                               day of the yeat      [C]
#         """
#         doy  = np.arange(1,366)[:, np.newaxis]
#         Temp = self.α*np.cos( 2*np.pi*(doy-self.T_p)/365 ) + \
#                self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
#                np.random.normal(0, 1.5, (365, 1000)) #* T_bias
#                # theano_random.normal(avg=0, std=1, size=z.shape,
#                #                      dtype=theano.config.floatX)
#
#         # print(tt.shape(z))
#         # print(type(z))
#
#         return Temp
#
#     def __tt_component(self, z, f_snow, C, grad_a, f_r, P_bias):
#         """Theano implementation of the forward model which supports shared
#            variables as input.
#
#            This is an intermidiate function which returns a list of the three
#            components of the PDD model (i.e. Accumulation, Refreezing and melt).
#
#         Inputs:
#             z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
#             f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
#             C      (float) --> factor relating f_ice to f_snow [-]
#                                 1 <= C < ?
#             grad_a (float) --> alititudinal precip factor      [m^-1]
#             f_r    (float) --> refreezing factor               [-]
#
#         Outputs:
#             [A_snow, R, M_melt]  ([theano.tt, theano.tt, theano.tt]) -->
#                 List of theano tensors without numeric values for each component
#                 of the mass balance model [m i.e yr^-1]
#         """
#         f_ice = C*f_snow
#
#
#         T   = self._air_temp(z)
#         DOY = tt.as_tensor(np.arange(1,366, dtype='float32')[:, np.newaxis], 'DOY')
#
#         #https://github.com/vitruvianscience/OpenDeep/blob/master/opendeep/utils/noise.py
#         # T = self.α*tt.cos(2*np.pi*(DOY-self.T_p)/365 ) + \
#         #     self.ΔTΔz*(self.ref_z-z)+self.T_ma + \
#         #     theano_random.normal(avg=0, std=0.5, size=z.shape, dtype=theano.config.floatX)
#
#         PDDs = tt.switch(tt.gt(T, self.T_m), T, 0.0).sum(axis=0)
#
#         accum_days = tt.switch(tt.lt(T, self.T_rs), 1/365., 0.0).sum(axis=0)
#
#         # # Rounce et al. 2020 / Huss and Hock 2014
#         # tmp = tt.switch(tt.lt(self.T_rs-1, T) & tt.gt(T, self.T_rs+1),
#         #                 0.5+(T-self.T_rs)/2, T)
#         # tmp = tt.switch(tt.le(tmp, self.T_rs - 1), 1 , tmp)
#         # tmp = tt.switch(tt.ge(tmp, self.T_rs + 1), 0 , tmp)
#         # accum_days = (tmp / 365.).sum(axis=0)
#
#         # calculate snow accumulation
#         # A_snow = tt.maximum((accum_days*self.A_mean)*(1+(z-self.ref_z)*grad_a),
#         #                 0.0)
#         A_snow = tt.maximum((accum_days*self.A_mean)*P_bias*(1+(z-self.ref_z)*grad_a),
#                             0.0)
#
#         # calculate local surface melt assuming f_m = f_snow
#         melt_local = PDDs * f_snow
#
#         # calculate refreezing
#         R = tt.minimum(f_r*A_snow, melt_local)
#
#
#         r_s2m = tt.switch(tt.eq(melt_local, 0.0), 1.0, A_snow/melt_local)
#
#         f_m = tt.switch(tt.ge(r_s2m, 1.), f_snow,
#                         f_ice - (f_ice - f_snow)*r_s2m)
#
#         # calculate surface melt [kg m^{-2} yr^{-1}] with f_m
#         M_melt = f_m*PDDs
#
#         # Return individual components of the mass balance model in [m i.e. / y]
#         return [A_snow * (1/910), R * (1/910), M_melt * (1/910)]
#
#     def tt_forward(self, z, f_snow, C, grad_a, f_r, P_bias):
#         """Wrapper for the individual component model. Computes the sum of the
#            individual components, resulting in the net balance (i.e. target value
#            in tunning procedure).
#
#         Inputs:
#             z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
#             f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
#             C      (float) --> factor relating f_ice to f_snow [-]
#                                 1 <= C < ?
#             grad_a (float) --> alititudinal precip factor      [m^-1]
#             f_r    (float) --> refreezing factor               [-]
#
#         Outputs:
#             MB (theano.tt) -->  theano tensor without numeric values [m i.e. / y]
#
#         """
#         # DOY = tt.drow('DOY')
#
#         A_snow_tt, R_tt, M_melt_tt = self.__tt_component(z, f_snow, C, grad_a, f_r, P_bias)
#
#         # Return the net balance
#         return A_snow_tt + R_tt - M_melt_tt
#
#     def __compile_forward(self):
#         """Function to compile the theano implementation of the forward model.
#         """
#         z      = tt.vector('z')
#         f_snow = tt.dscalar('f_s')
#         C      = tt.dscalar('C')
#         grad_a = tt.dscalar('grad_a')
#         f_r    = tt.dscalar('f_r')
#         # T_bias = tt.dscalar('T_bias')
#         P_bias = tt.dscalar('P_bias')
#         # DOY    = tt.drow('DOY')
#
#         # print(tt.shape(z))
#
#         compiled = theano.function([z, f_snow, C, grad_a, f_r, P_bias],
#                             self.__tt_component(z, f_snow, C, grad_a, f_r, P_bias))
#
#         return compiled
#
#     def eval_forward(self, z, f_snow, C, grad_a, f_r, P_bias, parts=True):
#         """Compiled version of the forward model which takes numeric inputs as
#            arguments and output numeric arguments.
#
#         Inputs:
#             z    (ndarray) --> Nx1 array of elevations         [m a.s.l.]
#             f_snow (float) --> degree-day factor for snow      [kg m^-2 yr^-1 K^-1]
#             C      (float) --> factor relating f_ice to f_snow
#                                 1 <= C < ?                     [-]
#             grad_a (float) --> alititudinal precip factor      [m^-1]
#             f_r    (float) --> refreezing factor               [-]
#             parts  (bool)  --> whether to return individual
#                                componetns of MB model.
#
#         Outputs:
#             MB     (ndarray) --> Nodal surface mass balance     [m i.e. / y]
#
#         Optional Outputs (Set by `parts` kwarg):
#             A_snow (ndarray) --> Snow accumulation (optional)   [m i.e. / y]
#             R      (ndarray) --> Refreezing                     [m i.e. / y]
#             M_melt (ndarray) --> Melting                        [m i.e. / y]
#         """
#
#         # DOY = np.arange(1,366)[:, np.newaxis]
#
#         # Evaluate individual components of MB model
#         A_snow, R, M_melt = self.__compiled(z, f_snow, C, grad_a, f_r, P_bias)
#         # Calculate the net balance from components
#         MB = A_snow + R - M_melt
#
#         if parts:
#             return MB, A_snow, R, M_melt
#         else:
#             return MB
