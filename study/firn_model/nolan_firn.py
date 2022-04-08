#!/usr/bin/env python3

import numpy as np
import xarray as xr
import scipy.linalg as LA
import matplotlib.pyplot as plt

from PDD import PositiveDegreeDayModel
SEED = 1234567

def tridiag(a,b,c,r):

    j   = a.size
    u   = np.zeros(j)
    gam = np.zeros(j)

    if b[0] == 0:
        raise Exception('Problem. See Press et al. pg 43')

    bet  = b[0]
    u[0] = r[0]/bet

    # Decompostion and forward subsitution
    for i in range(1,j-1):
        gam[i] = c[i-1]/bet
        bet    = b[i]-a[i]*gam[i]
        if bet == 0:
            raise Exception('TDMA failed')

        u[i] = (r[i]-a[i]*u[i-1])/bet

    # backsubsitution
    for i in range(j-2,0,-1):
        u[i] = u[i] - gam[i+1]*u[i+1]

    return u

Sr         = 0.03
rho_s      = 350.
rho_i      = 910.
rho_w      = 1000.
firn_param = 30                # y^-1
c_ice      = 2097              # [J kg-1 K-1]
spy        = 365.25*24*60*60.
L_heat     = 334000.0          # [J kg-1]

IC_fp = '../../initialization/coarse/result/glc1-b/nc/glc1-b_3000a_dt_1_dx_50_MB_-1.675_OFF_spline_k2.nc'

with xr.open_dataset(IC_fp) as src:
    # correct for minimum ice thickness
    src["height"] = xr.where(src.height <= 10, 0, src.height)
    # apply sigma coordinate transform for vertical coordinate
    src["Z"]      = src.zbed + src.Z * src.height
    # Calculate the magnitude of the velocity vectors
    src['vel_m']  = np.sqrt(src['velocity 1']**2 + src['velocity 2']**2)

# dictionary of PDD model parameters
p = [ 8.29376332e-05, -3.45256005e-02,  6.31076200e+00]
doy = np.arange(1,366)
fancy_std = np.polyval(p, doy)[:, np.newaxis]

const =  dict(T_m = 0.0, T_rs = 1.0, α = 10.8, T_ma = -9.02,
              ΔTΔz  = 6.5E-3, T_p = 196, ref_z = 2193, T_σ = fancy_std)

# initialize the PDD melt model classss
PDD   = PositiveDegreeDayModel(**const)

params = dict(f_s=4.1, f_i=8.2, f_r=0.1, grad_A=2.5e-4, A_mean=350)

mask = np.where(src.height.isel(coord_2=-1, t=-1).values > 0)
z_elmer  = src.Z.isel(t=-1, coord_2=-1).values[mask]
A, M = PDD.forward(z_elmer, **params)


# plt.plot(z_elmer, (A-M)/910); plt.show()

# Set some ~random firn thickness for now
FirnNodes = (A/910) / (firn_param/300)

#-------------------------------------------------------------------------------
# Begin 1D firn model Adrien's model
#-------------------------------------------------------------------------------
# run the test for just one node, the 100th surface node (should be accumulation zone)
surf_idx =  31 #81
# firn_idx = surf_idx - mask[0][0]
# surface elevation profile [m a.s.l.]
zprof_1D = src.Z.isel(t=-1, coord_1=surf_idx+mask[0][0]).values
# number of vertical elmer nodes (will be n+1 of the number of vertical layers)
n_elmer  = len(zprof_1D)
# calculate ice thickness from elmer
thick    = np.abs(zprof_1D[n_elmer-1] - zprof_1D[0])
# average vertical elmer gridcell spacing
dz_elmer = thick/(n_elmer-1)
# initial firn dz
dz_ini   = 0.06
# initial "firn"/seasonal snow profile thickness
thick1D  = 10.0
# firn thickness
firn_thick = FirnNodes[surf_idx]

print('Nodal Firn Thickness:')
# if elmer ice-thickness less than firn thickness
if thick < thick1D:
    n1D = int(np.floor(thick/dz_ini)+1)
    dz  = thick/(n1D-1)
# elemr thickness greater than firn column, set profile to 10 m.
else:
    n1D   = int(np.floor(thick1D/dz_ini)+1)
    dz    = thick1D/(n1D-1)
    thick = thick1D

# store thickness at beginning of time loop to keep track of
thick_ref = thick
dz_ref    = dz

# initialize empty arrays
dens       = np.zeros(n1D)    # [kg m-3]
dens_ref   = np.zeros(n1D)    # [kg m-3]
temp       = np.zeros(n1D)    # [K]
water      = np.zeros(n1D)    # [kg m-2]
LWC        = np.zeros(n1D)    # [m^3]
refreezing = np.zeros(n1D)    # [?]
# Matrix arrays for numerical solutions
alpha      = np.zeros(n1D)
beta       = np.zeros(n1D)
gamma      = np.zeros(n1D)
delta      = np.zeros(n1D)
diag       = np.zeros(n1D)
diag_L     = np.zeros(n1D)
diag_U     = np.zeros(n1D)
RHS        = np.zeros(n1D)

#-------------------------------------------------------------------------------
# Initial density
#-------------------------------------------------------------------------------
d = 0.0

for i in range(n1D):
    # if mb is postive then there is firn!
    if (A-M)[surf_idx] >= 0.0:
        dens_ref[i] = rho_s + d/(firn_thick*2.0*rho_w)*(rho_i**2-rho_s**2)
    # otherwise set density as ice
    else:
        dens_ref[i] = rho_i

    # density can't excceed density of ice
    if dens_ref[i] > rho_i:
        dens_ref[i] = rho_i

    d = d+dz

dens = dens_ref.copy()

# dens[:] = 910.0
# fig, ax = plt.subplots()
# ax.plot(dens, zprof_1D[0] - np.arange(n1D)*dz);
# plt.show()

#-------------------------------------------------------------------------------
# initialize firn/snow thickness and temp profile
# also calculate mean annual snow fall
#-------------------------------------------------------------------------------
# initialize null snow values
snow   = 0.0
T_mean = 0.0

# Calculate yearly values for temp and snow
for doy in range(0, 365):

    # Calculate surface air temperature for the day of year
    T = PDD.calc_air_temp(z_elmer)[doy, surf_idx]
    # Calculate mean annual ari temperature
    T_mean = T_mean + T/365.

    # Calculate melt (only using f_s??)
    melt = (T-PDD.T_m)*params['f_s']
    if melt < 0: melt=0.0

    # Calculate precip as snow
    if T <= PDD.T_rs:
        accu = params["A_mean"]/365.0 + (1 +(zprof_1D[0] - PDD.ref_z)*params["grad_A"])

    # mean annual snow accumulation for node
    snow = snow + accu - melt

    # snow can't be negative
    if snow < 0.0:
        snow = 0.0

# m / yr * 1 /yr , ???
snow = (snow)

# initialize temp array  with just surface temp
temp[:]  = T_mean + 273.15
temp_10m = 0.0
# days over which we are modeling
doy_i   = 0
doy_ip1 = 365


N_Y = 10
dens_mat = np.zeros((n1D,365*N_Y))
elev_mat = np.zeros((n1D,365*N_Y))
temp_mat = np.zeros((n1D,365*N_Y))
doy_mat  = (np.arange(0,365*N_Y)[None,:]*np.ones(n1D)[:,None])


dt = 1
# time integration loop!
for doy in range(doy_i, doy_ip1*N_Y, dt):

    doy_365 = int(doy - ((doy/365) - (doy/365)%1)*365)
    # Calculate surface air temperature for the day of year
    T = PDD.calc_air_temp(z_elmer)[doy_365, surf_idx]

    # Calculate melt (only using f_s??) [kg m^2 yr^-1]
    melt = (T-PDD.T_m)*params['f_s']
    if melt < 0: melt=0.0

    # Calculate precip as snow
    if T <= PDD.T_rs:
        accu = params["A_mean"]/365.0 + (1 +(zprof_1D[0] - PDD.ref_z)*params["grad_A"])

    # daily amount of snow on surface
    snow = snow + accu - melt - snow/firn_param/365.

    if snow < 0.0: snow = 0.0

    # calculate the surface temp for the current timestep
    tsurf = T + 273.15
    if tsurf > 273.15:
        tsurf = 273.15

    #---------------------------------------------------------------------------
    # Thickness change for seasonal snow
    #---------------------------------------------------------------------------
    # store previous timesteps surface elevation and dz
    old_zsurf      = zprof_1D[0]+(thick-thick_ref)
    old_dz         = dz
    old_temp       = temp.copy()
    old_water      = water.copy()
    old_refreezing = refreezing.copy()

    # Calculate the new thickness [m] based on the mass balance?
    thick = thick + (accu - melt - (snow/365.)) * (1/rho_s)

    # if column less than 10 m (i.e. net abalation has occured), set back to 10 m
    if thick < thick_ref: thick = thick_ref

    # current timesteps surface elevation and gridcell spacing
    dz    = thick/(n1D-1)
    zsurf = zprof_1D[0]+(thick-thick_ref)

    z=zsurf

    # linearly interpolating b/w grids time.
    for i in range(n1D):
    # for i in range(2):

        # interpolate from old grid to new grid
        #--------------------------------------
        if (old_zsurf >= z):
            # gridcell below the old free surface
            ii = int(np.floor((old_zsurf-z)/old_dz))+1
            z2 = old_zsurf - (ii)*old_dz

            if ((z2 - z)<0.0):
                # case where free surface has increased in elevation
                w1 = 1.0-(z-z2)/old_dz
            else:
                # case where free surface has decreased in elevation
                w1 = (z2-z)/old_dz

            if ii >= n1D:
                ii=n1D-1

            # interpolate onto new grids (indexing differs from Adiren b/c
            # python is zero indexed)
            temp[i]       = old_temp[ii-1]*(1.0-w1)+w1*old_temp[ii]
            water[i]      = old_water[ii-1]*(1.0-w1)+w1*old_water[ii]
            refreezing[i] = old_refreezing[ii-1]*(1.0-w1)+w1*old_refreezing[ii]

        else:
            # gridcells above the old freesurface, realy only the first node
            if i != 0:
                print(i)
            temp[i]       = tsurf
            water[i]      = 0.0
            refreezing[i] = 0.0

        # interpolate from elmer to new grid
        #-----------------------------------
        if (zprof_1D[0] >= z):
            ii_ref = int(np.floor((zprof_1D[0]-z)/dz_ref))+1
            z3     = zprof_1D[0] - (ii_ref)*dz_ref

            if ((z3 - z) > 0.0):
                w1 = (z3-z)/dz_ref
            else:
                w1 = 1.0-(z-z3)/dz_ref

            if ii_ref >= n1D:
                ii_ref=n1D-1

            dens[i] = dens_ref[ii_ref-1]*(1.0-w1)+w1*dens_ref[ii_ref]

        else:
            # pass
            dens[i] = rho_s

        # step to next surface elevation
        z = z - dz

    elev_mat[:,doy] = zsurf - np.arange(n1D)*dz
    dens_mat[:,doy] = dens.copy()


    # if doy == -1:
    #     # Conductivities
    #     K_im1 = 2.5e-6 * dens[0:-2]**2 - 1.23e-4 * dens[0:-2] + 0.024
    #     K_i   = 2.5e-6 * dens[1:-1]**2 - 1.23e-4 * dens[1:-1] + 0.024
    #     K_ip1 = 2.5e-6 * dens[2:  ]**2 - 1.23e-4 * dens[2:  ] + 0.024
    #
    #     # Diffusivities
    #     κ_im1 = K_im1 / (dens[0:-2] * c_ice) * 365.25
    #     κ_i   = K_i   / (dens[1:-1] * c_ice) * 365.25
    #     κ_ip1 = K_ip1 / (dens[2:  ] * c_ice) * 365.25
    #
    #     # variable Diffusivities
    #     κ_ipH = (κ_ip1 + κ_i)/2
    #     κ_imH = (κ_im1 + κ_i)/2
    #
    #     gamma[0]    = κ_im1[0] / dz**2
    #     gamma[1:-1] = κ_imH    / dz**2
    #     gamma[-1]   = (κ_ip1[-1]+κ_i[-1])/(2*dz**2)
    #     delta[0]    = (K_im1[0] + κ_i[1])/(2*dz**2)
    #     delta[1:-1] = κ_ipH    / dz**2
    #     delta[-1]   = κ_ipH[-1]/ dz**2
    #
    #     A = np.diagflat([-1*(gamma+delta)]) + \
    #         np.diagflat(gamma[1: ], k=-1) + \
    #         np.diagflat(delta[:-1], k= 1)
    #
    #     # Dirichlet BCs
    #     A[0,0] = 1
    #     A[0,1] = 0
    #
    #     # Neuman BCs
    #     A[-2,-1] =  gamma[-1] + delta[-1]
    #
    #     RHS[0] = 273.15
    #
    #     temp = LA.solve(A, RHS)
    #
    #     plt.plot(temp-273.15, zsurf - np.arange(n1D)*dz)
    #     plt.show()

    z = 0.0

    RHS = np.zeros(n1D)

    # Loop over interior layers
    for i in range(1,n1D-1):

        # Conductivities [W m^-1 K^-1]
        K_im1 = 2.5e-6 * dens[i-1]**2 - 1.23e-4 * dens[i-1] + 0.024
        K_i   = 2.5e-6 * dens[i]  **2 - 1.23e-4 * dens[i]   + 0.024
        K_ip1 = 2.5e-6 * dens[i+1]**2 - 1.23e-4 * dens[i+1] + 0.024

        # Diffusivities [m^2 / d] <---- [m^2 / s]
        κ_im1 = (K_im1 / (dens[i-1] * c_ice))*(60*60*24)
        κ_i   = (K_i   / (dens[i]   * c_ice))*(60*60*24)
        κ_ip1 = (K_ip1 / (dens[i+1] * c_ice))*(60*60*24)

        # variable Diffusivities
        κ_ipH = (κ_ip1 + κ_i)/2
        κ_imH = (κ_im1 + κ_i)/2

        # fourier number with variable kappas
        alpha[i] = (κ_ipH*dt)/(2*dz**2)
        beta[i]  = (κ_imH*dt)/(2*dz**2)

        # fill in componets of matrix A
        diag_L[i] =  -1 * beta[i]
        diag  [i] = ( 1 + alpha[i] + beta[i])
        diag_U[i] =  -1 * alpha[i]

        # fill componets of Force Vector (i.e. RHS)
        RHS[i] = beta[i]  * old_temp[i-1] + \
                 (1 - alpha[i] - beta[i]) * old_temp[i] + \
                 alpha[i] * old_temp[i+1]

    # Enforce surface boundary conditions
    #--------------------------------------
    alpha[0] = 0.0
    beta[0]  = 0.0
    # Dirichlet Boundary Condition at surface
    diag[0] = 1
    # Force Vector
    # RHS [0] = tsurf
    RHS[0] = tsurf

    # Enforce bottom of domain boundary conditions
    #---------------------------------------------
    alpha[-1] = (κ_i  *dt)/(2*dz**2)
    beta [-1] = (κ_ipH*dt)/(2*dz**2)
    # Neuman Boundary Condition at bottom of domain
    diag_L[-1] = -1 * beta[-1]
    diag  [-1] = (1 + beta[-1] )
    diag_U[-1] = 0.0
    # Force Vector
    RHS[-1]    = beta[-2] * old_temp[-2] + (1 - beta[-1]) * old_temp[-1]

    mat = np.diag(diag_U[0:-1], k=1 ) + \
          np.diag(diag,         k=0 ) + \
          np.diag(diag_L[1:  ], k=-1)


    temp = LA.solve(mat, RHS)

    # if temp.max() > 273.16:
    #     # print(dens[0], temp[0])
    #     print(f'{doy:>3}: {temp.max():.3f}, {temp.argmax()}')
    #     print(np.where(temp>273.16)[0])

    # print(doy, temp[0]- tsurf,)
#     temp = tridiag(diag_L, diag, diag_U, RHS)

    temp_mat[:,doy] = temp.copy()


    #---------------------------------------------------------------------------
    # CommunityFirnModel percolation / refreezing model
    #---------------------------------------------------------------------------

    # if melt > 0:
    #     print(doy)
    # LWC  = water * (1 / rho_w) * dz # [m]
    # pore close off
    rho_imp = 810                   # [kg m-3]

    # melt original calculated in kg m-2
    melt_vol_ie = melt / rho_i      # [m i.e.] <---- [kg m-2] / [kg m-3]
    melt_vol_we = melt / rho_w      # [m w.e.] <---- [kg m-2] / [kg m-3]

    # maybe this simplified to kg from kg m-2 b/c unit width and lenght
    melt_mass   = melt_vol_we * (rho_w)  # [kg m-2] of water i.e. og units

    # initialize arrays for calcs
    runofftot  = 0.
    LWCblocked = np.zeros(n1D)

    """
    TO DO: should claculate the cumulative mass in the column here to check
           for mass conservartion

           use cumsum??
    """
    # unit length [m], unit width [m], height = dz [m]
    mass = dens * dz * 1 * 1          # [kg] <--- [kg m^-3] * [ m * m * m]

    liq_in_mass = melt_mass + 0.0     # [kg]
    liq_in_vol  = liq_in_mass/rho_w   # [m3]

    # for mass conservation check
    liq_mc_init   = sum(LWC) + liq_in_vol #[m3]

    ###
    phi        = 1 - (dens / rho_i)      # porosity   [/]
    phi_vol    = phi * dz                # pore space [m]
    # i.e. tot potential pore space available for LWC; See Eq. 9 from Wever(2014)
    phi_vol_av = phi_vol * (rho_i/rho_w) # saturated water content [m w.e.]?
    # check for nodes were saturation could exceed ice density
    ilim       = np.where(dens + phi_vol_av * rho_w / dz  > rho_i)[0]

    if len(ilim) > 0: # limit pore space availability for storage in ilim nodes
        phi_vol_av[ilim] = np.maximum(dz*(916.99-dens[ilim])/rho_w, 0.)
        #https://github.com/UWGlaciology/CommunityFirnModel/blob/312fc30b62b7e36a609660e5b10e3269eb090bae/CFM_main/melt.py#L138

    LWCirr                 = Sr * phi_vol_av # volume of LWC that can be held as irreducible water [m w.e]
    LWCirr[dens > rho_imp] = 0.              # set irreducible water to zero where density > pore close off

    # LWC in excess of irreducible water conent, do be disrtributed
    LWC_excess = np.maximum(LWC-LWCirr, 0)

    cold_content = c_ice * mass * (273.15 - temp) # cold conten [J]
    # fix round_off errors
    cold_content = np.where(cold_content<1e-6, 0.0, cold_content)
    #
    refr_cap0      = (cold_content / L_heat) / rho_w   # refreezing capacity due to cold content, asuumes no existing water conent [m we]
    refr_cap0_supp = np.maximum(0, refr_cap0 - LWC)    # refreezing capacity for liquid water beyond currently present  [m we]
    refr_cap       = np.minimum(refr_cap0, phi_vol_av) # ?
    # (In theory, refr_cap should always be 0 for layers that have any LWC, except the upper most layer.)


    LWC_to_ice     = rho_w / rho_i * LWC                      # volume existing LWC would take if it froze [m we]
    refr_vol_supp  = np.maximum(0, phi_vol_av - LWC_to_ice)   # total volume available for additional LWC  [m we]
    refr_cap_supp  = np.minimum(refr_cap0_supp,refr_vol_supp) # refreezing capacity for additional LWC     [m we]

    rho_pot        = (mass + refr_cap * rho_w) / dz # potential denisty after refreezinf the maximum possible [kg m-3]
    phi_pot        = 1 - (rho_pot/rho_i)            # potential porosity  [/]
    cond1          = rho_pot > rho_i                # mask where potential density exceeds ice density
    phi_pot[cond1] = 0.0                            # [/]
    phi_vol_pot    = phi_pot * dz * 1 * 1           # potential pore space after refreezing [m^3]
    # i.e. tot potential pore space available for LWC; See Eq. 9 from Wever(2014)
    phi_vol_av_pot = phi_vol_pot * (rho_i/rho_w)   # saturated water content [m w.e.]?
    ilim           = np.where(rho_pot + phi_vol_av_pot * rho_w / dz > rho_i)[0] # nodes potentially exceeding 917 density

    if len(ilim) > 0:  # limit pore space availability for storage in ilim nodes
        phi_vol_av_pot[ilim] = np.maximum(dz * (916.99 - rho_pot[ilim]) / rho_w, 0.0)
        #https://github.com/UWGlaciology/CommunityFirnModel/blob/312fc30b62b7e36a609660e5b10e3269eb090bae/CFM_main/melt.py#L178


    LWCirr_pot = Sr * phi_vol_av_pot     # volume of LWC that can be held as irreducible water [m w.e]
    LWCirr_pot[rho_pot > rho_imp] = 0.   # set irreducible water to zero where density > pore close off


    LWCunf      = np.maximum(0,LWC-refr_cap)      # unfrozen LWC that will remain in each node after refreezing
    retcap_supp = np.maximum(0,LWCirr_pot-LWCunf) # retention capacity for additional LWC (assuming refreezing occurs first), potention minus how much there already is


    # Storage Capacity
    str_cap = refr_cap_supp + retcap_supp  # total storage capacity  for each node for additional LWC [m^3]

    ################################################
    # can add Conditions to work with icelenses here
    ################################################

    if np.any(dens < rho_imp):
        # the nodes below last rho<RhoImp are impermeable
        imp = np.arange(np.where(dens < rho_imp)[0][-1]+1,n1D,1).astype(int)
    else:
        # all nodes above impermeabilty threshold
        imp = np.arange(0,n1D,1).astype(int) # all nodes are impermeabless

    if len(imp) == 0:
        imp = [(len(dens) - 1)] # MS addition: If domain does not extend to full ice density, this will ensure the melt routine works (hack solution?)

    # https://fortran-lang.discourse.group/t/what-is-the-fastest-way-to-do-cumulative-sum-in-fortran-to-mimic-matlab-cumsum/1976/4
    str_cap[imp] = 0                  # set 0 storage capacity for impermeable nodes
    str_cap_sum  = np.cumsum(str_cap) # cumulative storage capacity, refreezing + irreducible


    ### Store surface melt according to storage capacity
    LWCblocked = np.zeros(n1D)

    # there is liquid water input at the surface
    if liq_in_vol > 0:
        # test if there is enough storage capacity for liquid input
        if str_cap_sum[-1] >= liq_in_vol:
            ii0 = np.where(str_cap_sum >= liq_in_vol)[0][0] # bottom most node storing surface melt
        else:
            ii0 = n1D - 1 # set to bottom node index


        if ii0 >= imp[0]:  # impermeable barrier or not enough pore space prevents full distribution of meltinput
            ii0             = max(0,imp[0]-1) # ii0 limited to node above impermeable barrier
            str_input       = np.concatenate((str_cap[0:ii0+1],np.zeros(n1D-ii0-1))) # each node above the barrier gets filled with its stcap
            LWCblocked[ii0] = liq_in_vol - sum(str_input) # volume of water that is excess due to blockage
        else:
            # all water input stored in surface node
            if ii0 == 0:
                str_input = np.concatenate(([liq_in_vol],np.zeros(n1D-1)))
            # water input is stored in several nodes ????
            else:
                str_input = np.concatenate((str_cap[0:ii0],[liq_in_vol-str_cap_sum[ii0-1]], np.zeros(n1D-ii0-1)))

    elif liq_in_vol == 0.0:
        str_input = np.zeros(n1D)


    str_cap = str_cap - str_input #update storage capcity

    ### Need to set LWC excess in imperiable nodes as blocked LWC for mass consv.
    # blocked_idxes = np.where(LWC_excess>0)

    LWC1    = np.copy(LWC)
    storage1 = np.zeros(n1D)

    if np.any(LWC_excess > 0.): # if there is some excess LWC
        tostore  = 0                             # LWC stock that must be stored
        idxs_exc = np.where(LWC_excess>0)[0]     # indices of nodes with excess LWC
        idx_f    = idxs_exc[-1]                  # bottom most node where LWC_excess exists
        idx_0    = idxs_exc[0]                   # top    most node where LWC_excess exists

        if np.any(str_cap > 0): # if any nodes have storage capacity
            idx_f_str = np.where(str_cap > 0)[0][-1] #bottom most node where LWC_excess can be stored
        else:
            idx_f_str = 0

        if ((np.any(str_cap[1:]>0) and (idx_f_str>idx_0))): #there is some storage capacity in the firn column, and it is deeper than upper most node where LWC_excess exists

            while ((idx_0 <= idx_f) or (tostore>0)):
                if (np.where(str_cap[idx_0:]>0)[0]).size > 0:
                    idx_1 = idx_0 + np.where(str_cap[idx_0:]>0)[0][0] # next node that can store some of the LWC_excess
                else: # all nodes with positive str_cap are shallower than idx_0, emulate the no storage capacity routine below
                    print("fuck")
                    break

                if imp[np.where(imp >= idx_0)[0][0]] > idx_1:  # idx_0 and idx_1 nodes not separated by an impermeable barrier
                    tostore         += np.sum(LWC_excess[idx_0:idx_1+1])  # all LWC_excess from idx_0 to idx_1 are subject to storage
                    LWC1[idx_0:idx_1+1] = np.minimum(LWC1[idx_0:idx_1+1],LWCirr[idx_0:idx_1+1]) # LWC_excess is evacuated
                    storage1[idx_1]  = np.minimum(str_cap[idx_1],tostore) # idx_1 node stores as much as possible
                    tostore         -= storage1[idx_1]             # tostore is reduced, idx_1 is filled
                    idx_0            = idx_1+1                     # go to next node with possible storage capacity

                    if idx_0 >= idx_f_str:  # no possible storage of LWC_excess
                        idx_1 = imp[np.where(imp>=idx_0)[0][0]] - 1 # find the next impermeable barrier
                        LWCblocked[idx_1] += tostore              # all LWC to be stored is blocked above the barrier
                        tostore = 0.                              # tostore is set to 0

                else:
                    print('THIS SHOULD BE THE CASE, BAD THINGS HAPPENIGN')

        else: # no storage capacity in the firn column
            for idx_0 in idxs_exc:  # find underlying impermeable barrier for each node with some LWC_excess
                idx_1 = imp[np.where(imp>=idx_0)[0][0]]-1   # idx_1 becomes index of node above the impermeable barrier
                LWCblocked[idx_1] += LWC_excess[idx_0]      # LWC_excess is blocked above the barrier
                LWC1[idx_0] = LWCirr[idx_0]                 # LWC of idx_0 is reduced to irreducible water content
    storagetot = str_input + storage1
    LWC1       = LWC1 + storagetot

    # Refreezing
    freeze = np.minimum(LWC1,refr_cap)            # refreezing in each individual node [m we]
    mass   = mass + rho_w*freeze                  # update mass [kg]
    LWC    = LWC1-freeze                          # update LWC  [m we]
    dens   = mass/dz                              # update density [kg m-3]
    lat_heat_relase = freeze * rho_w * L_heat     # latent heat released due to refreezing [J]
    cold_content = cold_content - lat_heat_relase # remaining cold content [J]
    refrozen_tot = sum(freeze)                    # total refreezing [m we]
    temp[freeze>0] = 273.15 - cold_content[freeze>0]/(c_ice*mass[freeze>0]) # update temp [K]


    ## CHECK FOR MASS conservation
    liq_mc_final = np.sum(LWC) + refrozen_tot + runofftot

    if np.abs(liq_mc_final - liq_mc_init) > 1e-4:
        print(f'Mass conservation error at step {doy}\n    Init:  {liq_mc_init} m\n    Final: {liq_mc_final} m')


    # ### Dry cold firn check ###
    # coldlayers = np.where(temp < 273.15)[0]
    # if np.all(LWC[coldlayers] < 1e-9):
    #     LWC[coldlayers] = 0.
    # if np.any(LWC[coldlayers] > 0.):
    #     print('Problem: water content in a cold layer')

    # Distrbute LWC_excess in the nodes support storage
    # #---------------------------------------------------------------------------
    # # percolation / refreezing model
    # #---------------------------------------------------------------------------
    # for i in range(n1D):
    #     """First loop accounts for the existing water content and the new solved
    #        for temperature. Make sure that with existing water content and new
    #        temperature we don't excceed 273.15 K
    #     """
    #     # account for water contents effect of enthalpy(or temp) of fusion
    #     # eqn 18 from gilbert et al. 2018, resolved for T from w
    #     temp[i] = temp[i]+water[i]*L_heat/(c_ice*dens[i])
    #
    #     if (temp[i] > 273.15):
    #         refreezing[i] = refreezing[i] + (water[i]-(temp[i]-273.15)*(c_ice*dens[i])/L_heat)*dz/rho_w
    #
    #         water[i] = (temp[i]-273.15)*(c_ice*dens[i])/L_heat
    #         temp[i]  = 273.15
    #     else:
    #         refreezing[i] = refreezing[i] + water[i]*dz/rho_w
    #
    #         water[i] = 0.0
    #
    # # melt [kg m-2], L_heat [J kg-1]
    # # convert to J m-2
    # cont = melt*L_heat #*rho_w
    #
    # j = 0
    #
    # while (cont != 0.0):
    #
    #     # pore close off
    #     if dens[j]>rho_imp:
    #         cont = 0.0
    #
    #     # effective temp
    #     #https://github.com/UWGlaciology/CommunityFirnModel/blob/312fc30b62b7e36a609660e5b10e3269eb090bae/CFM_main/melt.py#L803https://github.com/UWGlaciology/CommunityFirnModel/blob/312fc30b62b7e36a609660e5b10e3269eb090bae/CFM_main/melt.py#L803https://github.com/UWGlaciology/CommunityFirnModel/blob/312fc30b62b7e36a609660e5b10e3269eb090bae/CFM_main/melt.py#L803
    #     # cont [J m-2], convert to temp in [K]
    #     # assumes all melt refreezes and is converted to a temp change
    #     temp[j] = temp[j] + cont/(c_ice*dens[j]*dz)
    #
    #
    #     # temp has exceded 0 C,
    #     if temp[j] > 273.15:
    #         # units here seem wacky, but not sure if this is really needed
    #         refreezing[j] = refreezing[j] + (cont-(temp[j]-273.15)*(c_ice*dens[j]*dz)/(L_heat*rho_w))
    #
    #         # Amount of energy still available to melt [J m^-2]
    #         cont = (temp[j]-273.15)*(c_ice*dens[j]*dz)
    #
    #         # if water content has not reached saturation threshold
    #         if (water[j] < Sr*rho_w*(1-dens[j]/rho_i)):
    #
    #             # assume all melt is added to water conent [kg m-3]
    #             water[j] = water[j] + cont/dz/L_heat
    #
    #             # if newly updated water content exceed saturation threshold
    #             if (water[j] >= Sr*rho_w*(1-dens[j]/rho_i)):
    #                 # energy [J m-2] previous amount of energy minus maximum residual water conent
    #                 cont     = (water[j] - Sr*rho_w*(1-dens[j]/rho_i))*dz*L_heat
    #                 # set water conent to saturation threshold
    #                 water[j] = Sr*rho_w*(1-dens[j]/rho_i)
    #             else:
    #                 # all the melt has been consumed and there is no left to refreeze
    #                 cont = 0.0
    #
    #         #
    #         temp[j] = 273.15
    #     else:
    #         refreezing[j] = refreezing[j] + cont / (L_heat*rho_w)
    #         cont = 0.0
    #
    #
    #     # itterate to next vertical layer
    #     j = j+1
    #
    #     # reached the end of the vertical domain if true
    #     if (j==(n1D-1)):
    #         cont = 0.0



# fig, ax = plt.subplots(1,1, figsize=(6,3), constrained_layout=True)
#
# im = ax.contourf(doy_mat,
#                    elev_mat-zprof_1D[0],
#                    temp_mat-273.15,
#                    # shading='auto',
#                    extend='min',
#                    levels=np.linspace(-12,0,13),
#                    # clim=(-10,0)
#                    );
#
# ax.plot(doy_mat[33,:], elev_mat[33, :]-zprof_1D[0], lw=0.5, c='k')
# ax.plot(doy_mat[66,:], elev_mat[66, :]-zprof_1D[0], lw=0.5, c='k')
# ax.plot(doy_mat[99,:], elev_mat[99, :]-zprof_1D[0], lw=0.5, c='k')
# ax.plot(doy_mat[132,:], elev_mat[132, :]-zprof_1D[0], lw=0.5, c='k')
#
# fig.colorbar(im,ax=ax)
#
# ax.set_ylim(-thick1D, 1)
#
# fig.savefig(f'temp@node{surf_idx}.png', dpi=300, bbox_inches='tight', facecolor='w')


fig, ax = plt.subplots(2,1, sharey=True, sharex=True, constrained_layout=True)

im = ax[0].contourf(doy_mat,
                   elev_mat-zprof_1D[0],
                   temp_mat-273.15,
                   # shading='auto',
                   extend='min',
                   levels=np.linspace(-12,0,13),
                   # clim=(-10,0)
                   );
fig.colorbar(im,ax=ax[0])

im = ax[1].pcolormesh(doy_mat,
                   elev_mat-zprof_1D[0],
                   dens_mat,
                   cmap='plasma',
                   shading='auto');

fig.colorbar(im,ax=ax[1])

ax[0].set_ylim(-thick1D, 1)
#
#
plt.show()
