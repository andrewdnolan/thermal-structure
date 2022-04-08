#!/usr/bin/env python3

import numpy as np
import xarray as xr
import scipy.linalg as LA
import matplotlib.pyplot as plt


def CFM_bucket(melt, temp, LWC, dens):

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

    return temp,
    # ### Dry cold firn check ###
    # coldlayers = np.where(temp < 273.15)[0]
    # if np.all(LWC[coldlayers] < 1e-9):
    #     LWC[coldlayers] = 0.
    # if np.any(LWC[coldlayers] > 0.):
    #     print('Problem: water content in a cold layer')
