# thermal-structure

__To Do__:
  - Set up additional solvers, so that at each times step we have a record of the
    amount of heat contributed by each source term in the governing equation.
    - Diffusive Flux:
      - http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerModelsManual.pdf#page=229
      - How do I pass a variable diffusivity to the solver??
        - I can easily write another solver, but is there a way with the existing
          elmer variables?

  - Write a solver to calculate the peclet number as field variable for each timestep.
    - Should also write the "Brinkman Number" (see [Meyer and Minchew, 2018](https://www-sciencedirect-com.proxy.lib.sfu.ca/science/article/pii/S0012821X18303790?via%3Dihub#se0080) for example of it being used.)

  - Temperature Dependent Slip Coefficent?
    - Section 2.4.2 from [Cohen et al. 2018](https://tc.copernicus.org/articles/12/2515/2018/tc-12-2515-2018.pdf)


Notes from NWG:
  - from GEF during card ride:
    - we need to quantify how (and if) the changes in surge vigor during periodic surges are results of a less temperate area along the bed, or difference in driving stress

Timestepping:
  - Maybe we should be running all simulations with `dt <= 0.1` to account for seasonality.
    - Could have added benefit that noise from low convergence tolerances will superimposed on seasonal cycle, which should be generally be less noisy and more amenable to spectral filtering.


Maximum Water Content:
  - Allowing a maximum englacial water content of 3 % might help with some of the "noise" and instability looking things below the firn aquifer (in addition to the sub annual timesteps).
    - Also would allow for a more saturated colormap, and lessen the contrast to small variations in water content we are not comfortable interpreting given the discontinuous
      diffusivities and low numerical tolerances.
