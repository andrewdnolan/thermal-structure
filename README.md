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
