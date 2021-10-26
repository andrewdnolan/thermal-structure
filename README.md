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

_tests_
