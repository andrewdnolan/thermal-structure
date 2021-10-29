### Thursday (10/28) Research Meeting w/ G.F. & A.A.:

__Mass Balance Grid Searches__:
  - No need to be doing these at course mesh resolutions, `dx=200` or larger is fine.
  - Prioritize finding an initial mass balance offset for each glacier,
    different flowlines doesn't really matter right now (secondary priority).

### Wednesday (09/29) Research Meeting w/ G.F.:

__Enthalpy Surface B.C.__:
  - Applying the temperature forcing as the discrete temperature at the time step
    isn't super intuitive to Gwenn. Suggests average temperature between the previous
    time step and the current time step.  
  - Months physically mean something, but numerically for properly sampling the
    temperature forcing, they don't make that much sense. Suggests 9 time steps
    per year, you accurately capture all the signal.

__Corrugated Temperature__:
  - Present in prognostic simulations but not diagnostic.
  - Looks like a numerical problem to Gwenn.
  - Suggests turning off coupling (so heat transfer through diffusion only, no advection)
    to start de-bugging

### Thursday (09/23) Research Meeting w/ G.F.:

__Synthetic Geometries__:  
  - With the addition of Glacier #18, the "small" glacier class still only has 2 glacier, whereas "medium" and "large" classes have three each.

  - _To Do_:        
    - [ ] Add another small glacier to make it 3 in each class, maybe glacier #13 from Crompton field study.

    - [ ] Clip the fore-field topography in the "all glaciers" plot to see if it
    helps with the crazy aspect ratio needed to make the plot readable.

__Enthalpy Model__:
  - Current test runs use a 6-month time step for the enthalpy equation.  
    - __GF__: If we are actually trying to resolve seasonal forcing wouldn't a 1
    or 2 month time step make more sense?
  - _Experiment_:
    - [ ] Compare Enthalpy results for 1-year, 6-month, 2-month, and 1-month time steps.

__Mass Balance Model__:
  - Linear firn density profile is OK, but feels a little weird. An exponential
    profile like Nat used would make more sense physically.
      - Look into briefly, but linear is working and not worth sinking too much
        time into.

---

### Friday  (09/10):

__NetCDf SaveGridData__:
  - No need to extend the `SaveGridData` Solver yourself, this has already been
    done and was committed to the `elmerice` branch of the git repo.

      - Will need to make sure you are on the  `elmerice` branch and `git pull`
        the most recent commits.

      - Might need to remake `Elmer` so as the NetCDF library is properly included
        and linked. (Check README compilation instructions for linking NetCDF)

---

### Tuesday (09/07)

__Higher Order Elements__:
 - Was able to successfully create a mesh with 2nd order quadrilateral elements, but
   `.sif` file failed to execute. First and foremost the current mass balance model
   loops over the vertically aligned nodes assuming 1st order quadrilateral nodes.
   If we want to pursue higher order elements, we will need to update that code.

 - Adrian's fancy surface B.C. model from his 2020 paper requires a vertically
   structured mesh. Seems like an unstructured mesh is off the table, but need to
   look into indexing for higher order elements.

---

### Monday (08/30)  

__Synthetic Geometries__:
  - We need the flow following coordinate system to start (i.e. x=0) at the end of
    flowline. In the cases where there are multiple flowlines up different tributaries
    they will have different lengths, so for a uniform coordinate system (at least
    over the overlapping parts) the coordinates need to start at the end of flowline.

---

### Tuesday (08/24)  
__Bedrock Depressions at Terminus__:
  - Fisher Glacier: Terminus is at about ~18.5 km from the end of flowline.   

__Swath Profiling__:

  - On 08/17 discussed synthetic glacier profiles with Gwenn. Swath profiling came
    up as an option to generate more accurate glacier profiles, especially for
    glaciers with more complex topography in the foreground.  

    - [Hergarten et al. 2014](https://esurf.copernicus.org/articles/2/97/2014/)
      seems to be a widely cited method for doing this. `c++` code (~1000 lines)
      is included in supplemental material but fails to compile on both OSX and linux. Seems like a feasible thing to implement in `python`, but given the limited
      importance of this probably not worth the time.

    - [Billy Armstrong has some code](https://github.com/kbarnhart/pySwath) for
      swath profiling, but doesn't appear to be super user friendly. Seems like
      more of an archive of personal script than things meant for public use.
