# Tuesday (08/24)  
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
