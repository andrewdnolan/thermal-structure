program SurfaceMassBalance

  use netcdf
  use iso_fortran_env, only: dp => real64

  implicit none


  INTEGER :: i ! loop counter


  REAL(KIND=dp) :: &
  T,                       & ! nodal surface air temp.                [K]
  z,                       & ! surface elevation of current node      [m a.s.l.]
  ref_z         = 2300.0,  & ! reference surface elevation            [m a.s.l.]
  T_mean        = 264.15,  & ! Mean annual air temp @ ref_z           [K]
  A_mean        = 675.0,   & ! Mean annual accum.  @ ref_z            [kg m^-2 yr^-1]
  alpha         = 12.0,    & ! anual air temp. amplitude              [K]
  grad_A        = 0.333,   & ! accum. gradient                        [(kg m^{-2} yr^{-1}) m^{-1}]
  grad_T        = 6.5E-3,  & ! air temp lapse rate                    [K m^-1]
  thresh_melt   = 273.15,  & ! threshold temperature for melting      [K]
  thresh_precip = 271.15,  & ! rain-snow percipitation threshold      [K]
  temp_peak     = 365./2., & ! DOY of annual temp peak                [doy]
  PDDs,                    & ! cumulative positive degrees per year   [K]
  accu_days,               & ! Days w/ accumulation    in frac. years [a]
  A_snow,                  & ! surface accumulation                   [kg m^{-2} yr^{-1}]
  R,                       & ! rate of superimposed ice formation     [kg m^{-2} yr^{-1}]
  M_melt,                  & ! surface abalation                      [kg m^{-2} yr^{-1}]
  melt_local                 ! intermediate melt calculation          [kg m^{-2} yr^{-1}]

  logical :: found,GotIt,first_time=.true.,TransientSimulation

  ! Get Surface Elevation
  z = 2300

  accu_days = 0.0 ! Days where snow accumulation occured
  PDDs      = 0.0 ! Positive Degree per year (PDD) at node n

  ! Itterate over the julian calendar days
  do i=1,365
      ! Find surface temp for day(i)
      T=alpha*cos( 2*3.14*(i-temp_peak)/365 )+grad_T*(ref_z-z)+T_mean

      ! If temp. above then calculate the positive degrees for that day
      if (T>thresh_melt) then
        ! Equations (7) and (8) from Gilbert et al. 2016
        PDDs = PDDs + (T-thresh_melt)
      endif
      ! If temp. below the add one to the day tally
      if (T<thresh_precip) then
        accu_days=accu_days+1.0/365.0
      endif
  enddo

  !if accu_days (in fractional years) exceeds one there is problem
  if (accu_days > 1.0) write(*,*) accu_days

  ! calculate snow accumulation
  A_snow=(accu_days*A_mean)*(1 + (z-ref_z)*grad_A)
  ! calculate local surface melt assuming f_m = f_snow
  melt_local = PDDs * f_snow
  

end program SurfaceMassBalance
