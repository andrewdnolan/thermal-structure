! ******************************************************************************
! *
! *  Authors: Andrew Nolan, Adrian Gilbert
! *  Email:   anolan@sfu.ca
! *  Github:  andrewdnolan
! *
! *  Date Written:
! *   2021/07/01
! *
! *  Note:
! *   This script was adapted from material provided at the Elmer/Ice Beginners
! *   course at UiO in 2016. The original material can be found on Documentation
! *   page of the Elmer/Ice Wiki under "Course Material."
! ******************************************************************************

program SurfaceMassBalance
  ! *****************************************************************************
  !> SurfaceBoundary.f90, function SurfaceMassBalance
  !>
  !> Solve for the surface mass balance using a simple degree day approach which
  !> explicitly accounts for snow accumulation, melting, and refreezing (Gilbert et al. 2016).
  !>
  ! *****************************************************************************

  IMPLICIT NONE

  !--------------------
  ! Internal Variables
  !--------------------
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14)
  INTEGER :: n,i,j,k,cont,Day,ierr=0,year
  INTEGER :: nb_surf                        ! number of surface  nodes
  INTEGER :: nb_vert                        ! number of vertical nodes
  REAL(KIND=dp) :: z                        ! surface elevation of current node  [m a.s.l.]
  REAL(KIND=dp) :: ref_z         = 2300.0   ! reference surface elevation        [m a.s.l.]
  REAL(KIND=dp) :: z_precip      = 2300.0   ! reference surf. elev. for precip   [m a.s.l.]
  REAL(KIND=dp) :: rho_w         = 1000.0   ! water        density               [kg m^{-3}]
  REAL(KIND=dp) :: rho_i         = 910.0    ! ice          density               [kg m^{-3}]
  REAL(KIND=dp) :: rho_s         = 350.0    ! snow/surface density               [kg m^{-3}]
  REAL(KIND=dp) :: T                        ! surface air temperature            [K]
  REAL(KIND=dp) :: T_mean        = 265.97   ! mean annual surf. air. temp        [K]
  REAL(KIND=dp) :: A_mean        = 1100     ! mean annual accumulation           [kg m^{-2} yr^{-1}]
  REAl(KIND=dp) :: alpha         = 9.3      ! annual surf. air temp. amplitude   [K]
  REAL(KIND=dp) :: grad_T        = 0.006    ! air temp lapse rate                [K m^{-1}]
  REAL(KIND=dp) :: grad_accu     = 4e-4     ! precipitation lapse rate           [(kg m^{-2} yr^{-1}) m^{-1}]
  REAL(KIND=dp) :: thresh_melt   = 273.15   ! threshold temperature for melting  [K]
  REAL(KIND=dp) :: thresh_precip = 275.15   ! rain-snow percipitation threshold  [K]
  REAL(KIND=dp) :: r_s2m                    ! ratio of snow accum. to melt       [-]
  REAL(KIND=dp) :: f_r           = 2.0e-1   ! refreezing factor                  [-]
  REAL(KIND=dp) :: f_m           = 0.0      ! degree-day factor                  [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: f_ice         = 4.52e-3  ! degree-day factor for ice          [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: f_snow        = 1.56e-3  ! degree-day factor for snow         [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: PDDs                     ! positive degrees per day           [K]
  REAL(KIND=dp) :: accu_days                ! days where accum took place        [ ]
  REAL(KIND=dp) :: A_snow                   ! surface accumulation               [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: R                        ! rate of superimposed ice formation [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: M_melt                   ! surface abalation                  [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: melt_local               ! intermediate melt calculation      [kg m^{-2} yr^{-1}]

  logical :: found,GotIt,first_time=.true.,TransientSimulation

  ! Get Surface Elevation
  z = 2000

  accu_days = 0.0 ! Days where snow accumulation occured
  PDDs      = 0.0 ! Positive Degree per Day (PDD) at node n

  ! Itterate over the julian calendar days
  DO i=1,365
      ! Find surface temp for day(i)
      T=alpha*sin(2*3.14*i/365)+grad_T*(ref_z-z)+T_mean

      ! If temp. above then calculate the positive degrees for that day
      if (T>thresh_melt) then
        ! Equations (7) and (8) from Gilbert et al. 2016
        PDDs = PDDs + (T-thresh_melt)
      endif
      ! If temp. below the add one to the day tally
      if (T<thresh_precip) then
        accu_days=accu_days+1.0/365.0
      endif
  ENDDO

  ! calculate snow accumulation
  A_snow=accu_days*(A_mean+(z-z_precip)*grad_accu)

  write(*,*) A_snow,           "kg m^-2 a^-1"
  write(*,*) A_snow*(1/rho_i), "m a^-1"
  ! ! calculate local surface melt assuming f_m = f_snow
  ! melt_local = PDDs * f_snow
  ! ! calculate refreezing
  ! R = min(f_r*A_snow, melt_local)
  ! ! compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
  ! r_s2m = (A_snow - R*(1 + rho_w/((1-rho_s/rho_i)*rho_s)) ) / (melt_local)
  !
  ! ! compute the degree-day factor
  ! if (r_s2m >= 1) then
  !   f_m = f_snow
  ! else
  !   f_m = f_ice - (f_ice - f_snow)*r_s2m
  ! endif
  !
  ! ! calculate surface melt [kg m^{-2} yr^{-1}] with f_m
  ! M_melt = f_m*PDDs
  !
  ! ! Set the mass balance [m yr^{-1}]
  ! MB % values (MB % perm(n)) = (A_snow + R - M_melt) * (1 / rho_i)
  ! ! Set the surface melt [kg m^{-2} yr^{-1}]
  ! Melting % values (Melting % perm(n)) = M_melt
  !
  ! ! calculate the time dependent firn thickness
  ! ! NOTE: UNITS Inconsitent here
  ! Firn % values (Firn % perm(n)) = Firn % values (Firn % perm(n)) &
  ! + (r_s2m*melt_local-M_melt)*dt - Firn % values (Firn % perm(n))*0.3*dt/10.0
  ! ! fix negative firn thickness if ablation occurs
  ! if ( Firn % values (Firn % perm(n))  < 0.0) then
  !     Firn % values (Firn % perm(n)) = 0.0
  ! endif
end program SurfaceMassBalance
