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

SUBROUTINE SurfaceMassBalance( Model,Solver,dt,TransientSimulation )
  ! *****************************************************************************
  !> SurfaceBoundary.f90, function SurfaceMassBalance
  !>
  !> Solve for the surface mass balance using a simple degree day approach which
  !> explicitly accounts for snow accumulation, melting, and refreezing (Gilbert et al. 2016).
  !>
  ! *****************************************************************************
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  IMPLICIT NONE
  TYPE(Model_t)              :: Model
  TYPE(Solver_t),    POINTER :: Solver
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Element_t),   POINTER :: Element
  INTEGER,           POINTER :: NodeIndexes(:)

  !-----------------
  ! Field Variables
  !-----------------
  TYPE(Variable_t),  POINTER :: MB             ! Mass balance         [m yr^{-1}]
  TYPE(Variable_t),  POINTER :: Dens           ! Density              []
  TYPE(Variable_t),  POINTER :: Firn           ! Firn thickness       [m]
  TYPE(Variable_t),  POINTER :: Depth          ! Ice Depth            [m]
  TYPE(Variable_t),  POINTER :: Melting        ! Surface Melting      [kg m^{-2} yr^{-1}]

  !--------------------
  ! Internal Variables
  !--------------------
  INTEGER :: n,i,j,k,cont,Day,ierr=0,year
  INTEGER :: nb_surf             ! number of surface  nodes
  INTEGER :: nb_vert             ! number of vertical nodes
  REAL(KIND=dp) :: dt            ! timestep size                      [yr]
  REAL(KIND=dp) :: z             ! surface elevation of current node  [m a.s.l.]
  REAL(KIND=dp) :: ref_z         ! reference surface elevation        [m a.s.l.]
  REAL(KIND=dp) :: z_precip      ! reference surf. elev. for precip   [m a.s.l.]
  REAL(KIND=dp) :: rho_w         ! water        density               [kg m^{-3}]
  REAL(KIND=dp) :: rho_i         ! ice          density               [kg m^{-3}]
  REAL(KIND=dp) :: rho_s         ! snow/surface density               [kg m^{-3}]
  REAL(KIND=dp) :: d_firn        ! firn regularization                [kg m^{-2} yr^{-2}]
  REAL(KIND=dp) :: T             ! surface air temperature            [K]
  REAL(KIND=dp) :: T_mean        ! mean annual surf. air. temp        [K]
  REAL(KIND=dp) :: A_mean        ! mean annual accumulation           [kg m^{-2} yr^{-1}]
  REAl(KIND=dp) :: alpha         ! annual surf. air temp. amplitude   [K]
  REAL(KIND=dp) :: grad_T        ! air temp lapse rate                [K m^{-1}]
  REAL(KIND=dp) :: grad_accu     ! precipitation lapse rate           [(kg m^{-2} yr^{-1}) m^{-1}]
  REAL(KIND=dp) :: thresh_melt   ! threshold temperature for melting  [K]
  REAL(KIND=dp) :: thresh_precip ! rain-snow percipitation threshold  [K]
  REAL(KIND=dp) :: r_s2m         ! ratio of snow accum. to melt       [-]
  REAL(KIND=dp) :: f_r           ! refreezing factor                  [-]
  REAL(KIND=dp) :: f_m           ! degree-day factor                  [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: f_ice         ! degree-day factor for ice          [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: f_snow        ! degree-day factor for snow         [kg m^{-2} yr^{-1} K^{-1}]
  REAL(KIND=dp) :: PDDs          ! positive degrees per day           [K]
  REAL(KIND=dp) :: accu_days     ! days where accum took place        [ ]
  REAL(KIND=dp) :: A_snow        ! surface accumulation               [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: R             ! rate of superimposed ice formation [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: M_melt        ! surface abalation                  [kg m^{-2} yr^{-1}]
  REAL(KIND=dp) :: melt_local    ! intermediate melt calculation      [kg m^{-2} yr^{-1}]

  logical :: found,GotIt,first_time=.true.,TransientSimulation

  !-------------------------
  ! First time loop
  !-------------------------
  save first_time,nb_surf,nb_vert

  if (first_time) then
    first_time=.false.

    do i=1,model % NumberOfNodes
      if (model % nodes % x(i+1)<model % nodes % x(i)) then
        exit
      endif
    enddo

    nb_surf=i
    nb_vert=model % NumberOfNodes/i

  endif
  !-------------------------

  ! Get Model Variables
  Firn    => VariableGet( Model % Variables, 'Firn')
  Dens    => VariableGet( Model % Variables, 'Densi')
  Depth   => VariableGet( Model % Variables, 'Depth')
  MB      => VariableGet( Model % Variables, 'Mass Balance')
  Melting => VariableGet( Model % Variables, 'Melting')

  ! Get Model Constants
  rho_i         = GetConstReal(Model % Constants, "rho_i")         ![kg m^{-3}]
  rho_w         = GetConstReal(Model % Constants, "rho_w")         ![kg m^{-3}]
  rho_s         = GetConstReal(Model % Constants, "rho_s")         ![kg m^{-3}]
  alpha         = GetConstReal(Model % Constants, "delta_T")       ![K]
  T_mean        = GetConstReal(Model % Constants, "T_mean" )       ![K]
  A_mean        = GetConstReal(Model % Constants, "A_mean" )       !TBD
  ref_z         = GetConstReal(Model % Constants, "ref_z"  )       ![m a.s.l.]
  d_firn        = GetConstReal(Model % Constants, "d_firn" )       ![kg m^{-2} y^{-1}]
  grad_T        = GetConstReal(Model % Constants, "grad_T" )       ![K m^{-1}]
  grad_accu     = GetConstReal(Model % Constants, "grad_accu" )    ![(kg m^{-2} yr^{-1}) m^{-1}]
  thresh_precip = GetConstReal(Model % Constants, "thresh_precip") ![K]
  thresh_melt   = GetConstReal(Model % Constants, "thresh_melt")   ![K]
  f_r           = GetConstReal(Model % Constants, "f_r")           ![-]
  f_ice         = GetConstReal(Model % Constants, "f_ice")         ![kg m^{-2} yr^{-1}]
  f_snow        = GetConstReal(Model % Constants, "f_snow")        ![kg m^{-2} yr^{-1}]
  z_precip      = GetConstReal(Model % Constants, "z_precip")      ![m a.s.l.]

  ! Itterate over all the nodes
  DO n=1,model % NumberOfNodes

    ! Depth is set to zero along free surface, so test if node is surface node
    if (Depth % Values (Depth % perm (n))==0.0) then
      ! Get Surface Elevation
      z = model % nodes % y(n)

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
      ! calculate local surface melt assuming f_m = f_snow
      melt_local = PDDs * f_snow
      ! calculate refreezing
      R = min(f_r*A_snow, melt_local)
      ! compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
      r_s2m = (A_snow - R*(1 + rho_w/((1-rho_s/rho_i)*rho_s)) ) / (melt_local)

      ! compute the degree-day factor
      if (r_s2m >= 1) then
        f_m = f_snow
      else
        f_m = f_ice - (f_ice - f_snow)*r_s2m
      endif

      ! calculate surface melt [kg m^{-2} yr^{-1}] with f_m
      M_melt = f_m*PDDs

      ! Set the mass balance [m yr^{-1}]
      MB % values (MB % perm(n)) = (A_snow + R - M_melt) * (1 / rho_i)
      ! Set the surface melt [kg m^{-2} yr^{-1}]
      Melting % values (Melting % perm(n)) = M_melt

      ! calculate the time dependent firn thickness
      ! NOTE: UNITS Inconsitent here
      Firn % values (Firn % perm(n)) = Firn % values (Firn % perm(n)) &
      + (r_s2m*melt_local-M_melt)*dt - Firn % values (Firn % perm(n))*0.3*dt/10.0
      ! fix negative firn thickness if ablation occurs
      if ( Firn % values (Firn % perm(n))  < 0.0) then
          Firn % values (Firn % perm(n)) = 0.0
      endif

      ! iterate over vertically aligned nodes
      do i=1,nb_vert
        ! index of ith vertically aligned node
        cont=n-(i-1)*nb_surf

        ! Check if within firn aquifer, if so set linear density profile
        if (  Firn % values (Firn % perm(n)) > 1.0) then
          ! NOTE: UNITS Inconsitent here
          Dens % values ( Dens % perm(cont)) = rho_s + Depth%values(Depth%perm(cont))/Firn%values(Firn%perm(n))*(rho_i - rho_s)
        else
          Dens % values ( Dens % perm(cont)) = rho_i
        endif

        ! Make sure density isn't greater than the ice density
        if (Dens % values ( Dens % perm(cont)) > rho_i) then
          Dens % values ( Dens % perm(cont)) = rho_i
        endif
      end do

    endif
  END DO
END SUBROUTINE SurfaceMassBalance

! function GetModelConstant(constant_name, GotIt)
!   constant = GetConstReal(Model % Constants, constant_name, GotIt)
!   if (.not. GotIt) then
!     call fatal('SurfaceMassBalance ---> GetModelConstant', 'Could not find ')
!   end if
! end function GetModelConstant


! FUNCTION getSurfaceEnthalpy(Model, Node, InputArray) RESULT(Enthalpy)
!   ! provides you with most Elmer functionality
!   USE DefUtils
!   ! saves you from stupid errors
!   IMPLICIT NONE
!   ! the external variables
!   !----------------------------------------------------------------------------
!   TYPE(Model_t) :: Model         ! the access point to everything about the model
!   INTEGER       :: Node          ! the current Node number
!   REAL(KIND=dp) :: InputArray(1) ! Contains the argument passed to the function
!   REAL(KIND=dp) :: Temp          ! intermediate result
!   REAL(KIND=dp) :: Enthalpy      ! the final result
!   !----------------------------------------------------------------------------
!   ! internal variables
!   !----------------------------------------------------------------------------
!   REAL(KIND=dp) :: lapserate
!   REAL(KIND=dp) :: intercept
!   REAL(KIND=dp) :: elevation
!   REAL(KIND=dp) :: T_ref, CapA, CapB
!
!   ! Enthalpy related constants
!   T_ref     = 200.0               ! [K]
!   CapA      = 7.253               ! [J kg-1 K-2]
!   CapB      = 146.3               ! [J kg-1 K-1]
!
!   ! lets hard-code our values (if we have time we can later make them being read from SIF)
!   lapserate = -0.0065_dp          ! [K m^{-1}] atmospheric lapse rate
!   intercept = 285.15              ! [K] temp at z=0
!   elevation = InputArray(1)       ! [m] elevation of current surface node
!
!   ! Calculate the temperature
!   Temp      = lapserate*elevation + intercept
!
!   ! calculate mean annual surface temperature
!   T_surf = grad*(alt-z)+Tmean+273.15
!   ! calculate the seasonal (6-month) surface temp
!   if (mod(floor(t_simu/dt), 2) == 1) then
!     T_surf = T_surf + alpha/2.0
!   else
!     T_surf = T_surf - alpha/2.0
!   endif
!   ! ! NOTE: Make sure this is correct .......
!   ! ! Temperature can't exced the melting point at the surface
!   ! if (T_surf > 273.15) then
!   !   T_surf = 273.15
!   ! endif
!
!
!   Enthalpy  = (CapA/2*(Temp**2 - T_ref**2) + CapB*(Temp-T_ref)) ! [J kg-1]
!
!   RETURN
!
! END FUNCTION getSurfaceEnthalpy
