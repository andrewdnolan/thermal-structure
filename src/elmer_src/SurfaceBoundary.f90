! ******************************************************************************
! *
! *  Authors: Andrew Nolan
! *  Email:   anolan@sfu.ca
! *  Github:  andrewdnolan
! *
! *  Date Written:
! *   2021/07/01
! *
! ******************************************************************************

SUBROUTINE SurfaceMassBalance( Model,Solver,dt,TransientSimulation )
  ! *****************************************************************************
  !> SurfaceBoundary.f90, function SurfaceMassBalance
  !>
  !> Solve for the surface mass balance using a simple degree day approach which
  !> explicitly accounts for snow accumulation, melting, and refreezing (Gilbert et al. 2016).
  !>
  ! *****************************************************************************

  ! Elmer modules
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  ! Local modules
  USE SurfaceTemperature

  IMPLICIT NONE

  TYPE(Model_t)              :: Model
  TYPE(Solver_t),    POINTER :: Solver
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Element_t),   POINTER :: Element
  INTEGER,           POINTER :: NodeIndexes(:)

  !-----------------
  ! Field Variables
  !-----------------
  TYPE(Variable_t),  POINTER :: MB,      & ! Mass balance    [m yr^{-1}]
                                Dens,    & ! Density         []
                                Firn,    & ! Firn thickness  [m]
                                Depth,   & ! Ice Depth       [m]
                                Melting    ! Surface Melting [kg m^{-2} yr^{-1}]

  !--------------------
  ! Internal Variables
  !--------------------
  INTEGER       :: n,i,j,k,cont,ierr=0
  INTEGER       :: nb_surf,    & ! number of surface  nodes
                   nb_vert       ! number of vertical nodes
  ! reference model constants
  REAL(KIND=dp) :: dt,         & ! timestep size                [yr]
                   z,          & ! surf. elev. of current node  [m a.s.l.]
                   rho_w,      & ! water        density         [kg m-3]
                   rho_i,      & ! ice          density         [kg m-3]
                   rho_s         ! snow/surface density         [kg m-3]

  ! air temp related
  INTEGER       :: T_peak        ! DOY of annual temp peak      [DOY]
  REAL(KIND=dp) :: alpha,      & ! Anual air temp. amp          [K]
                   dTdz,       & ! air temp lapse rate          [K m-1]
                   T_mean,     & ! mean annual surf. air. temp  [K]
                   T(365),     & ! surface air temperature      [K]
                   std_c0,     &
                   std_c1,     &
                   std_c2

  ! degree day model related
  REAL(KIND=dp) :: z_ref,      & ! reference surface elevation  [m a.s.l.]
                   A_mean,     & ! mean annual accumulation     [kg m-2 yr-1]
                   dPdz,       & ! precipitation lapse rate     [(kg m-2 yr-1) m-1]
                   T_melt,     & ! melting temp. thresh         [K]
                   T_r2s,      & ! rain 2 snow precip threshold [K]
                   r_s2m,      & ! ratio of snow accum. to melt [-]
                   f_r,        & ! refreezing factor            [-]
                   f_m,        & ! degree-day factor            [kg m-2 yr-1 K-1]
                   f_ice,      & ! degree-day factor for ice    [kg m-2 yr-1 K-1]
                   f_snow        ! degree-day factor for snow   [kg m-2 yr-1 K-1]
  ! firn model related parameters
  REAL(KIND=dp) :: d_firn        ! firn regularization          [kg m-2 yr-2]

  ! values calculated within subroutine
  REAL(KIND=dp) :: PDDs,       & ! positive degrees per day     [K]
                   accu_days,  & ! days where accum took place  [y]
                   A_snow,     & ! surface accumulation         [kg m-2 yr-1]
                   R,          & ! superimposed ice formation   [kg m-2 yr-1]
                   M_melt,     & ! surface abalation            [kg m-2 yr-1]
                   melt_local    ! intermediate melt calc.      [kg m-2 yr-1]

  logical :: found,GotIt,TransientSimulation,first_time=.true.

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
  rho_i  = GetConstReal(Model % Constants, "rho_i")    ![kg m-3]
  rho_w  = GetConstReal(Model % Constants, "rho_w")    ![kg m-3]
  rho_s  = GetConstReal(Model % Constants, "rho_s")    ![kg m-3]
  alpha  = GetConstReal(Model % Constants, "alpha")    ![K]
  T_mean = GetConstReal(Model % Constants, "T_mean" )  ![K]
  T_peak = INT(ANINT(GetConstReal(  Model % Constants, "T_peak" )))  ![DOY]
  A_mean = GetConstReal(Model % Constants, "A_mean" )  !TBD
  z_ref  = GetConstReal(Model % Constants, "z_ref"  )  ![m a.s.l.]
  d_firn = GetConstReal(Model % Constants, "d_firn" )  ![kg m-2 y-1]
  dTdz   = GetConstReal(Model % Constants, "dTdz" )    ![K m-1]
  dPdz   = GetConstReal(Model % Constants, "dPdz" )    ![(kg m-2 yr-1) m-1]
  T_r2s  = GetConstReal(Model % Constants, "T_r2s")    ![K]
  T_melt = GetConstReal(Model % Constants, "T_melt")   ![K]
  f_r    = GetConstReal(Model % Constants, "f_r")      ![-]
  f_ice  = GetConstReal(Model % Constants, "f_ice")    ![kg m-2 yr-1]
  f_snow = GetConstReal(Model % Constants, "f_snow")   ![kg m-2 yr-1]
  std_c0 = GetConstReal(Model % Constants, "std_c0")   ![K?]
  std_c1 = GetConstReal(Model % Constants, "std_c1")   ![K?]
  std_c2 = GetConstReal(Model % Constants, "std_c2")   ![K?]

  ! Itterate over all the nodes
  DO n=1,model % NumberOfNodes

    ! Depth is set to zero along free surface, so test if node is surface node
    if (Depth % Values (Depth % perm (n))==0.0) then
      ! Get Surface Elevation
      z = model % nodes % y(n)

      accu_days = 0.0 ! Days where snow accumulation occured
      PDDs      = 0.0 ! Positive Degree per Day (PDD) at node n

      ! Calculate air temperature curve
      call SurfTemp(z, T, alpha, dTdz, z_ref, T_mean, T_peak)

      ! Itterate over the julian calendar days
      DO i=1,365
          ! If temp. above then calculate the positive degrees for that day
          if ( T(i) > T_melt ) then
            ! Equations (7) and (8) from Gilbert et al. 2016
            PDDs = PDDs + ( T(i) - T_melt )
          endif
          ! If temp. below the add one to the day tally
          if ( T(i) < T_r2s ) then
            accu_days=accu_days + 1.0_dp / 365.0_dp
          endif
      ENDDO

      ! if accu_days (in fractional years) exceeds one there is problem
      if (accu_days > 1.0) write(*,'(a)') "Warning: accu_days exceeds 1"

      ! calculate snow accumulation
      A_snow = max((accu_days*A_mean)*(1 + (z-z_ref)*dPdz), 0.0)
      ! calculate local surface melt assuming f_m = f_snow
      melt_local = PDDs * f_snow
      ! calculate refreezing
      R = min(f_r*A_snow, melt_local)
      ! compute the ratio b/w accumulated snow and total melt assuming f_m = f_snow
      r_s2m = A_snow / melt_local

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
      Melting % values (Melting % perm(n)) = M_melt / rho_w

      ! calculate the time dependent firn thickness
      ! NOTE: assuming all new firn is added with the desnity of ice. Inconsistent
      !       with Gilbert et al. 2020, check with Gwenn about this
      Firn % values (Firn % perm(n)) = Firn % values (Firn % perm(n)) &
            + ( (A_snow - R  - M_melt) * (1/rho_i) - 0.03 * Firn % values (Firn % perm(n)) )*dt

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
          Dens % values ( Dens % perm(cont)) = rho_s + (Depth%values(Depth%perm(cont))/Firn%values(Firn%perm(n)))*(rho_i - rho_s)
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

FUNCTION getSurfaceEnthalpy(Model, Node, InputArray) RESULT(Enthalpy)
  ! Elmer modules
  USE DefUtils
  ! Local modules
  USE SurfaceTemperature

  IMPLICIT NONE

  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model         ! the access point to everything about the model
  INTEGER       :: Node          ! the current Node number
  REAL(KIND=dp) :: InputArray(1) ! Contains the argument passed to the function
  REAL(KIND=dp) :: Enthalpy      ! the final result
  TYPE(Variable_t), POINTER :: TimeVar
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: Time,        & ! simulation time                    [a]
                   Timestep(1), & ! vector of timestep size            [a]
                   z              ! surface elevation of current node  [m a.s.l.]
  integer       :: DOY            ! number of timesteps in a year
  logical       :: GotIt,       &
                   Found,       &
                   Transient

  ! air temp related
  INTEGER       :: T_peak        ! DOY of annual temp peak      [DOY]
  REAL(KIND=dp) :: z_ref,      & ! reference surface elevation  [m a.s.l.]
                   alpha,      & ! Anual air temp. amp          [K]
                   dTdz,       & ! air temp lapse rate          [K m-1]
                   T_mean,     & ! mean annual surf. air. temp  [K]
                   T(365),     & ! surface air temperature      [K]
                   T_surf,     & ! intermediate result
                   std_c0,     &
                   std_c1,     &
                   std_c2
  ! Enthalpy related params
  REAL(KIND=dp) :: T_ref, CapA, CapB

  ! Figure out if simulation is diagnostic or prognostic
  Transient = GetString(GetSimulation(), "Simulation type", Found)=='transient'

  ! Enthalpy related constants
  T_ref  = GetConstReal(Model % Constants, "t_ref_enthalpy")            ! [K]
  CapA   = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A")  ! [J kg-1 K-2]
  CapB   = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")  ! [J kg-1 K-1]
  ! Air temperature related constants
  alpha  = GetConstReal(Model % Constants, "alpha")                 ! [K]
  dTdz   = GetConstReal(Model % Constants, "dTdz" )                 ! [K m^{-1}]
  T_mean = GetConstReal(Model % Constants, "T_mean")                ! [K]
  T_peak = INT(ANINT(GetConstReal( Model % Constants, "T_peak" )))  ! [DOY]
  z_ref  = GetConstReal(Model % Constants, "z_ref" )                ! [m a.s.l.]
  std_c0 = GetConstReal(Model % Constants, "std_c0")                ! [K?]
  std_c1 = GetConstReal(Model % Constants, "std_c1")                ! [K?]
  std_c2 = GetConstReal(Model % Constants, "std_c2")                ! [K?]
  ! Firn Heating Parameters


  ! unpack surface elevation of current surface node [m]
  z = InputArray(1)

  ! Calculate air temperature curve
  call SurfTemp(z, T, alpha, dTdz, z_ref, T_mean, T_peak)

  if (Transient) then
    ! if transient get air temp for doy of timestep
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    Time     =  TimeVar % Values(1)
    ! Convert real time into approximate DOY
    DOY    = NINT((Time - floor(Time))*365.0)
    ! Get surface temp for DOY
    T_surf = T(DOY) + 273.15 ! K <-- C
  else
    ! if steady state return mean annual air temp
    T_surf = SUM(T)/365 + 273.15 ! K <-- C
  endif

  ! Temperature can't exced the melting point at the surface
  if (T_surf > 273.15) then
    T_surf = 273.15
  endif

  Enthalpy  = (CapA/2.0*(T_surf**2 - T_ref**2) + CapB*(T_surf-T_ref)) ! [J kg-1]

  RETURN

  contains

  ! subfunction for reading constants and error checking
  function GetModelConstant(Model, constant_name, GotIt) result(constant)
    USE DefUtils
    implicit none

    real :: constant
    logical, intent(in) :: GotIt
    TYPE(Model_t), intent(in) :: Model
    character(len=2000), intent(in) :: constant_name

    constant = GetConstReal(Model % Constants, trim(constant_name), GotIt)

    if (.not. GotIt) then
      call fatal('getSurfaceEnthalpy ---> GetModelConstant', &
                 'Could not find '//trim(constant_name) )
    end if

  end function GetModelConstant
END FUNCTION getSurfaceEnthalpy

SUBROUTINE Surface_Processes( Model, Solver, dt, TransientSimulation)
  !------------------------------------------------------------------------------
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

  ! Local modules
  USE SurfaceTemperature
  USE PositiveDegreeDays

  IMPLICIT NONE      ! saves you from stupid errors

  ! Solver name
  CHARACTER(*), PARAMETER :: Caller = 'Surface_Processes'
  !----------------------------------------------------------------------------
  ! external variables
  !----------------------------------------------------------------------------
  LOGICAL                    :: TransientSimulation
  REAL(KIND=dp)              :: dt
  TYPE(Model_t)              :: Model
  TYPE(Solver_t)             :: Solver
  TYPE(Element_t),   POINTER :: Element
  TYPE(ValueList_t), POINTER :: Params

  TYPE(Variable_t), POINTER  :: Surf_Enthalpy, &  ! Surface  Enthalpy     [J kg-1]
                                Enth,          &  ! Enthalpy              [J kg-1]
                                Q_lat_vol,     &  ! Latent (volumetric) Heat source  [J m-3 y-1]
                                H_f,           &  ! Phase Change Enthalpy [J kg-1]
                                Depth,         &  ! Depth below surf      [m]
                                Height,        &  ! Height above bed      [m]
                                TimeVar,       &  ! Time                  [yr]
                                MB,            &  ! Mass Balance (i.e.)   [m yr]
                                RunOff,        &  ! fraction of melt that runsoff [-]
                                surf_melt,     &  ! surface melt from PDD [m]
                                Dens              ! Ice Depth             [m]

  !----------------------------------------------------------------------------
  ! Local variables
  !----------------------------------------------------------------------------
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Misc. internal integers
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  INTEGER       :: i, n, k,     &  ! index counter
                   cont,        &  ! vertically alligned node index counter
                   e,           &  ! element counter
                   d,           &  ! julian calendar day  counter
                   N_n,         &  ! number of model    nodes
                   N_s,         &  ! number of surface  nodes
                   N_v,         &  ! number of vertical nodes
                   N_years,     &  ! number of years in timsetep      [a]
                   doy_i,       &  ! day of year of current timestep  [doy]
                   doy_ip1         ! day of year of next timesteps    [doy]
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! air temp related
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  INTEGER       :: T_peak,      & ! DOY of annual temp peak      [DOY]
                   seed=123456789 ! seed for random num. generator
  REAL(KIND=dp) :: z_ref,       & ! reference surface elevation  [m a.s.l.]
                   alpha,       & ! Anual air temp. amp          [K]
                   dTdz,        & ! air temp lapse rate          [K m-1]
                   T_mean,      & ! mean annual surf. air. temp  [C]
                   T(365),      & ! surface air temperature      [K]
                   std(365),    & ! std. of surf air temp        [K]
                   T_surf,      & ! intermediate result          [K]
                   H_surf,      & ! intermediate result          [J kg-1]
                   std_c0,      &
                   std_c1,      &
                   std_c2
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Surface heating (i.e. latent heat) related variables
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  REAL(KIND=dp) :: z,           &  ! nodal surface elev.   [m a.s.l.]
                   Time,        &  ! Simulation time       [a]
                   T_melt,      &  ! Melting threshold     [K]
                   L_heat,      &  ! Latent heat of fusion [J kg-1]
                   f_dd,        &  ! Degree-day factor     [m K^-1 a^-1]
                   h_aq,        &  ! Firn aquifer thick    [m]
                   rho_i,       &  ! denisty of ice        [Kg m^-3]
                   rho_w,       &  ! denisty of water      [Kg m^-3]
                   rho_s,       &  ! denisty of snow       [Kg m^-3]
                   rho_f,       &  ! pore close off dens.  [kg m^-3]
                   rho,         &  ! nodal density         [kg m^-3]
                   C_firn,      &  ! Surf. Dens. const.    [-]
                   z_trans,     &  ! firn/ice trans. depth [m]
                   Melt,        &  ! surf. (snow) melting  [m s.e.]
                   water,       &  ! water content         [-]
                   w_res,       &  ! residule water cont.  [-]
                   dz,          &  ! vert. layer thickess  [m]
                   Q_lat,       &  ! latent heat source.   [J m-3 yr-1]
                   Q_max,       &  ! Max. possible Q_lat   [J m-3 yr-1]
                   pump,        &  !
                   PDD(365) = 0    ! + degrees. for DOY    [K]
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Enthalpy related params
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  REAL(KIND=dp) :: T_ref,       &  ! reference temp.      [K]
                   CapA,        &  ! Heat cap. const 1    [J kg-1 K-2]
                   CapB,        &  ! Heat cap. const 1    [J kg-1 K-1]
                   w_max_aq,    &  ! max water content in frin aq  [-]
                   w_max_en,    &  ! max water content in glac ice [-]
                   Enthalpy_max    ! Maximum Enthalpy     [J kg-1]
  ! REAL(KIND=dp), allocatable    &
  !               :: Enthalpy_max(:) ! Maximum Enthalpy     [J kg-1]

  ! How the heat souce is treated. Either Mass (J/kg) or Volumetric (J/m)
  character(len=max_name_len) :: Heat_Source ! How the heat souce is treated

  ! how to set the maximum englacial water content
  CHARACTER(*), PARAMETER :: omega_maximum = 'denisty'

  LOGICAL :: GotIt, Seasonality, Source_at_Surface, first_time=.true., percolate

  ! Variables to keep track of between calls to the solver
  save first_time,N_v,N_s,N_n


  !-------------------------
  ! First time loop
  !-------------------------
  ! Variable declaration that only needs to be run once (i.e. first time)
  IF (first_time) THEN
    ! set first_time to false for all subsequent time steps
    first_time=.false.

    ! Loop over all model nodes
    DO i=1,model % NumberOfNodes
      IF (model%nodes%x(i+1) < model%nodes%x(i)) THEN
        EXIT
      ENDIF
    ENDDO
    N_n = Model % NumberOfNodes ! Number of Nodes in Models
    N_s = i                     ! Number of surface nodes
    N_v = N_n/i                 ! Number of vertical nodes
  ENDIF
  ! End "first_time" loop

  ! Solver passed params
  Params => GetSolverParams()

  ! Pointer to the model variables
  Depth         => VariableGet( Model % Variables, "Depth")  ! [m]
  Height        => VariableGet( Model % Variables, "Height") ! [m]
  Dens          => VariableGet( Model % Variables, "Densi")  ! [kg m-3]
  MB            => VariableGet( Model % Variables, "mass balance") ! [m y-1]
  Enth          => VariableGet( Model % Variables, "Enthalpy_h")   ! [J kg-1]
  H_f           => VariableGet( Model % Variables, "Phase Change Enthalpy") ! [J kg-1]
  Surf_Enthalpy => VariableGet( Model % Variables, "Surface_Enthalpy") ![J kg-1]
  Q_lat_vol     => VariableGet( Model % Variables, "Q_lat")       ![J m-3 yr-1]
  surf_melt     => VariableGet( Model % Variables, "surf_melt")   ![m]
  RunOff        => VariableGet( Model % Variables, "runoff_frac") ![-]

  ! Physical Parameters
  rho_i       = GetParam(Model,  "rho_i")                   ![kg m-3]
  rho_s       = GetParam(Model,  "rho_s")                   ![kg m-3]
  rho_w       = GetParam(Model,  "rho_w")                   ![kg m-3]
  ! Physical constants
  L_heat      = GetParam(Model, "L_heat")                   ![J kg-1]
  T_ref       = GetParam(Model, "t_ref_enthalpy")           ! [K]
  CapA        = GetParam(Model, "Enthalpy Heat Capacity A") ! [J kg-1 K-2]
  CapB        = GetParam(Model, "Enthalpy Heat Capacity B") ! [J kg-1 K-1]
  ! Melt Parameters
  f_dd        = GetParam(Model, "f_dd")                     ![m K-1 yr-1]
  T_melt      = GetParam(Model, "T_melt")                   ![K]
  ! Firn Aquifer Parameters
  rho_f       = GetParam(Model, "rho_f")                    ![kg m-3]
  h_aq        = GetParam(Model, "h_aq")                     ![m]
  C_firn      = GetParam(Model, "C_firn")                   ![-]
  w_max_aq    = GetParam(Model, "w_max_aq")
  w_max_en    = GetParam(Model, "w_max_en")
  ! Air temperature related Parameters
  alpha       = GetParam(Model, "alpha")                    ! [K]
  dTdz        = GetParam(Model, "dTdz" )                    ! [K m^{-1}]
  T_mean      = GetParam(Model, "T_mean")                   ! [K]
  T_peak      = INT(ANINT(GetParam(Model, "T_peak" )))      ! [DOY]
  z_ref       = GetParam(Model, "z_ref" )                   ! [m a.s.l.]
  std_c0      = GetParam(Model, "std_c0")                   ! [K?]
  std_c1      = GetParam(Model, "std_c1")                   ! [K?]
  std_c2      = GetParam(Model, "std_c2")                   ! [K?]

  ! Determine if there should be seasonality in airtemp forcing
  Seasonality = ListGetLogical(Model % Constants, "Seasonality", GotIt )
  IF ( .NOT.GotIt ) then
    ! default to true
    Seasonality = .true.
    CALL INFO(Caller,"No entry for >Seasonality<. in Constant Section", Level=6)
    CALL INFO(Caller,"Setting to .true. ",     Level=6)
  end if

  ! Determine how to treat the latent heat from meltwater refreezing
  Heat_Source = GetString( Params, 'Latent Heat Source', GotIt )
  IF (GotIt) then
    ! in this case, whatever was passed was not a valid option, so just using default value
    if ((TRIM(Heat_Source) .ne. "volumetric") .and. &
        (TRIM(Heat_Source) .ne. "mass")     ) then
      Heat_Source = "volumetric"
      CALL INFO(Caller,"Entry found for >Heat_Source<. is not a valid option.", Level=6)
      CALL INFO(Caller,"Setting to 'volumetric'",     Level=6)
    endif
  ELSE
    ! in this case, no option was passed at all, so just using default value
    Heat_Source = "volumetric" ! or "mass"
    CALL INFO(Caller, "No entry found for >Heat_Source<. in Solver Section", Level=6)
    CALL INFO(Caller, "Setting to 'volumetric'",     Level=6)
  END IF

  ! Determine where to add the volumetric heat source, keyword really only
  ! matters when Heat_Source = "volumetric"
  Source_at_Surface = GetLogical( Params, 'Source at Surface', GotIt )
  IF ( .NOT.GotIt ) then
    ! default to one below the free surface
    Source_at_Surface = .false.
    CALL INFO(Caller,"No entry for >Source at Surface<. in Solver Section", Level=6)
    CALL INFO(Caller,"Setting to .false. ",     Level=6)
  end if

  if (TransientSimulation) then
    ! if transient get current timestep, which is really time at end of timestep
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    ! Get current time
    Time    =   TimeVar % Values(1)

    ! Check that a valid timestep is passed
    if ((dt .ge. 1.0) .and. (MOD(dt, 1.0) /= 0.0)) then
      call fatal(Caller, 'For dt > 1, only full year strides are supported. E.g. modulo(dt, 1.0) == 0')
    end if

    ! Find the number of years in the timestep
    if (dt .ge. 1.0) then
      N_years = floor(dt)
    else
      N_years = 1.0
    end if

    ! find DOY at begining of timestep [d]
    doy_i = NINT(MOD(time-dt, 1.0) * 365.0)

    ! B/C round off error ith doy is sometimes 365 or 0 instead of 1. Check and fix
    if ((doy_i==365) .or. (doy_i==0)) then
      doy_i = 1
    end if

    ! Find DOY at end of timestep [d]
    doy_ip1   = NINT(MOD(time, 1.0) * 365.0)

    ! B/C round off error ith+1 doy is sometimes 0 instead of 365. Check and fix
    if (doy_ip1==0) then
      doy_ip1 = 365
    end if

    ! yearly timesteps give inf mean annual air temp due to division by zero
    ! also don't prodcuce any melt by default
    if ((doy_ip1-doy_i .eq. 0) .and. (MOD(dt, 1.0) .eq. 0.0)) then
      doy_i   = 1
      doy_ip1 = 365
    end if
  else
    ! set the "timestep" as one year
    dt = 1.0
    ! if steady state get yearly amount of melt
    doy_i   = 1
    doy_ip1 = 365
  endif

  ! Create std(doy) array (only needs to be done once)
  call temp_std(std, (/ std_c0, std_c1, std_c2 /))

  ! Outter Most Loop: Itterate over model nodes
  DO n=1,N_n

    ! Check if depth == 0.0 for node, i.e. is it a surface node
    IF (Depth%Values(Depth%perm(n)) == 0.0) THEN

      ! Get nodal surface elevation [m]
      z = model%nodes%y(n)

      ! Calculate nodal air temperature curve as function of day of year
      call SurfTemp(z, T, alpha, dTdz, z_ref, T_mean, T_peak)

      ! Calculte PDDs from semi-analytical solution from Calvo and Greve 2005
      PDD = Cavlo_Greve_PDDs(T,std) * N_years    ! [K d]

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Calculate surface heating term
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Test if mass balance is positive (i.e. above the ELA)
      IF (MB % values (MB % perm(n)) .ge. 0.0) THEN
        if ( Seasonality ) then
          ! Calculate surace melt in meters of snow equivalent
          Melt = f_dd * SUM(PDD(doy_i:doy_ip1)) * 1.0 ![m] <= [m K-1 d-1] * [K d]
        else
          ! since no seasonality sum up all the melt for the year
          Melt = f_dd * SUM(PDD)    ![m] <= [m K-1 d-1] * [K d]
          ! then partition the melt equally (via fraction of year)
          Melt = Melt * (dt/1.0)    ![m] <= [m] * [a] / [a]
        end if
      ElSE
        ! below the ELA so no latent heat source available
        Q_lat = 0.0
        ! for our purposes here, melt is zero below the ELA
        Melt  = 0.0
      ENDIF

      ! export the melt within the timestep [m]
      surf_melt % Values (surf_melt % perm(n)) = Melt
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Convert surface air temperature to enthalpy
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Calculate mean annual airtemp over timestep and convert from [C] to [K]
      if (Seasonality) then
        ! average airtemp over timestep
        T_surf = SUM(T(doy_i:doy_ip1))/(doy_ip1-doy_i) + 273.15
      else
        ! b/c no seasonality just return mean annual air temp
        T_surf = SUM(T)/365.0 + 273.15
      end if

      ! print info for the headwall, just as a debugging measure
      if (n .eq. N_n) write(*,*) doy_i, doy_ip1, doy_ip1-doy_i, T_surf-273.15, N_v

      ! Temperature can't exced the melting point at the surface
      if (T_surf > 273.15_dp) then
        T_surf = 273.15_dp
      endif

      ! Convert surface air temp to enthalpy
      H_surf = (CapA/2.0_dp*(T_surf**2.0 - T_ref**2.0) + CapB*(T_surf-T_ref)) ! [J kg-1]
      ! Set the surface enthalpy for the Dirichlet Surface B.C.
      Surf_Enthalpy % values (Surf_Enthalpy % perm(n)) = H_surf
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Set variable (w/ depth) surface density in the accumulation zone
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! iterate over vertically aligned nodes
      do i=1,N_v

        ! index of ith vertically aligned node
        ! only works with vertically structured mesh
        cont=n-(i-1)*N_s

        ! Check if mass balance of surface node is postive, i.e. in accumulation zone
        if (MB % values (MB % perm(n)) .ge. 0.0) THEN

          ! approximate firn/ice transiton depth [Cuffey and Paterson, pg. 19]
          z_trans = (1.00 / C_firn) * 1.9

          ! EQN (2.2) from Cuffey and Paterson [kg m-3]
          Dens % values ( Dens % perm(cont) ) = &
          rho_i - (rho_i - rho_s) * exp(-C_firn * Depth % values (Depth % perm(cont)))

        else
          ! if below the ELA, set all nodes to density of ice [kg m-3]
          Dens % values ( Dens % perm(cont)) = rho_i
        endif

        ! Make sure density isn't greater than the ice density
        if (Dens % values ( Dens % perm(cont)) > rho_i) then
          Dens % values ( Dens % perm(cont)) = rho_i
        endif
      end do
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Pseudo percolation scheme for latent heat release
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      ! surface melt available for refreezing [kg m-2]
      pump = Melt * rho_w
      ! start at the bottom of the first layer
      i = 2

      ! if no melt within timestep, then no runoff has occured either [-]
      if (pump .eq. 0.0_dp) RunOff % values ( RunOff % perm(n)) = 0.0_dp

      ! loop over vertical nodes as long as theres still meltwater to refreeze
      do while (pump .ne. 0.0_dp)
        ! index of kth vertical node, only works with vertically structured mesh
        k   = n-(i-1)*N_s
        ! kth layer thicnkess [m]
        dz  = Depth % values (Depth % perm(k-N_s)) - Depth % values (Depth % perm(k))
        ! denisty of current node, as set above
        rho = Dens % values ( Dens % perm(k) )

        ! check if percolation can still occur
        select case (trim(omega_maximum))
          case ("depth")
            ! check if within the firn aquifer
            if (Depth % values (Depth % perm(k)) .le. h_aq) then
              percolate = .true.
            else
              percolate = .false.
            end if

          case ("denisty")

            ! check if nodal density is less than pore close off
            if (rho .le. rho_f) then
              percolate = .true.
            else
              percolate = .false.
            end if
        end select

        if (percolate) then
          ! current water content [-]
          water = (Enth % values (Enth % perm(k)) - H_f % values (H_f % perm(k))) / L_heat

          ! residual water content [kg m-3], i.e. water content left to be filled:
          ! based on:
          !     1. density i.e. porosity
          !     2. difference b/w current and maximum water content
          w_res = (w_max_aq-water) * (1.0_dp - rho/rho_i) * rho_w

          ! residual water content [kg m-3] is bounded by zero, enforce that bound
          w_res = max(w_res, 0.0_dp)

          ! Upper limit for heat source [J m-3 yr-1] based on residual water content
          Q_max = w_res * L_heat * (1.0_dp/dt)
          ! assume all available meltwater refreezes [J m-3 yr-1]
          Q_lat = (pump / (dz * dt)) * L_heat

          if ( Q_lat .gt. Q_max ) then
            ! more meltwater was refrozen than the gridcell can accomodate

            ! amount of excess surface melt still available [kg m-2]
            pump = (Q_lat - Q_max) * dz * dt * (1.0_dp/L_heat)

            ! set the volumetric heat source variable [J m-3 yr-1]
            Q_lat_vol % values ( Q_lat_vol % perm(k)) = Q_max
            ! Note: if the gridcell is fully staurated (i.e. omega == Sr) then
            !       both w_res and Q_max are zero, i.e. no heat is added, but the
            !       melt can still percolate to the next cell
          else
            ! remaining meltwater has been refrozen within the current timestep
            pump = 0.0_dp
            ! set the volumetric heat source variable [J m-3 yr-1]
            Q_lat_vol % values ( Q_lat_vol % perm(k)) = Q_lat
            ! All melt has been consumed, so no-runoff occurs [-]
            RunOff % values ( RunOff % perm(n)) = 0.0_dp
          end if

        else
          ! refreezing can no longer occur. either have reached the bottom of the
          ! firn aquifer or the pore close off density

          ! whatever melt is left is instaneously runoff.
          ! runoff value is specified at the free surface (n) for convinence
          RunOff % values ( RunOff % perm(n)) = pump / (Melt * rho_w)

          ! Set the pump to zero, to stop while loop
          pump = 0.0_dp
        endif

        ! procede to the next vertical layer
        i = i + 1

        ! We've reached the bottom of the vertical column, cut the itteration
        if (i==N_v) pump=0.0_dp
      end do !percolation while loop

      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! Set a seasonal snow layer, a crude approximation
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      ! if mean air temp less than 0 C set surface denisty to that of snow
      if ( T_surf .le. 273.15_dp ) then
        Dens % values ( Dens % perm(n)) = rho_s
      end if
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ENDIF ! if at free sufrace
  ENDDO ! loop over nodes

  contains

  subroutine temp_std(std, coefs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time dependent standard deviation in air temperature, which facilitates
    ! greater varaince in winter than summer, but is also neccisary for
    ! recreating melt curve
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    real(dp), intent(in)                  :: coefs(3) ! Coefs for T_std(doy)
    real(dp), intent(out), dimension(365) :: std      ! std(doy) array
    ! internal parameters
    real(dp), dimension(365)              :: d        ! Day of year array  [DOY]
    integer                               :: i
    ! populate day of year array
    d = [ (real(i, dp), i = 1, 365) ]
    ! evaluate the quadaratic function for the daily std
    std = coefs(1) * d**2  + coefs(2) * d + coefs(3)
  end subroutine temp_std

  ! subfunction for reading constants and error checking
  function GetParam(Model, constant_name) result(constant)
    real :: constant
    logical :: GotIt
    TYPE(Model_t) :: Model
    character(len=*) :: constant_name

    constant = GetConstReal(Model % Constants, trim(constant_name), GotIt)

    if (.not. GotIt) then
      call fatal(Caller//'---> GetParam', &
                 'Could not find '//trim(constant_name) )
    end if

  end function GetParam
END SUBROUTINE Surface_Processes
