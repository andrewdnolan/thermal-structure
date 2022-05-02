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
      call SurfTemp(z, T, alpha, dTdz, &
                    z_ref, T_mean, T_peak, (/ std_c0, std_c1, std_c2 /))

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

! function GetModelConstant(constant_name, GotIt)
!   constant = GetConstReal(Model % Constants, constant_name, GotIt)
!   if (.not. GotIt) then
!     call fatal('SurfaceMassBalance ---> GetModelConstant', 'Could not find ')
!   end if
! end function GetModelConstant

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
  logical       :: Found,       &
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

  Transient = GetString(GetSimulation(), "Simulation type", Found)=='transient'

  ! Enthalpy related constants
  T_ref     = GetConstReal(Model % Constants, "t_ref_enthalpy")            ! [K]
  CapA      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A")  ! [J kg-1 K-2]
  CapB      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")  ! [J kg-1 K-1]
  ! Air temperature related constants
  alpha     = GetConstReal(Model % Constants, "alpha")    ! [K]
  dTdz      = GetConstReal(Model % Constants, "dTdz" )    ! [K m^{-1}]
  T_mean    = GetConstReal(Model % Constants, "T_mean")  ! [K]
  T_peak    = INT(ANINT(GetConstReal( Model % Constants, "T_peak" )))  ! [DOY]
  z_ref     = GetConstReal(Model % Constants, "z_ref" )   ! [m a.s.l.]
  std_c0    = GetConstReal(Model % Constants, "std_c0")   ! [K?]
  std_c1    = GetConstReal(Model % Constants, "std_c1")   ! [K?]
  std_c2    = GetConstReal(Model % Constants, "std_c2")   ! [K?]

  ! unpack surface elevation of current surface node [m]
  z = InputArray(1)

  ! Calculate air temperature curve
  call SurfTemp(z, T, alpha, dTdz, &
                z_ref, T_mean, T_peak, (/ std_c0, std_c1, std_c2 /))

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
END FUNCTION getSurfaceEnthalpy

SUBROUTINE SurfaceMelt( Model, Solver, dt, TransientSimulation )
  !-----------------------------------------------------------------------------
  !> Implementation of surface heating ($Q_m$) scheme (Eqn. 9) from
  !> Wilson and Flowers (2013).
  !-----------------------------------------------------------------------------
  USE DefUtils       ! provides you with most Elmer functionality
  USE SolverUtils
  USE ElementUtils

  ! Local modules
  USE SurfaceTemperature

  IMPLICIT NONE      ! saves you from stupid errors

  !----------------------------------------------------------------------------
  ! external variables
  !----------------------------------------------------------------------------
  LOGICAL                   :: TransientSimulation  ! Should be .FAlSE. always for now
  REAL(KIND=dp)             :: dt                   ! Should be 0 for now
  TYPE(Model_t)             :: Model
  TYPE(Solver_t)            :: Solver
  TYPE(Element_t),  POINTER :: Element

  TYPE(Variable_t), POINTER :: Melt,    &
                               Depth,   &
                               TimeVar, &
                               MB

  !----------------------------------------------------------------------------
  ! Local variables
  !----------------------------------------------------------------------------
  ! air temp related
  INTEGER       :: T_peak         ! DOY of annual temp peak      [DOY]
  REAL(KIND=dp) :: z_ref,       & ! reference surface elevation  [m a.s.l.]
                   alpha,       & ! Anual air temp. amp          [K]
                   dTdz,        & ! air temp lapse rate          [K m-1]
                   T_mean,      & ! mean annual surf. air. temp  [K]
                   T(365),      & ! surface air temperature      [K]
                   std_c0,      &
                   std_c1,      &
                   std_c2
  ! internal integers
  INTEGER       :: i, n,        &  ! index counter
                   cont,        &  ! vertically alligned node index counter
                   d,           &  ! julian calendar day  counter
                   N_n,         &  ! number of model    nodes
                   N_s,         &  ! number of surface  nodes
                   N_v,         &  ! number of vertical nodes
                   doy_i,       &
                   doy_ip1

  REAL(KIND=dp) :: z,           &
                   Time,        &
                   T_melt,      &
                   f_snow,      &  ! Degree-day factor  [m K^-1 a^-1]
                   rho_i,       &  ! denisty of ice     [Kg m^-3]
                   PDD(365) = 0    ! + degrees. for DOY [K]

  LOGICAL :: GotIt,first_time=.true.

  !-------------------------
  ! First time loop
  !-------------------------
  save first_time,N_v,N_s,N_n


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

  ! Pointer to the model variables
  Melt   => VariableGet( Model % Variables, "Melting")
  Depth  => VariableGet( Model % Variables, "Depth")
  MB     => VariableGet( Model % Variables, "mass balance")

  ! Physical params
  rho_i  = GetConstReal(Model % Constants, "rho_i")    ![kg m-3]
  ! Melt params
  f_snow = GetConstReal(Model % Constants, "f_snow")   ![kg m-2 yr-1]
  T_melt = GetConstReal(Model % Constants, "T_melt")   ![K]
  ! Air temperature related constants
  alpha  = GetConstReal(Model % Constants, "alpha")    ! [K]
  dTdz   = GetConstReal(Model % Constants, "dTdz" )    ! [K m^{-1}]
  T_mean = GetConstReal(Model % Constants, "T_mean")   ! [K]
  T_peak = INT(ANINT(GetConstReal( Model % Constants, "T_peak" )))  ! [DOY]
  z_ref  = GetConstReal(Model % Constants, "z_ref" )   ! [m a.s.l.]
  std_c0 = GetConstReal(Model % Constants, "std_c0")   ! [K?]
  std_c1 = GetConstReal(Model % Constants, "std_c1")   ! [K?]
  std_c2 = GetConstReal(Model % Constants, "std_c2")   ! [K?]


  if (TransientSimulation) then
    ! if transient get current timestep
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    ! Get current time
    Time     =  TimeVar % Values(1)
    ! Find DOY of current timestep
    doy_i   = NINT((Time - floor(Time)) * 365.0)
    ! find DOY of next    timestep
    doy_ip1 = NINT(((Time+dt) - floor(Time)) * 365.0)

    write(*,*) doy_ip1-doy_i
  else
    ! if steady state get yearly amount of melt
    doy_i   = 1
    doy_ip1 = 365
  endif

  ! Outter Most Loop: Itterate of model nodes
  DO n=1,N_n
    ! Check if depth == 0.0 for node, i.e. is it a surface node
    IF (Depth%Values(Depth%perm(n))==0.0) THEN

      z = model%nodes%y(n)

      ! Calculate nodal air temperature curve
      call SurfTemp(z, T, alpha, dTdz, &
                    z_ref, T_mean, T_peak, (/ std_c0, std_c1, std_c2 /))

      ! only loop over DOY within current timestep
      DO d=doy_i,doy_ip1
        ! Subtract the melting temperature to find the daily positive degrees
        PDD(d) = MAX(T(d)-T_melt, 0.0) ! [K]
      END DO

      ! Test if 2nd coord (z) is greater then Z_ELA
      IF (MB % values (MB % perm(n)) > 0.0) THEN
        ! Calculate the amount of melt
        Melt % values (Melt % perm(n)) = f_snow*SUM(PDD) * (1/rho_i) ! [m yr-1]
      ElSE
        ! Set the firn thickness and melt to 0 everywhere else
        Melt % values (Melt % perm(n)) = 0.0
      ENDIF
    ENDIF
  ENDDO
END SUBROUTINE SurfaceMelt

SUBROUTINE UpdateFirn( Model, Solver, dt, TransientSimulation )
  !------------------------------------------------------------------------------
  USE DefUtils       ! provides you with most Elmer functionality
  USE SolverUtils    ! not totally sure about this, but adrian does it so....
  USE ElementUtils   ! not totally sure about this, but adrian does it so....
  !------------------------------------------------------------------------------


  IMPLICIT NONE      ! saves you from stupid errors

  !----------------------------------------------------------------------------
  ! external variables
  !----------------------------------------------------------------------------
  LOGICAL                   :: TransientSimulation  ! Should be .FAlSE. always for now
  REAL(KIND=dp)             :: dt                   ! Should be 0 for now
  TYPE(Model_t)             :: Model
  TYPE(Solver_t)            :: Solver
  TYPE(Element_t),  POINTER :: Element

  TYPE(Variable_t), POINTER :: Firn,    &           ! firn thickness        [m]
                               Depth,   &
                               Height,  &
                               TimeVar, &
                               MB,      &
                               Dens   ! Ice Depth       [m]

  !----------------------------------------------------------------------------
  ! Local variables
  !----------------------------------------------------------------------------
  ! internal integers
  INTEGER       :: i, n,        &  ! index counter
                   cont,        &  ! vertically alligned node index counter
                   d,           &  ! julian calendar day  counter
                   N_n,         &  ! number of model    nodes
                   N_s,         &  ! number of surface  nodes
                   N_v,         &  ! number of vertical nodes
                   doy_i,       &
                   doy_ip1
  real(kind=dp) :: test

  REAL(KIND=dp) :: z,           &
                   Time,        &
                   T_melt,      &
                   f_snow,      &  ! Degree-day factor  [m K^-1 a^-1]
                   rho_i,       &  ! denisty of ice     [Kg m^-3]
                   rho_s,       &  ! denisty of water   [Kg m^-3]
                   PDD(365) = 0    ! + degrees. for DOY [K]

  LOGICAL :: GotIt,first_time=.true.

  !-------------------------
  ! First time loop
  !-------------------------
  save first_time,N_v,N_s,N_n


  ! Variable declaration that only needs to be run once (i.e. first time)
  IF (first_time) THEN
    ! set first_time to false for all subsequent time steps
    first_time=.false.


    ! Loop over all model nodes
    DO i=1,model % NumberOfNodes
      IF (model%nodes%x(i+1) < model%nodes%x(i)) THEN
        EXIT
      ENDIF

    N_n = Model % NumberOfNodes ! Number of Nodes in Models
    N_s = i                     ! Number of surface nodes
    N_v = N_n/i                 ! Number of vertical nodes
    ENDDO
  ENDIF
  ! End "first_time" loop

  ! Pointer to the model variables
  Firn   => VariableGet( Model % Variables, "Firn")
  Height => VariableGet( Model % Variables, "Height")
  Depth  => VariableGet( Model % Variables, "Depth")
  Dens   => VariableGet( Model % Variables, "Densi")
  MB     => VariableGet( Model % Variables, "mass balance")

  ! Physical params
  rho_i  = GetConstReal(Model % Constants, "rho_i")    ![kg m-3]
  rho_s  = GetConstReal(Model % Constants, "rho_s")    ![kg m-3]


  ! Outter Most Loop: Itterate of model nodes
  DO n=1,N_n

    ! Check if depth == 0.0 for node, i.e. is it a surface node
    IF (Depth%Values(Depth%perm(n))==0.0) THEN

      z = model%nodes%y(n)

      ! Test if 2nd coord (z) is greater then Z_ELA
      IF (MB % values (MB % perm(n)) > 0.0) THEN
        ! Set the firn thickness
        Firn % values (Firn % perm(n)) = MB % values (MB % perm(n)) /  0.05! [m]
        !Firn % values (Firn % perm(n)) = 50 ! [m]

        ! If firn thickness exceeds ice-thickness, then set firn thickness to
        ! ice thickness
        IF (Firn%values(Firn%perm(n)) > Height%values(Height%perm(n))) THEN
          Firn % values (Firn % perm(n)) = Height%values(Height%perm(n))
        ENDIF

      ElSE
        ! Set the firn thickness and melt to 0 everywhere else
        Firn % values (Firn % perm(n)) = 0.0
      ENDIF

      ! iterate over vertically aligned nodes to set firn density
      do i=1,N_v

        ! index of ith vertically aligned node
        cont=n-(i-1)*N_s

        ! Check if within firn aquifer, if so set linear density profile
        if (Firn % values (Firn % perm(n)) > 1.0) then

          ! Cuffey and paterson EQN 2.2
          Dens % values ( Dens % perm(cont)) = &
          rho_i - (rho_i - rho_s) * exp( -1.0 * Depth % values (Depth % perm(cont)) / &
                                         (Firn % values (Firn  % perm(n)) / 2.1))

        else
          Dens % values ( Dens % perm(cont)) = rho_i
        endif

        ! Make sure density isn't greater than the ice density
        if (Dens % values ( Dens % perm(cont)) > rho_i) then
          Dens % values ( Dens % perm(cont)) = rho_i
        endif

      end do
    ENDIF
  ENDDO
END SUBROUTINE UpdateFirn


SUBROUTINE Wilson_and_Flowers_2013( Model, Solver, dt, TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils       ! provides you with most Elmer functionality
  USE SolverUtils    ! not totally sure about this, but adrian does it so....
  USE ElementUtils   ! not totally sure about this, but adrian does it so....

  ! Local modules
  USE SurfaceTemperature

  IMPLICIT NONE      ! saves you from stupid errors

  !----------------------------------------------------------------------------
  ! external variables
  !----------------------------------------------------------------------------
  LOGICAL                   :: TransientSimulation  ! Should be .FAlSE. always for now
  REAL(KIND=dp)             :: dt                   ! Should be 0 for now
  TYPE(Model_t)             :: Model
  TYPE(Solver_t)            :: Solver
  TYPE(Element_t),  POINTER :: Element

  TYPE(Variable_t), POINTER :: Q_lat,    &
                               Firn,     &  ! firn thickness        [m]
                               Depth,    &
                               TimeVar,  &
                               MB,       &
                               Dens         ! Ice Depth       [m]

  !----------------------------------------------------------------------------
  ! Local variables
  !----------------------------------------------------------------------------
  ! air temp related
  INTEGER       :: T_peak         ! DOY of annual temp peak      [DOY]
  REAL(KIND=dp) :: z_ref,       & ! reference surface elevation  [m a.s.l.]
                   alpha,       & ! Anual air temp. amp          [K]
                   dTdz,        & ! air temp lapse rate          [K m-1]
                   T_mean,      & ! mean annual surf. air. temp  [K]
                   T(365),      & ! surface air temperature      [K]
                   std_c0,      &
                   std_c1,      &
                   std_c2
  ! internal integers
  INTEGER       :: i, n,        &  ! index counter
                   cont,        &  ! vertically alligned node index counter
                   d,           &  ! julian calendar day  counter
                   N_n,         &  ! number of model    nodes
                   N_s,         &  ! number of surface  nodes
                   N_v,         &  ! number of vertical nodes
                   doy_i,       &
                   doy_ip1
  real(kind=dp) :: test

  REAL(KIND=dp) :: z,           &
                   Time,        &
                   T_melt,      &
                   L_heat,      &  ! Latnet heat of fusion [J kg-1]
                   f_snow,      &  ! Degree-day factor     [m K^-1 a^-1]
                   rho_i,       &  ! denisty of ice        [Kg m^-3]
                   rho_w,       &  ! denisty of water      [Kg m^-3]
                   rho_s,       &  ! denisty of snow       [Kg m^-3]
                   PDD(365) = 0    ! + degrees. for DOY    [K]

  LOGICAL :: GotIt,first_time=.true.

  !-------------------------
  ! First time loop
  !-------------------------
  save first_time,N_v,N_s,N_n


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

  ! Pointer to the model variables
  Q_lat    => VariableGet( Model % Variables, "Q_lat")
  Firn     => VariableGet( Model % Variables, "Firn")
  Depth    => VariableGet( Model % Variables, "Depth")
  Dens     => VariableGet( Model % Variables, "Densi")
  MB       => VariableGet( Model % Variables, "mass balance")

  ! Physical Parameters
  rho_i  = GetConstReal(Model % Constants,  "rho_i")   ![kg m-3]
  rho_s  = GetConstReal(Model % Constants,  "rho_s")   ![kg m-3]
  rho_w  = GetConstReal(Model % Constants,  "rho_w")   ![kg m-3]
  ! Physical constants
  L_heat = GetConstReal( Model % Constants, "L_heat")  ![J kg-1]
  ! Melt Parameters
  f_snow = GetConstReal(Model % Constants, "f_snow")   ![kg m-2 yr-1]
  T_melt = GetConstReal(Model % Constants, "T_melt")   ![K]
  ! Air temperature related Parameters
  alpha  = GetConstReal(Model % Constants, "alpha")    ! [K]
  dTdz   = GetConstReal(Model % Constants, "dTdz" )    ! [K m^{-1}]
  T_mean = GetConstReal(Model % Constants, "T_mean")   ! [K]
  T_peak = INT(ANINT(GetConstReal( Model % Constants, "T_peak" )))  ! [DOY]
  z_ref  = GetConstReal(Model % Constants, "z_ref" )   ! [m a.s.l.]
  std_c0 = GetConstReal(Model % Constants, "std_c0")   ! [K?]
  std_c1 = GetConstReal(Model % Constants, "std_c1")   ! [K?]
  std_c2 = GetConstReal(Model % Constants, "std_c2")   ! [K?]



  if (TransientSimulation) then
    ! if transient get current timestep
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    ! Get current time
    Time    =   TimeVar % Values(1)
    ! Find DOY of current timestep
    doy_i   = NINT(MOD(time, 1.0) * 365.0)

    ! B/C round off error ith doy sometimes is 365 instead of 0. Check and fix
    if (doy_i==365) then
      doy_i = 0
    end if

    ! find DOY of next    timestep
    doy_ip1 = NINT(MOD(time+dt, 1.0) * 365.0)
    ! write(*,*) Time,  Time+dt, dt
    ! write(*,*) doy_i, doy_ip1, doy_ip1-doy_i
    
  else
    ! set the "timestep" as one year
    dt = 1
    ! if steady state get yearly amount of melt
    doy_i   = 1
    doy_ip1 = 365
  endif

  ! Outter Most Loop: Itterate of model nodes
  DO n=1,N_n
    ! Check if depth == 0.0 for node, i.e. is it a surface node
    IF (Depth%Values(Depth%perm(n))==0.0) THEN

      z = model%nodes%y(n)

      ! Calculate nodal air temperature curve
      call SurfTemp(z, T, alpha, dTdz, &
                    z_ref, T_mean, T_peak, (/ std_c0, std_c1, std_c2 /))

      ! only loop over DOY within current timestep
      DO d=doy_i,doy_ip1
        ! Subtract the melting temperature to find the daily positive degrees
        PDD(d) = MAX(T(d)-T_melt, 0.0) ! [K]
      END DO

      ! Test if 2nd coord (z) is greater then Z_ELA
      IF (MB % values (MB % perm(n)) > 0.0) THEN
        ! Set the firn thickness
        Firn % values (Firn % perm(n)) = 3.0 ! [m]

        ! set the heat source one layer below surface so as not to conflict with
        ! direchlet boundary condition
        cont=n-N_s
        ! Calculate the volumtric latent heat source [J m-3 yr-1]
        Q_lat % values (Q_lat % perm(cont)) = (1 - 0.3) * &
                                           (rho_w/3.0) * L_heat * (1/rho_s) * f_snow*SUM(PDD)
      ElSE
        ! Set the firn thickness and melt to 0 everywhere else
        Firn  % values (Firn  % perm(n)) = 0.0
        Q_lat % values (Q_lat % perm(n)) = 0.0
      ENDIF

      ! iterate over vertically aligned nodes to set firn density
      do i=1,N_v

        ! index of ith vertically aligned node
        cont=n-(i-1)*N_s

        ! Check if within firn aquifer, if so set linear density profile
        if (Firn % values (Firn % perm(n)) > 1.0) then

          ! Cuffey and paterson EQN 2.2
          Dens % values ( Dens % perm(cont)) = &
          rho_i - (rho_i - rho_s) * exp( -0.05 * Depth % values (Depth % perm(cont)))

        else
          Dens % values ( Dens % perm(cont)) = rho_i
        endif

        ! Make sure density isn't greater than the ice density
        if (Dens % values ( Dens % perm(cont)) > rho_i) then
          Dens % values ( Dens % perm(cont)) = rho_i
        endif
      end do

      ! if no melt occurs set the surface denisty to that of snow
      if ( SUM(PDD) == 0 ) then
        Dens % values ( Dens % perm(n)) = rho_s
      end if

    ENDIF
  ENDDO
END SUBROUTINE Wilson_and_Flowers_2013


! SUBROUTINE Wilson_and_Flowers_2013( Model, Solver, dt, TransientSimulation )
! !------------------------------------------------------------------------------
!   USE DefUtils       ! provides you with most Elmer functionality
!   USE SolverUtils    ! not totally sure about this, but adrian does it so....
!   USE ElementUtils   ! not totally sure about this, but adrian does it so....
!
!   ! Local modules
!   USE SurfaceTemperature
!
!   IMPLICIT NONE      ! saves you from stupid errors
!
!   !----------------------------------------------------------------------------
!   ! external variables
!   !----------------------------------------------------------------------------
!   LOGICAL                   :: TransientSimulation  ! Should be .FAlSE. always for now
!   REAL(KIND=dp)             :: dt                   ! Should be 0 for now
!   TYPE(Model_t)             :: Model
!   TYPE(Solver_t)            :: Solver
!   TYPE(Element_t),  POINTER :: Element
!
!   TYPE(Variable_t), POINTER :: SurfEnth, &
!                                H_f,      &  ! Enthalpy of fusion (phase change enth)
!                                Firn,     &  ! firn thickness        [m]
!                                Depth,    &
!                                TimeVar,  &
!                                MB,       &
!                                Dens         ! Ice Depth       [m]
!
!   !----------------------------------------------------------------------------
!   ! Local variables
!   !----------------------------------------------------------------------------
!   ! air temp related
!   INTEGER       :: T_peak         ! DOY of annual temp peak      [DOY]
!   REAL(KIND=dp) :: z_ref,       & ! reference surface elevation  [m a.s.l.]
!                    alpha,       & ! Anual air temp. amp          [K]
!                    dTdz,        & ! air temp lapse rate          [K m-1]
!                    T_mean,      & ! mean annual surf. air. temp  [K]
!                    T(365),      & ! surface air temperature      [K]
!                    std_c0,      &
!                    std_c1,      &
!                    std_c2
!   ! internal integers
!   INTEGER       :: i, n,        &  ! index counter
!                    cont,        &  ! vertically alligned node index counter
!                    d,           &  ! julian calendar day  counter
!                    N_n,         &  ! number of model    nodes
!                    N_s,         &  ! number of surface  nodes
!                    N_v,         &  ! number of vertical nodes
!                    doy_i,       &
!                    doy_ip1
!   real(kind=dp) :: test
!
!   REAL(KIND=dp) :: z,           &
!                    Time,        &
!                    T_melt,      &
!                    L_heat,      &  ! Latnet heat of fusion [J kg-1]
!                    f_snow,      &  ! Degree-day factor     [m K^-1 a^-1]
!                    rho_i,       &  ! denisty of ice        [Kg m^-3]
!                    rho_w,       &  ! denisty of water      [Kg m^-3]
!                    rho_s,       &  ! denisty of snow       [Kg m^-3]
!                    PDD(365)=0,  &  ! + degrees. for DOY    [K]
!                    T_ref,       &
!                    CapA,        &
!                    CapB,        &
!                    H_surf,      &  ! surface enthalpy
!                    Q_m             ! latent heat souce
!
!   LOGICAL :: GotIt,first_time=.true.
!
!   !-------------------------
!   ! First time loop
!   !-------------------------
!   save first_time,N_v,N_s,N_n
!
!
!   ! Variable declaration that only needs to be run once (i.e. first time)
!   IF (first_time) THEN
!     ! set first_time to false for all subsequent time steps
!     first_time=.false.
!
!
!     ! Loop over all model nodes
!     DO i=1,model % NumberOfNodes
!       IF (model%nodes%x(i+1) < model%nodes%x(i)) THEN
!         EXIT
!       ENDIF
!     ENDDO
!     N_n = Model % NumberOfNodes ! Number of Nodes in Models
!     N_s = i                     ! Number of surface nodes
!     N_v = N_n/i                 ! Number of vertical nodes
!   ENDIF
!   ! End "first_time" loop
!
!   ! Pointer to the model variables
!   SurfEnth => VariableGet( Model % Variables, "Surface Enthalpy")
!   H_f      => VariableGet( Model % Variables, "Phase Change Enthalpy")
!   Firn     => VariableGet( Model % Variables, "Firn")
!   Depth    => VariableGet( Model % Variables, "Depth")
!   Dens     => VariableGet( Model % Variables, "Densi")
!   MB       => VariableGet( Model % Variables, "mass balance")
!
!   ! Physical Parameters
!   rho_i  = GetConstReal(Model % Constants,  "rho_i")   ![kg m-3]
!   rho_s  = GetConstReal(Model % Constants,  "rho_s")   ![kg m-3]
!   rho_w  = GetConstReal(Model % Constants,  "rho_w")   ![kg m-3]
!   ! Physical constants
!   L_heat = GetConstReal( Model % Constants, "L_heat")  ![J kg-1]
!   T_ref  = GetConstReal(Model % Constants, "t_ref_enthalpy")            ! [K]
!   CapA   = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A")  ! [J kg-1 K-2]
!   CapB   = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")  ! [J kg-1 K-1]
!   ! Melt Parameters
!   f_snow = GetConstReal(Model % Constants, "f_snow")   ![kg m-2 yr-1]
!   T_melt = GetConstReal(Model % Constants, "T_melt")   ![K]
!   ! Air temperature related Parameters
!   alpha  = GetConstReal(Model % Constants, "alpha")    ! [K]
!   dTdz   = GetConstReal(Model % Constants, "dTdz" )    ! [K m^{-1}]
!   T_mean = GetConstReal(Model % Constants, "T_mean")   ! [K]
!   T_peak = INT(ANINT(GetConstReal( Model % Constants, "T_peak" )))  ! [DOY]
!   z_ref  = GetConstReal(Model % Constants, "z_ref" )   ! [m a.s.l.]
!   std_c0 = GetConstReal(Model % Constants, "std_c0")   ! [K?]
!   std_c1 = GetConstReal(Model % Constants, "std_c1")   ! [K?]
!   std_c2 = GetConstReal(Model % Constants, "std_c2")   ! [K?]
!
!
!
!   if (TransientSimulation) then
!     ! if transient get current timestep
!     TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
!     ! Get current time
!     Time     =  TimeVar % Values(1)
!     ! Find DOY of current timestep
!     doy_i   = NINT((Time - floor(Time)) * 365.0)
!     ! find DOY of next    timestep
!     doy_ip1 = NINT(((Time+dt) - floor(Time)) * 365.0)
!
!     write(*,*) doy_ip1-doy_i
!   else
!     ! set the "timestep" as one year
!     dt = 1
!     ! if steady state get yearly amount of melt
!     doy_i   = 1
!     doy_ip1 = 365
!   endif
!
!   ! Outter Most Loop: Itterate of model nodes
!   DO n=1,N_n
!     ! Check if depth == 0.0 for node, i.e. is it a surface node
!     IF (Depth%Values(Depth%perm(n))==0.0) THEN
!
!       z = model%nodes%y(n)
!
!       ! Calculate nodal air temperature curve
!       call SurfTemp(z, T, alpha, dTdz, &
!                     z_ref, T_mean, T_peak, (/ std_c0, std_c1, std_c2 /))
!
!       ! only loop over DOY within current timestep
!       DO d=doy_i,doy_ip1
!         ! Subtract the melting temperature to find the daily positive degrees
!         PDD(d) = MAX(T(d)-T_melt, 0.0) ! [K]
!       END DO
!
!       ! Calculate mean air temp of timestep [K]
!       T_mean = SUM( T(doy_i:doy_ip1) ) / (doy_ip1-doy_i) + 273.15
!       ! temperature of ice cannot exceed 0 C
!       if (T_mean > 273.15) then
!         T_mean = 273.15
!       end if
!       ! Calculate surface enthalpy from just direchlet condition [J kg-1]
!       H_surf = (CapA/2.0*(T_mean**2 - T_ref**2) + CapB*(T_mean-T_ref))
!
!       ! Test if 2nd coord (z) is greater then Z_ELA
!       IF (MB % values (MB % perm(n)) > 0.0) THEN
!         ! Set the firn thickness
!         Firn % values (Firn % perm(n)) = 3.0 ! [m]
!
!         ! Calculate the volumtric latent heat source [J m-3 yr-1]
!         Q_m = (1 - 0.3)*(rho_w/3.0)*L_heat*(1/rho_s)*f_snow*SUM(PDD)
!         ! Convert from volumetric heat source to enthalpy [J kg-1]
!         Q_m = Q_m * (1/rho_s) * dt
!
!         ! Add the enthalpy released from meltwater refreezing
!         H_surf = H_surf + Q_m
!
!         ! If the surface enthalpy from Dirichlet condition + heat source exceeds
!         ! the enthalpy of fusion, then set surface enthalpy to H_f
!         if (H_surf > H_f % values ( H_f % perm(n) ) ) then
!           H_surf =  H_f % values ( H_f % perm(n) )
!         end if
!
!         ! Set the value for the surface enthalpy
!         SurfEnth % values ( SurfEnth % perm (n) ) = H_surf
!
!         ! ! set the heat source one layer below surface so as not to conflict with
!         ! ! direchlet boundary condition
!         ! cont=n-N_s
!         ! ! Calculate the volumtric latent heat source [J m-3 yr-1]
!         ! Q_lat % values (Q_lat % perm(cont)) = (1 - 0.3) * &
!         !                                    (rho_w/3.0) * L_heat * (1/rho_s) * f_snow*SUM(PDD)
!       ElSE
!         ! Set the firn thickness and melt to 0 everywhere else
!         Firn  % values (Firn  % perm(n)) = 0.0
!         SurfEnth % values (SurfEnth % perm(n)) = H_surf
!       ENDIF
!
!       ! iterate over vertically aligned nodes to set firn density
!       do i=1,N_v
!
!         ! index of ith vertically aligned node
!         cont=n-(i-1)*N_s
!
!         ! Check if within firn aquifer, if so set linear density profile
!         if (Firn % values (Firn % perm(n)) > 1.0) then
!
!           ! Cuffey and paterson EQN 2.2
!           Dens % values ( Dens % perm(cont)) = &
!           rho_i - (rho_i - rho_s) * exp( -0.05 * Depth % values (Depth % perm(cont)))
!
!         else
!           Dens % values ( Dens % perm(cont)) = rho_i
!         endif
!
!         ! Make sure density isn't greater than the ice density
!         if (Dens % values ( Dens % perm(cont)) > rho_i) then
!           Dens % values ( Dens % perm(cont)) = rho_i
!         endif
!       end do
!     ENDIF
!   ENDDO
! END SUBROUTINE Wilson_and_Flowers_2013
