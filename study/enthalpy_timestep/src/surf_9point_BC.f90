FUNCTION SE_discrete(Model, Node, InputArray) RESULT(Enthalpy)
  USE DefUtils
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
                   dt             !
  REAL(KIND=dp) :: T_surf         ! intermediate result
  REAL(KIND=dp) :: z              ! surface elevation of current node  [m a.s.l.]
  REAl(KIND=dp) :: alpha          ! annual surf. air temp. amplitude   [K]
  REAL(KIND=dp) :: T_mean         ! mean annual surf. air. temp        [K]
  REAL(KIND=dp) :: ref_z          ! reference surface elevation        [m a.s.l.]
  REAL(KIND=dp) :: grad_T         ! air temp lapse rate                [K m^{-1}]
  REAL(KIND=dp) :: T_ref, CapA, CapB
  REAL(KIND=dp) :: year


  integer :: DOY                  ! number of timesteps in a year
  logical :: Found, Transient

  Transient = GetString(GetSimulation(), "Simulation type", Found)=='transient'

  ! Enthalpy related constants
  T_ref     = GetConstReal(Model % Constants, "t_ref_enthalpy")            ! [K]
  CapA      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A ") ! [J kg-1 K-2]
  CapB      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")  ! [J kg-1 K-1]
  alpha     = GetConstReal(Model % Constants, "delta_T")                   ! [K]
  grad_T    = GetConstReal(Model % Constants, "grad_T" )                   ! [K m^{-1}]
  T_mean    = GetConstReal(Model % Constants, "T_mean" )                   ! [K]
  ref_z     = GetConstReal(Model % Constants, "ref_z" )                    ! [m a.s.l.]
  z         = InputArray(1)              ! elevation of current surface node [m]


  if (Transient) then
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    Time     =  TimeVar % Values(1)
    ! Convert real time into approximate DOY
    DOY    = NINT((Time - floor(Time))*365.0)
    ! Find surface temp for DOY
    T_surf=alpha*sin(2*3.14*DOY/365)+grad_T*(ref_z-z)+T_mean
  else
    ! calculate mean annual surface temperature
    T_surf = grad_T*(ref_z-z)+T_mean
  endif

  ! Temperature can't exced the melting point at the surface
  if (T_surf > 273.15) then
    T_surf = 273.15
  endif

  Enthalpy  = (CapA/2.0*(T_surf**2 - T_ref**2) + CapB*(T_surf-T_ref)) ! [J kg-1]

  RETURN

END FUNCTION SE_discrete

FUNCTION SE_mean(Model, Node, InputArray) RESULT(Enthalpy)
  USE DefUtils
  IMPLICIT NONE

  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model         ! the access point to everything about the model
  INTEGER       :: Node          ! the current Node number
  REAL(KIND=dp) :: InputArray(1) ! Contains the argument passed to the function
  REAL(KIND=dp) :: Enthalpy      ! the final result
  TYPE(Variable_t), POINTER :: TimeVar
  REAL(KIND=dp), POINTER :: TimestepSizes(:,:)

  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: Time_n,      & ! simulation time                    [a]
                   Time_nm1   , & ! Time_n minus 1                     [a]
                   dt_sizes(1), & ! vector of timesteps, but assumed of size one
                   dt             ! timestep size
  REAL(KIND=dp) :: T_surf         ! intermediate result
  REAL(KIND=dp) :: z              ! surface elevation of current node  [m a.s.l.]
  REAl(KIND=dp) :: alpha          ! annual surf. air temp. amplitude   [K]
  REAL(KIND=dp) :: T_mean         ! mean annual surf. air. temp        [K]
  REAL(KIND=dp) :: ref_z          ! reference surface elevation        [m a.s.l.]
  REAL(KIND=dp) :: grad_T         ! air temp lapse rate                [K m^{-1}]
  REAL(KIND=dp) :: T_ref, CapA, CapB
  REAL(KIND=dp) :: year, sum


  integer :: i, DOY_n, DOY_nm1       ! number of timesteps in a year
  logical :: Found, Transient

  Transient = GetString(GetSimulation(), "Simulation type", Found)=='transient'

  ! Enthalpy related constants
  T_ref     = GetConstReal(Model % Constants, "t_ref_enthalpy")            ! [K]
  CapA      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A ") ! [J kg-1 K-2]
  CapB      = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")  ! [J kg-1 K-1]
  alpha     = GetConstReal(Model % Constants, "delta_T")                   ! [K]
  grad_T    = GetConstReal(Model % Constants, "grad_T" )                   ! [K m^{-1}]
  T_mean    = GetConstReal(Model % Constants, "T_mean" )                   ! [K]
  ref_z     = GetConstReal(Model % Constants, "ref_z" )                    ! [m a.s.l.]
  z         = InputArray(1)              ! elevation of current surface node [m]


  if (Transient) then
    TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
    Time_n   =  TimeVar % Values(1)
    TimestepSizes => ListGetConstRealArray(Model%Simulation, "Timestep Sizes", Found)
    dt       = TimestepSizes(1,1)

    !write(*,*) dt
    Time_nm1 = Time_n - dt

    ! Convert real time into approximate DOY
    DOY_n    = NINT((Time_n - floor(Time_n))*365.0)
    DOY_nm1  = NINT((Time_nm1 - floor(Time_nm1))*365.0)

    ! Find the average surface temp, between current and last timesteps
    do i = DOY_nm1, DOY_n, 1
      sum = sum + alpha*sin(2*3.14*i/365)+grad_T*(ref_z-z)+T_mean
    end do
    ! Find surface temp for DOY
    T_surf= sum / real(DOY_n - DOY_nm1)
  else
    ! calculate mean annual surface temperature
    T_surf = grad_T*(ref_z-z)+T_mean
  endif

  ! Temperature can't exced the melting point at the surface
  if (T_surf > 273.15) then
    T_surf = 273.15
  endif

  Enthalpy  = (CapA/2.0*(T_surf**2 - T_ref**2) + CapB*(T_surf-T_ref)) ! [J kg-1]

  RETURN

END FUNCTION SE_mean
