!==============================================================================
function periodic_surge(Model, Node, omega) result(beta)
!==============================================================================
  ! This function returns the time dependent slip coefficent for periodic
  ! surge simulations.

  use DefUtils
  USE SolverUtils
  USE ElementUtils

  implicit none

  ! Solver name
  CHARACTER(*), PARAMETER :: Caller = 'periodic_surge'

  type(Model_t) :: Model
  logical       :: GotIt
  integer       :: Node            ! current node number
  real(kind=dp) :: omega,        & ! water content fraction [--]
                   omega_thresh, & ! w.c. threshold to be considered temperate [--]
                   beta,         & ! nodal slip coefficent value returned [???]
                   slip_coef,    & ! slip coefficent for temperate beds [???]
                   S_period,     & ! duration of the active phase    [a]
                   Q_period,     & ! duration of the quiescent phase [a]
                   C_period,     & ! surge cyle duration             [a]
                   Time

  ! pointers to model variables
  type(Variable_t), POINTER  :: TimeVar

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Read parameters values from the constants section.
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  S_period = GetConstReal(Model % Constants, trim('surge_period'), GotIt)
  if (.not. GotIt) then
    call fatal(Caller, ' ---> Could not find \"surge_period\" in Constants Section')
  end if

  Q_period = GetConstReal(Model % Constants, trim('quiescent_period'), GotIt)
  if (.not. GotIt) then
    call fatal(Caller, ' ---> Could not find \"quiescent_period\" in Constants Section')
  end if

  slip_coef = GetConstReal(Model % Constants, trim('beta'), GotIt)
  if (.not. GotIt) then
    call fatal(Caller, ' ---> Could not find \"beta\" in Constants Section')
  end if

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! Do the simple calculation
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! Calculte the full sugre cycle lenght [yr]
  C_period = S_period + Q_period

  ! Water content threshold for ice to be "temperate", hardcoded for now
  omega_thresh = 1.0e-3

  ! Get current timestep, which is really time at end of timestep
  TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
  ! Get current time [a]
  Time    =   TimeVar % Values(1)

  ! Check if within surge period and that the bed is temperate
  if (( MOD(Time, C_period) .lt. S_period ) .and. (omega .ge. omega_thresh)) then
    beta = slip_coef
  else
    beta = 1.0_dp
  end if

end function periodic_surge
