!==============================================================================
function periodic_surge(Model, Node, omega) result(beta)
!==============================================================================
  ! This function returns the time dependent slip coefficent for periodic
  ! surge simulations.
  !
  !

  use DefUtils
  USE SolverUtils
  USE ElementUtils

  implicit none

  ! Solver name
  CHARACTER(*), PARAMETER :: Caller = 'periodic_surge'

  type(Model_t) :: Model
  integer       :: Node            ! current node number
  real(kind=dp) :: InputArray(2),& ! arguments passed to the function
                   omega,        & ! water content fraction [--]
                   omega_thresh, & ! w.c. threshold to be considered temperate [--]
                   beta,         & ! slip coefficent [???]
                   S_period,     & ! duration of the active phase    [a]
                   Q_period,     & ! duration of the quiescent phase [a]
                   C_period,     & ! surge cyle duration             [a]
                   Time
  ! pointers to model variables
  type(Variable_t), POINTER  :: TimeVar

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! TO DO: should read parameters values from the constants section.
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! ! copy input (should match the arguments!)
  ! omega = InputArray(1)
  ! time  = InputArray(2)

  ! 5 year surge, 45 year quiescesnce, for a 50 year surge period
  S_period = 5.0_dp
  Q_period = 45.0_dp
  C_period = S_period + Q_period

  ! Water content threshold for ice to be "temperate"
  omega_thresh = 1.0e-3

  ! Get current timestep, which is really time at end of timestep
  TimeVar  => VariableGet( Model % Mesh % Variables, "Time" )
  ! Get current time [a]
  Time    =   TimeVar % Values(1)

  ! Check if within surge period and that the bed is temperate
  if (( MOD(Time, C_period) .lt. S_period ) .and. (omega .ge. omega_thresh)) then
    beta = 1.0e-3
  else
    beta = 1.0_dp
  end if

end function periodic_surge
