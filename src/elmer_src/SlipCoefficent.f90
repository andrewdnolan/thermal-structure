module periodic_timesteps

  use DefUtils

  implicit none

  contains

    subroutine read_params(S_P, Q_P, C_P, T_f, SD_dt, QD_dt, ST_dt, QT_dt)

      ! Let this be a dummy function which "reads" our hard coded params for now

      implicit none

      real(dp), intent(out) :: S_P,   & ! (s)urge  (p)eriod
                               Q_P,   & ! (q)uies. (p)eriod
                               C_P,   & ! (c)ycle  (p)eriod
                               T_f,   & ! (t)ime (f)inal
                               SD_dt, & ! (s)urge (d)ynamic timestep
                               QD_dt, & ! (q)uise (d)ynamic timestep
                               ST_dt, & ! (s)urge (t)hermal timestep
                               QT_dt    ! (q)uise (t)heraml timestep

      S_P   = 5.00_dp
      Q_P   = 45.0_dp
      C_P   = S_P + Q_P
      T_f   = 500.0_dp
      SD_dt = 0.05_dp
      QD_dt = 1.00_dp
      ST_dt = 0.05_dp
      QT_dt = 0.10_dp
    end subroutine  read_params

    subroutine fill_timestep_sizes(timestepsizes, surging_dt, quiescesnce_dt)

      implicit none

      real(dp), intent(in)    :: surging_dt, quiescesnce_dt
      real(dp), intent(inout) :: timestepsizes(:)

      integer :: i

      do i = 1, size(timestepsizes), 2
        timestepsizes(i)   = surging_dt
        timestepsizes(i+1) = quiescesnce_dt
      end do
    end subroutine fill_timestep_sizes

    subroutine fill_timestep_intervals(timestep_intervals, &
                                       surging_dt, quiescesnce_dt, S_P, Q_P)

      implicit none

      real(dp), intent(in)    :: surging_dt, quiescesnce_dt, S_P, Q_P
      integer,  intent(inout) :: timestep_intervals(:)

      integer :: i

      do i = 1, size(timestep_intervals), 2
        timestep_intervals(i)   = floor(S_P / surging_dt)
        timestep_intervals(i+1) = floor(Q_P / quiescesnce_dt)
      end do
    end subroutine fill_timestep_intervals

    subroutine fill_dynamic_exec_intervals(dynamic_exec_intervals, &
                                           SD_dt, QD_dt, ST_dt, QT_dt)

      implicit none

      real(dp), intent(in)    :: SD_dt, QD_dt, ST_dt, QT_dt
      integer,  intent(inout) :: dynamic_exec_intervals(:)


      integer :: i

      do i = 1, size(dynamic_exec_intervals), 2
        dynamic_exec_intervals(i)   = floor(SD_dt / ST_dt)
        dynamic_exec_intervals(i+1) = floor(QD_dt / QT_dt)
      end do
    end subroutine fill_dynamic_exec_intervals

end module periodic_timesteps



function get_timestep_sizes(Model, Node, var) result(timestep_sizes)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! returns an array timesteps which is passed to the "Timestep Sizes" keyword
  ! in the Simulation section of the .sif file.
  !
  ! the timesteps are the smallest timesteps used by the model at any given
  ! point. In practice that usually means the thermal timesteps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use DefUtils
  use SolverUtils
  use ElementUtils
  use periodic_timesteps
  implicit none

  ! Solver name
  character(*), parameter :: caller = 'get_timestep_sizes'
  type(Model_t) :: Model
  integer       :: Node            ! current node number

  ! parameters which are read from the sif header
  real(dp) :: S_P,   & ! (s)urge  (p)eriod
              Q_P,   & ! (q)uies. (p)eriod
              C_P,   & ! (c)ycle  (p)eriod
              T_f,   & ! (t)ime   (f)inal
              SD_dt, & ! (s)urge  (d)ynamic timestep
              QD_dt, & ! (q)uise  (d)ynamic timestep
              ST_dt, & ! (s)urge  (t)hermal timestep
              QT_dt, & ! (q)uise  (t)heraml timestep
              var      ! dummy variable

  ! internal parameters for declaring timestep arrays
  integer  :: n_cycles,  & ! number of (surge) cycles
              M,         & ! lenght of the arrays
              i            ! loop counter

  ! resulting arrays
  real(dp), allocatable :: timestep_sizes(:)          !

  ! get the neccessary parameter values from the .sif file header
  call read_params(S_P, Q_P, C_P, T_f, SD_dt, QD_dt, ST_dt, QT_dt)
  ! find the number of surge cycle for the simulation
  n_cycles = floor(T_f/C_P)
  ! timestep array are twice as long as the number of surge-cycle
  M  = 2*n_cycles
  ! allocate the timestep arrays
  allocate(timestep_sizes(M))
  ! fill array with the smallest timestep for each phase of surge cycles
  call fill_timestep_sizes(timestep_sizes, min(SD_dt, ST_dt), min(QD_dt, QT_dt))
  ! return timesteps
  !
end function get_timestep_sizes

function get_timestep_intervals(Model, Node, var) result(timestep_intervals)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! returns an array timesteps which is passed to the "Timestep Sizes" keyword
  ! in the Simulation section of the .sif file.
  !
  ! the timesteps are the smallest timesteps used by the model at any given
  ! point. In practice that usually means the thermal timesteps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use DefUtils
  use SolverUtils
  use ElementUtils
  use periodic_timesteps
  implicit none

  ! Solver name
  character(*), parameter :: caller = 'get_timestep_intervals'
  type(Model_t) :: Model
  integer       :: Node            ! current node number

  ! parameters which are read from the sif header
  real(dp) :: S_P,   & ! (s)urge  (p)eriod
              Q_P,   & ! (q)uies. (p)eriod
              C_P,   & ! (c)ycle  (p)eriod
              T_f,   & ! (t)ime   (f)inal
              SD_dt, & ! (s)urge  (d)ynamic timestep
              QD_dt, & ! (q)uise  (d)ynamic timestep
              ST_dt, & ! (s)urge  (t)hermal timestep
              QT_dt, & ! (q)uise  (t)heraml timestep
              var      ! dummy variable

  ! internal parameters for declaring timestep arrays
  integer  :: n_cycles,  & ! number of (surge) cycles
              M,         & ! lenght of the arrays
              i            ! loop counter

  ! resulting arrays
  integer, allocatable :: timestep_intervals(:)

  ! get the neccessary parameter values from the .sif file header
  call read_params(S_P, Q_P, C_P, T_f, SD_dt, QD_dt, ST_dt, QT_dt)
  ! find the number of surge cycle for the simulation
  n_cycles = floor(T_f/C_P)
  ! timestep array are twice as long as the number of surge-cycle
  M  = 2*n_cycles
  ! allocate the timestep interval array
  allocate(timestep_intervals(M))
  ! fill array with number of intervals to execute corresponding timestep
  call fill_timestep_intervals(timestep_intervals, min(SD_dt, ST_dt), &
                               min(QD_dt, QT_dt), S_P, Q_P)
  ! return timesteps
  !
end function get_timestep_intervals

function get_dynamic_exec_intervals(Model, Node, var) result(dynamic_exec_intervals)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! returns an array timesteps which is passed to the "Timestep Sizes" keyword
  ! in the Simulation section of the .sif file.
  !
  ! the timesteps are the smallest timesteps used by the model at any given
  ! point. In practice that usually means the thermal timesteps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use DefUtils
  use SolverUtils
  use ElementUtils
  use periodic_timesteps
  implicit none

  ! Solver name
  character(*), parameter :: caller = 'get_dynamic_exec_intervals'
  type(Model_t) :: Model
  integer       :: Node            ! current node number

  ! parameters which are read from the sif header
  real(dp) :: S_P,   & ! (s)urge  (p)eriod
              Q_P,   & ! (q)uies. (p)eriod
              C_P,   & ! (c)ycle  (p)eriod
              T_f,   & ! (t)ime   (f)inal
              SD_dt, & ! (s)urge  (d)ynamic timestep
              QD_dt, & ! (q)uise  (d)ynamic timestep
              ST_dt, & ! (s)urge  (t)hermal timestep
              QT_dt, & ! (q)uise  (t)heraml timestep
              var    ! dummy variable

  ! internal parameters for declaring timestep arrays
  integer  :: n_cycles,  & ! number of (surge) cycles
              M,         & ! lenght of the arrays
              i            ! loop counter

  ! resulting arrays
  integer, allocatable :: dynamic_exec_intervals(:)

  ! get the neccessary parameter values from the .sif file header
  call read_params(S_P, Q_P, C_P, T_f, SD_dt, QD_dt, ST_dt, QT_dt)
  ! find the number of surge cycle for the simulation
  n_cycles = floor(T_f/C_P)
  ! timestep array are twice as long as the number of surge-cycle
  M  = 2*n_cycles
  ! allocate the timestep interval array
  allocate(dynamic_exec_intervals(M))
  ! fill array with number of intervals to execute corresponding timestep
  call fill_dynamic_exec_intervals(dynamic_exec_intervals, &
                                   SD_dt, QD_dt, ST_dt, QT_dt)
  ! return timesteps
  !
end function get_dynamic_exec_intervals

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

  ! 5 year surge, 25 year quiescesnce, for a 50 year surge period
  S_period = 5.0_dp
  Q_period = 25.0_dp
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
