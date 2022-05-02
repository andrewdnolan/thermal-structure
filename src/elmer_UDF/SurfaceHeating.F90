!------------------------------------------------------------------------------
!> Implementation of surface heating ($Q_m$) scheme (Eqn. 9) from
!> Wilson and Flowers (2013).
!------------------------------------------------------------------------------
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

  TYPE(Variable_t), POINTER :: Melt,    &
                               Firn,    &           ! firn thickness        [m]
                               Depth,   &
                               Height,  &
                               TimeVar, &
                               H_f,     &
                               MB,      &
                               Dens   ! Ice Depth       [m]

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
    ENDDO
    N_n = Model % NumberOfNodes ! Number of Nodes in Models
    N_s = i                     ! Number of surface nodes
    N_v = N_n/i                 ! Number of vertical nodes
  ENDIF
  ! End "first_time" loop

  ! Pointer to the model variables
  Melt   => VariableGet( Model % Variables, "Melting")
  Firn   => VariableGet( Model % Variables, "Firn")
  Height => VariableGet( Model % Variables, "Height")
  Depth  => VariableGet( Model % Variables, "Depth")
  Dens   => VariableGet( Model % Variables, "Densi")
  MB     => VariableGet( Model % Variables, "mass balance")

  ! Physical params
  rho_i  = GetConstReal(Model % Constants, "rho_i")    ![kg m-3]
  rho_s  = GetConstReal(Model % Constants, "rho_s")    ![kg m-3]
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
        ! Set the firn thickness
        !Firn % values (Firn % perm(n)) = MB % values (MB % perm(n)) /  0.03! [m]
        Firn % values (Firn % perm(n)) = 50 ! [m]


        ! If firn thickness exceeds ice-thickness, then set firn thickness to
        ! ice thickness
        IF (Firn%values(Firn%perm(n)) > Height%values(Height%perm(n))) THEN
          Firn % values (Firn % perm(n)) = Height%values(Height%perm(n))
        ENDIF

        ! Calculate the amount of melt
        Melt % values (Melt % perm(n)) = f_snow*SUM(PDD) * (1/rho_i) ! [m yr-1]
      ElSE
        ! Set the firn thickness and melt to 0 everywhere else
        Firn % values (Firn % perm(n)) = 0.0
        Melt % values (Melt % perm(n)) = 0.0
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
                                         (Firn % values (Firn  % perm(n)) / 1.1))

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
END SUBROUTINE Wilson_and_Flowers_2013
