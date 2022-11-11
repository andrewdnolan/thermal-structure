module Thermodynamics

  use DefUtils
  implicit none

contains

  !=============================================================================
  function Conductivity_Van_Dunsen_1929(rho) result(cond)
  !=============================================================================

    implicit none

    real(kind=dp) :: rho, &  ! Ice Density          [kg m-3]
                     cond    ! Thermal Conductivity [W K-1 m-1]

    cond = 2.1d-2 + 4.2d-4 * rho + 2.2d-9 * rho**3

  end function Conductivity_Van_Dunsen_1929


  !=============================================================================
  function Conductivity_Strum_etal_1997(rho) result(cond)
  !=============================================================================

    implicit none

    real(kind=dp) :: rho, &  ! Ice Density          [kg m-3]
                     cond    ! Thermal Conductivity [W K-1 m-1]

    cond = 0.138_dp - 1.01d-3 * rho + 3.233d-6 * rho**2

  end function Conductivity_Strum_etal_1997


  !=============================================================================
  function Conductivity_Calonne_etal_2011(rho) result(cond)
  !=============================================================================

    implicit none

    real(kind=dp) :: rho, &  ! Ice Density          [kg m-3]
                     cond    ! Thermal Conductivity [W K-1 m-1]

    cond = 2.5d-6 * rho**2 - 1.23d-4 * rho + 0.024_dp

  end function Conductivity_Calonne_etal_2011


  !=============================================================================
  function Conductivity_Oster_Albert_2022(rho) result(cond)
  !=============================================================================

    implicit none

    real(kind=dp) :: rho, &  ! Ice Density          [kg m-3]
                     cond    ! Thermal Conductivity [W K-1 m-1]

    ! Conductivities of snow/firn using Oster and Albert 2022. [W K-1 m-1]
    cond = 0.144_dp * exp( 3.08d-3 * rho )

  end function Conductivity_Oster_Albert_2022

  !=============================================================================
  function Conductivity_Yen_1981(temp) result(cond)
  !=============================================================================

    implicit none

    real(kind=dp) :: temp, &  ! Ice temperature      [K]
                     cond     ! Thermal Conductivity [W K-1 m-1]

    ! Conductivity of pure ice, Eqn. 9.2 Cuffey and Paterson.  [W K-1 m-1]
    cond = 9.828_dp * exp( -5.7d-3 * Temp )


  end function Conductivity_Yen_1981

end module Thermodynamics

!==============================================================================
function  Diffusivity(Model, Node, Temp) result(EnthalpyDiffusivity)
!==============================================================================
  use Thermodynamics
  implicit none

  type(Model_t) :: Model
  type(Variable_t), pointer :: Density
  integer :: Node                        ! current node number

  real(kind=dp) :: Temp,    &             ! [K]
                   T_ptr,   &             ! [K]
                   rho,     &             ! [kg m^-3]
                   rho_i,   &             ! [kg m^-3]
                   Heat_Capacity, &       ! [J kg^-1 K^-1
                   EnthalpyDiffusivity, & ! [kg m^-1 a^-1]
                   HeatConductivity, &    ! [W K^-1 m^-1] thermal conductivity
                   K_rho,   &             ! [W K^-1 m^-1] desnity dependence of ""
                   K_ptr,   &             ! [W K^-1 m^-1] "" of ice at PMP
                   K_ice,   &             ! [W K^-1 m^-1] "" of ice
                   K_rho_i, &             ! [W K^-1 m^-1] "" at rho of ice
                   CapA,    &             ! [J kg-1 K-2]
                   CapB                   ! [J kg-1 K-1]

  character(len=*), parameter :: Caller="Diffusivity"
  character(len=256) :: dens_cond, temp_and_dens_cond
  logical :: GotIt


  ! Read from the constant section, which conductivity formula to use
  CapA    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A") ! [J kg-1 K-2]
  CapB    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B") ! [J kg-1 K-1]
  rho_i   = GetConstReal(Model % Constants, "rho_i")                    ! [kg m^-3]

  ! Find which published formula for thermal conductivity of snow/firn to use
  dens_cond = GetString( Model % Constants, "density dependent conductivity", GotIt )
  IF (GotIt) then
    ! in this case, whatever was passed was not a valid option, so just using default value
    if ((TRIM(dens_cond) .ne. "Van_Dunsen_1929")    .and. &
        (TRIM(dens_cond) .ne. "Strum_etal_1997")    .and. &
        (TRIM(dens_cond) .ne. "Calonne_etal_2011")  .and. &
        (TRIM(dens_cond) .ne. "Oster_Albert_2022")) then
      dens_cond = "Oster_Albert_2022"
      CALL INFO(Caller,"Entry found for >density dependent conductivity<. is not a valid option.", Level=6)
      CALL INFO(Caller,"Setting to 'Oster_Albert_2022'",     Level=6)
    endif
  ELSE
    ! in this case, no option was passed at all, so just using default value
    dens_cond = "Oster_Albert_2022" ! or "Gilbert_etal_2014"
    CALL INFO(Caller, "No entry found for >density dependent conductivity<.", Level=6)
    CALL INFO(Caller, "Setting to 'Oster_Albert_2022'",     Level=6)
  END IF

  ! Find how to treat the simultanoues temperature and denisty dependence
  ! of thermal conductivity
  temp_and_dens_cond = GetString( Model % Constants, "temperature and density dependent conductivity", GotIt )
  IF (GotIt) then
    ! in this case, whatever was passed was not a valid option, so just using default value
    if ((TRIM(temp_and_dens_cond) .ne. "Zwinger_etal_2007") .and. &
        (TRIM(temp_and_dens_cond) .ne. "Gilbert_etal_2014")) then
      temp_and_dens_cond = "Zwinger_etal_2007"
      CALL INFO(Caller, "Entry found for >temperature and density dependent conductivity<. is not a valid option.", Level=6)
      CALL INFO(Caller, "Setting to 'Zwinger_etal_2007'",     Level=6)
    endif
  ELSE
    ! in this case, no option was passed at all, so just using default value
    temp_and_dens_cond = "Zwinger_etal_2007" ! or "Gilbert_etal_2014"
    CALL INFO(Caller, "No entry found for >temperature and density dependent conductivity<.", Level=6)
    CALL INFO(Caller, "Setting to 'Zwinger_etal_2007'",     Level=6)
  END IF


  Density => VariableGet(Model % Variables, 'Densi')                    ! [kg m^-3]
  CapA    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A") ! [J kg-1 K-2]
  CapB    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B") ! [J kg-1 K-1]
  rho_i   = GetConstReal(Model % Constants, "rho_i")                    ! [kg m^-3]
  T_ptr   = ListGetConstReal( Model % Constants, 'T_triple', GotIt )
  IF (.NOT. GotIt) THEN
     CALL FATAL(Caller, 'T_triple not found in constant section')
  END IF

  Temp = Temp + 273.15_dp                    ! [K] <-- [C]
  rho  = Density%values(Density%perm(Node))  ! [kg m^-3]

  ! Heat capacity of ice [J kg-1 K-2]
  Heat_Capacity = CapA * Temp + CapB

  ! use the specified method to calculate the density dependence of
  ! the thermal conductivity
  SELECT CASE(TRIM(dens_cond))
    CASE("Oster_Albert_2022")
      K_rho   = Conductivity_Oster_Albert_2022(rho)   ! [W K-1 m-1]
      K_rho_i = Conductivity_Oster_Albert_2022(rho_i) ! [W K-1 m-1]
    CASE("Van_Dunsen_1929")
      K_rho   = Conductivity_Van_Dunsen_1929(rho )    ! [W K-1 m-1]
      K_rho_i = Conductivity_Van_Dunsen_1929(rho_i)   ! [W K-1 m-1]
    CASE("Calonne_etal_2011")
      K_rho   = Conductivity_Calonne_etal_2011(rho)   ! [W K-1 m-1]
      K_rho_i = Conductivity_Calonne_etal_2011(rho_i) ! [W K-1 m-1]
    CASE("Strum_etal_1997")
      K_rho   = Conductivity_Strum_etal_1997(rho)     ! [W K-1 m-1]
      K_rho_i = Conductivity_Strum_etal_1997(rho_i)   ! [W K-1 m-1]
  END SELECT

  ! using Yen 1981 formula, calculate temperature dependence of
  ! thermal conductivity
  K_ice = Conductivity_Yen_1981(temp)                   ! [W K-1 m-1]
  K_ptr = Conductivity_Yen_1981(T_ptr)                  ! [W K-1 m-1]

  ! use the specified method to account for both temp and denisty dependence
  ! of the thermal conductivity
  SELECT CASE(TRIM(temp_and_dens_cond))
    CASE("Zwinger_etal_2007")
      HeatConductivity = (K_rho / K_rho_i) * K_ice      ! [W K^-1 m^-1]
    CASE("Gilbert_etal_2014")
      HeatConductivity = (K_ice / K_ptr)   * K_rho      ! [W K^-1 m^-1]
  END SELECT

  ! Convert to units consistent with the rest of the solvers
  ! [J a^-1 kg^-1 K^-1] <-- [W K^-1 m^-1] == [J s^-1 K^-1 m^-1]
  HeatConductivity = HeatConductivity * 3600.0_dp*24.0_dp*365.25_dp

  ! Calculte diffusivity from  conductivity [kg m^-1 a^-1] <-- [J a^-1 kg^-1 K^-1]
  EnthalpyDiffusivity = HeatConductivity / Heat_Capacity
end function Diffusivity


!==============================================================================
function Lliboutry_and_Duval_Enhancment(Model, Node, omega) result(E)
!==============================================================================
  ! This function enhancment factor for the flow law when water content (omega)
  ! is present.
  !
  ! Parameterization following Greve and Blatter (2016) Eq. 14
  ! and  Wilson et al. (2013) Eqn. 5, which come from
  ! Lliboutry and Duval (1985):
  ! A = A(H)*Upsilon*(1.18125*omega*100)

  use DefUtils
  implicit none

  type(Model_t) :: Model
  integer       :: Node               ! current node number
  real(kind=dp) :: omega,  &          ! Water content fraction [--]
                   E,      &          ! Enhancment factor []
                   Upsilon = 1.0_dp   ! Switches enhancment on and off, hard coded for now

  if (omega .gt. 0.03) then
    ! We only consider enhancment up to a fractional water content of 0.03
    ! anything above that threshold is "floored" to 0.03
    omega = 0.03_dp
  elseif (omega .lt. 0.0) then
    ! If some numerical reason the fractional water content is less than 0.0
    ! set to 0.0 so there is no enhancment
    omega = 0.00_dp
  endif


  if (omega .gt. 0.00) then
    ! Calculte enhancment factor [-] for a given fractional water content
    E = Upsilon * (1.0_dp + 1.18125_dp * omega * 100.00_dp)
  else
    E = 1.0
  end if

end function Lliboutry_and_Duval_Enhancment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Limit_Enthalpy:
! Function to set the Enthalpy maxium, based on maximum water content. Following
! Wilson and Flowers (2013) we allow for a higher max. water content in the near
! surface aquifer to account for the difference in porosity between firn/snow and
! ice. As an example:
! Body Force 1
!    Enthalpy_h Upper Limit = Variable Densi, Phase Change Enthalpy
!                             Real Procedure "./bin/Thermodynamic" "Limit_Enthalpy"
! End
!
! The "Apply Limiter" keyword  must be set to true in the Enthalpy solver for the
! soft limiter to be applied.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Limit_Enthalpy_rho(Model, Node, InputArray) result(H_max)

  use DefUtils
  implicit none

  integer       :: Node
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: InputArray(2), &
                   ! depth,         &  ! depth below free sufrace     [m]
                   rho,           &  ! nodal density                [kg m-3]
                   rho_i,         &  ! ice density                  [kg m-3]
                   rho_f,         &  ! pore close off density       [kg m-3]
                   H_f,           &  ! enthalpy of fusion           [J kg-1]
                   H_max,         &  ! Limit enthalpy (returned)    [J kg-1]
                   h_aq,          &  ! firn aquifer thickness       [m]
                   w_max_en,      &  ! max englacial water content  [-]
                   w_max_aq,      &  ! max water content in frin aq [-]
                   L_heat            ! Latnet heat of fusion        [J kg-1]

  ! Read constants from sif file
  !-----------------------------
  L_heat   = GetParam(Model, "L_heat")                   ![J kg-1]
  rho_i    = GetParam(Model, "rho_i")                    ![kg m-3]
  rho_f    = GetParam(Model, "rho_f")                    ![kg m-3]
  h_aq     = GetParam(Model, "h_aq")                     ![m]
  w_max_en = GetParam(Model, "w_max_en")                 ![-]
  w_max_aq = GetParam(Model, "w_max_aq")                 ![-]
  ! Unpack input array
  !-------------------
  rho   = InputArray(1) ! Density            [kg m-3]
  H_f   = InputArray(2) ! Enthalpy of fusion [J kg-1]

  ! If density less than poreclose off, allow more water content
  if (rho .lt. rho_f) then
    ! Eqn. (10) from Aschwanden et al. 2012
    H_max = H_f + w_max_aq*L_heat
  else
    H_max = H_f + w_max_en*L_heat
  endif

  contains

  ! subfunction for reading constants and error checking
  function GetParam(Model, constant_name) result(constant)
    USE DefUtils
    implicit none

    real :: constant
    logical :: GotIt
    TYPE(Model_t) :: Model
    character(len=*) :: constant_name

    constant = GetConstReal(Model % Constants, trim(constant_name), GotIt)

    if (.not. GotIt) then
      call fatal('Limit_Enthalpy ---> GetParam', &
                 'Could not find '//trim(constant_name) )
    end if

  end function GetParam
end function Limit_Enthalpy_rho

function Limit_Enthalpy_surf(Model, Node, InputArray) result(H_max)

  use DefUtils
  implicit none

  integer       :: Node, i, N_n, N_s, N_v
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: InputArray(2), &
                   depth,         &  ! depth below free sufrace     [m]
                   H_f,           &  ! enthalpy of fusion           [J kg-1]
                   H_max,         &  ! Limit enthalpy (returned)    [J kg-1]
                   w_max_en,      &  ! max englacial water content  [-]
                   w_max_aq,      &  ! max water content in frin aq [-]
                   L_heat            ! Latnet heat of fusion        [J kg-1]
  logical :: first_time=.true.

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

  ! Read constants from sif file
  !-----------------------------
  L_heat   = GetParam(Model, "L_heat")                   ![J kg-1]
  w_max_en = GetParam(Model, "w_max_en")                 ![-]
  w_max_aq = GetParam(Model, "w_max_aq")                 ![-]

  ! Unpack input array
  !-------------------
  depth = InputArray(1) ! depth from surface [m]
  H_f   = InputArray(2) ! Enthalpy of fusion [J kg-1]

  ! If density less than poreclose off, allow more water content
  if (depth .le. 8.0) then
    ! Eqn. (10) from Aschwanden et al. 2012
    H_max = H_f + w_max_aq*L_heat
  else
    H_max = H_f + w_max_en*L_heat
  endif

  contains

  ! subfunction for reading constants and error checking
  function GetParam(Model, constant_name) result(constant)
    USE DefUtils
    implicit none

    real :: constant
    logical :: GotIt
    TYPE(Model_t) :: Model
    character(len=*) :: constant_name

    constant = GetConstReal(Model % Constants, trim(constant_name), GotIt)

    if (.not. GotIt) then
      call fatal('Limit_Enthalpy ---> GetParam', &
                 'Could not find '//trim(constant_name) )
    end if

  end function GetParam
end function Limit_Enthalpy_surf
