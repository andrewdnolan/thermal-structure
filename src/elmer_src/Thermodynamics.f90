function  Diffusivity(Model, Node, Temp) result(EnthalpyDiffusivity)
  use DefUtils
  implicit none

  type(Model_t) :: Model
  type(Variable_t), pointer :: Density
  integer :: Node                        ! current node number

  real(kind=dp) :: Temp,    &             ! [K]
                   rho,     &             ! [kg m^-3]
                   rho_i,   &             ! [kg m^-3]
                   Heat_Capacity, &       ! [J kg^-1 K^-1
                   EnthalpyDiffusivity, & ! [kg m^-1 a^-1]
                   HeatConductivity, &    ! [W K^-1 m^-1] thermal conductivity
                   K_rho,   &             ! [W K^-1 m^-1] desnity dependence of ""
                   K_ice,   &             ! [W K^-1 m^-1] "" of ice
                   K_rho_i, &             ! [W K^-1 m^-1] "" at rho of ice
                   CapA,    &             ! [J kg-1 K-2]
                   CapB                   ! [J kg-1 K-1]

  Density => VariableGet(Model % Variables, 'Densi')                    ! [kg m^-3]
  CapA    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A") ! [J kg-1 K-2]
  CapB    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B") ! [J kg-1 K-1]
  rho_i   = GetConstReal(Model % Constants, "rho_i")                    ! [kg m^-3]

  Temp = Temp + 273.15                       ! [K] <-- [C]
  rho  = Density%values(Density%perm(Node))  ! [kg m^-3]

  ! Heat capacity of ice [J kg-1 K-2]
  Heat_Capacity = CapA * Temp + CapB

  ! Conductivities of snow/firn using Oster and Albert 2022. [W K-1 m-1]
  K_rho   = 0.144 * exp( 3.08e-3 * rho )
  K_rho_i = 0.144 * exp( 3.08e-3 * rho_i )
  ! Conductivity of pure ice, Eqn. 9.2 Cuffey and Paterson.  [W K-1 m-1]
  K_ice   = 9.828 * exp( -5.7e-3 * Temp )

  ! Conductivity  [J a^-1 kg^-1 K^-1] <-- [W K^-1 m^-1] == [J s^-1 K^-1 m^-1]
  ! accounting for both temp and denisty dependence [Eqn. 18 Zwinger et al. 2007]
  HeatConductivity = (K_rho / K_rho_i) * K_ice * 3600.0*24.0*365.25

  ! Diffusivity [kg m^-1 a^-1]
  EnthalpyDiffusivity = HeatConductivity / Heat_Capacity
end function Diffusivity

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

  integer       :: Node
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: InputArray(2), &
                   depth,         &  ! depth below free sufrace     [m]
                   H_f,           &  ! enthalpy of fusion           [J kg-1]
                   H_max,         &  ! Limit enthalpy (returned)    [J kg-1]
                   w_max_en,      &  ! max englacial water content  [-]
                   w_max_aq,      &  ! max water content in frin aq [-]
                   L_heat            ! Latnet heat of fusion        [J kg-1]

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
  if (depth .eq. 0.0) then
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
