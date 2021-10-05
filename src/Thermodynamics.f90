function  Diffusivity(Model, Node, Temp) result(EnthalpyDiffusivity)
  use DefUtils
  implicit none

  type(Model_t) :: Model
  type(Variable_t), pointer :: Density
  integer :: Node                        ! current node number

  real(kind=dp) :: Temp, &                ! [K]
                   T_ptr, &               ! [K] Temp. (of water) at triple point
                   rho, &                 ! [kg m^-3]
                   Heat_Capacity, &       ! [J kg^-1 K^-1
                   EnthalpyDiffusivity, & ! [kg m^-1 a^-1]
                   HeatConductivity, &    ! [W K^-1 m^-1] thermal conductivity
                   K_rho, &               ! [W K^-1 m^-1] desnity dependence of ""
                   K_ice, &               ! [W K^-1 m^-1] "" of ice
                   K_ptr, &               ! [W K^-1 m^-1] "" at triple point of water
                   CapA,  &               ! [J kg-1 K-2]
                   CapB,  &               ! [J kg-1 K-1]
                   CondA, &               ! [W m^5 K^-1 kg^-2]
                   CondB, &               ! [W m^2 K^-1 kg^-1]
                   CondC                  ! [W kg-1 m-1]

  Density => VariableGet(Model % Variables, 'Densi')                        ! [kg m^-3]
  CapA    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity A")     ! [J kg-1 K-2]
  CapB    = GetConstReal(Model % Constants, "Enthalpy Heat Capacity B")     ! [J kg-1 K-1]
  CondA   = GetConstReal(Model % Constants, "Enthalpy Heat Conductivity A") ! [W m^5 K^-1 kg^-2]
  CondB   = GetConstReal(Model % Constants, "Enthalpy Heat Conductivity B") ! [W m^2 K^-1 kg^-1]
  CondC   = GetConstReal(Model % Constants, "Enthalpy Heat Conductivity C") ! [W kg-1 m-1]
  T_ptr   = GetConstReal(Model % Constants, "T_triple")                     ! [K]

  Temp = Temp + 273.15                       ! [K] <-- [C]
  rho  = Density%values(Density%perm(Node))  ! [kg m^-3]

  ! Heat capacity of ice [J kg-1 K-2]
  Heat_Capacity = CapA * Temp + CapB

  ! Intermediate conductivity calcs. [W kg-1 m-1]]
  K_rho = CondA * rho**2 - CondB * rho + CondC
  K_ice = 9.828*exp(-5.7e-3 *(Temp))
  K_ptr = 9.828*exp(-5.7e-3 *(T_ptr))

  ! Conductivity  [J a^-1 kg^-1 K^-1] <-- [W K^-1 m^-1] == [J s^-1 K^-1 m^-1]
  HeatConductivity = (K_ice / K_ptr * K_rho) * 3600.0*24.0*365.25

  ! Diffusivity [kg m^-1 a^-1]
  EnthalpyDiffusivity = HeatConductivity / Heat_Capacity
end function Diffusivity
