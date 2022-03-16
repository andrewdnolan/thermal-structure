! ******************************************************************************
! *
! *  Authors: Andrew Nolan
! *  Email:   anolan@sfu.ca
! *  Github:  andrewdnolan
! *
! *  Date Written:
! *   2021/09/21
! *
! *  Note:
! *   This subroutine was adapted from material provided at the Elmer/Ice
! *   Beginners course at UiO in 2016. The original material can be found on the
! *   Documentation page of the Elmer/Ice Wiki under "Course Material."
! * https://github.com/ElmerCSC/elmerfem/blob/d2cfc4b8ffc6061c8088b85283d5b925ea325b6f/elmerice/Solvers/SurfaceBoundaryEnthalpy.F90#L964
! ******************************************************************************

subroutine Percolation(Model, Solver, dt, Transient)
  ! *****************************************************************************
  !> Percolation.f90, subroutine Percolation
  !>
  !> Solve for latent heat released by meltwater refreezing in the firn aquifer
  !> following Gilbert et al. 2014/2020.
  !>
  !> Surface melt ("Surf_melt") must be specified as boundary condition along the
  !> free surface for this subroutine to work correctly.
  ! *****************************************************************************
  use DefUtils
  implicit none

  type(Model_t)  :: Model
  type(Solver_t) :: Solver
  type(Element_t),   pointer :: Element
  type(Variable_t),  pointer :: Depth, Enthalpy, H_f
  type(ValueList_t), pointer :: Material, BC
  integer,           pointer :: NodeIndexes(:)

  integer ::       i, j,         & ! Loop indexes
                   k,            & ! Vertical node indexes
                   N,            & ! Total number of nodes within mesh
                   NN,           & ! number of nodes for an individual element
                   Nsurf           ! number of surface nodes

  real(kind=dp) :: dt,           & ! timestep size         [a]
                   dz,           & ! layer thickness       [m]
                   Sr,           & ! residual saturation   [-]
                   wres,         & !
                   rho_w,        & ! water desnity         [kg m^-3]
                   rho_i,        & ! ice   density         [kg m^-3]
                   L_heat,       & ! Latent heat of fusion [J kg^-1]
                   melt_surf,    & ! Nodal surface melting [kg m^-2]
                   H_limit,      & ! Enthalpy of fusion    [J kg^-1]
                   filler

  real(kind=dp), allocatable ::  &
                   Dens(:),      & !
                   Density(:),   & !
                   Melt(:),      & !
                   Melting(:)      !

  logical :: Transient, GotIt, first_time=.true.

  character(len=*), parameter :: Caller = "Percolation"

  save first_time, N, Nsurf

  ! Pointers to variable fields
  ! ----------------------------------------------------------------------------
  Depth    => VariableGet( Model % Variables, "Depth")
  Enthalpy => VariableGet( Model % Variables, "enthalpy_h")
  H_f      => VariableGet( Model % Variables, "Phase Change Enthalpy")

  ! First time loop:  find number of nodes and read constants from .sif
  ! ----------------------------------------------------------------------------
  if (first_time) then
    first_time=.false.

    N     = Model % NumberOfNodes
    Nsurf = 0

    ! Find the number of surface nodes
    do i = 1, N
      if ( Depth % Values ( Depth % perm(i) ) == 0.0) Nsurf = Nsurf + 1
    end do

    ! Read constants from .sif file
    Sr     = GetConstReal( Model % Constants, "Sr")
    rho_w  = GetConstReal( Model % Constants, "rho_w"  )
    rho_i  = GetConstReal( Model % Constants, "rho_ice")
    L_heat = GetConstReal( Model % Constants, "L_heat" )

  ! end first_time loop
  end if

  ! Density in SI units [kg m^-3] is defined as a material parameter, assigned
  ! to the variable "Enthalpy Density". Therefore, to access the values we need
  ! read the material values, instead of the normal (and more compact) pointer
  ! "VariableGet(...)" method
  ! ----------------------------------------------------------------------------
  allocate(Density(N))
  allocate(Dens(Model % MaxElementNodes))

  do i = 1, Solver % NumberOfActiveElements
    Element     => GetActiveElement(i)
    Material    => GetMaterial(Element)
    NodeIndexes => Element % NodeIndexes

    NN  =  GetElementNOFNodes(Element)
    Dens(1:NN) = ListGetReal( Material, "Enthalpy Density", NN, NodeIndexes, GotIt )
    if (.not.GotIt) then
      call WARN(Caller, "Problem reading ""Enthalpy Denisty"" from Material Section")
      return
    endif

    Density(NodeIndexes) = Dens(1:NN)
  end do

  ! Melting is solved for as part of the mass balance model, but specified along
  ! the top surface as a boundary value. We follow a similar approach as we did
  ! for "Enthalpy Density", but instead of looping over the material we will loop
  ! over the boundary elements
  ! ----------------------------------------------------------------------------
  allocate(Melt(N))
  allocate(Melting(Model % MaxElementNodes))
  do i = 1, Model % NumberOfBoundaryElements
    Element     => GetBoundaryElement(i)
    BC          => GetBC(Element)
    NodeIndexes => Element % NodeIndexes

    NN = GetElementNOFNodes(Element)
    Melting(1:NN) = ListGetReal( BC, "Surf_melt", NN, NodeIndexes, GotIt)
    if (.not.GotIt) then
      call WARN(Caller, "Problem reading ""Surf_melt"" from Boundary Condition Section")
      return
    endif

    Melt(NodeIndexes) = Melting(1:NN)
  end do

  ! ----------------------------------------------------------------------------
  ! Actually Do Calculation:
  ! ----------------------------------------------------------------------------

  ! Loop over surface nodes
  do i = 1, Nsurf
    ! vertical node loop number
    j = 0
    ! Surface Melting [kg m-2]
    melt_surf = Melt(N+1-i) * rho_w

    ! Continue loop while there is still surface melt to refreeze, once all the
    ! melt has refrozen the loop breaks
    do while (melt_surf .ne. 0.0)
      ! get vertically aligned node number
      k  = N+1-i-j*Nsurf
      ! Get current vertical layer thickness
      dz = Depth % Values(Depth % perm(k-Nsurf))-Depth % Values(Depth % perm(k))

      ! pore close off at ~800 kg m^-3 ??
      if (Density(k) < 800.0) then

        ! Maximum residual water content in kg m^-3
        wres = Sr*(1.0 - Density(k)/rho_i)*rho_w

        write(*,*) wres
        ! Enthalpy of fusion for maximum residual water content
        H_limit = Enthalpy % values(Enthalpy % perm(k)) + wres/Density(k)*L_heat

        ! Latent heat released from meltwater refreezing
        Enthalpy % values( Enthalpy % perm(k)) = Enthalpy % values( Enthalpy % perm(k)) &
                                                 + L_heat*melt_surf/dz/Density(k)

        ! Enthalpy is limited by the enthalpy of fusion (calculated above).
        ! If we've reached the enthalpy of fusion then reduce the amount of
        ! meltwater available for refreezing
        if (Enthalpy % values( Enthalpy % perm(k)) > H_limit) then
          melt_surf = (Enthalpy % values( Enthalpy % perm(k)) - H_limit) * &
                       dz*Density(k)/L_heat
          Enthalpy % values( Enthalpy % perm(k)) = H_limit
        else
          melt_surf = 0.0
        endif
      else
        melt_surf = 0.0
      end if

      ! Go onto next vertical node
      j = j + 1

      ! If we've iterated through all the vertical nodes, then break the while loop
      if ( j == N/Nsurf-1 ) then
        melt_surf = 0.0
      end if

    end do
  end do
end subroutine Percolation
