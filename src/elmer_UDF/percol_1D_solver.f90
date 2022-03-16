SUBROUTINE percol_1D_solver(Model, Solver, dt, Transient)
Use DefUtils
IMPLICIT NONE
TYPE(Solver_t) :: Solver
TYPE(Model_t) :: Model
TYPE(ValueList_t), POINTER :: BC,Material,SolverParams
TYPE(Element_t), POINTER :: Element
INTEGER, POINTER :: NodeIndexes(:)
TYPE(variable_t), POINTER :: Enth,Hs,Depth
REAL(KIND=dp) :: dt,dz,Sr,rho_w,rho_ice,L_heat,wres,H_lim,cont
LOGICAL :: Transient,GotIt,Found
INTEGER :: N,i,j,nb,k,Nsurf
real,dimension(:),allocatable :: dens,fonte
REAL(KIND=dp), ALLOCATABLE :: Density(:),fonte_loc(:)

logical :: first_time=.true.

SAVE first_time,Sr,rho_w,rho_ice,L_heat,N,Nsurf

if (first_time) then
  first_time=.false.

N = Model % NumberOfNodes

!========Get number of surface node===================================================

Depth => VariableGet( Model % Variables, 'Depth')
Nsurf=0
DO i=1,N
IF (Depth % values (Depth % perm (i))==0.0) THEN
Nsurf=Nsurf+1
ENDIF
ENDDO


rho_w=GetConstReal(Model % Constants, "rho_w")
rho_ice=GetConstReal(Model % Constants, "rho_ice")
Sr=GetConstReal(Model % Constants, "Sr")
L_heat=GetConstReal(Model % Constants, "L_heat")

endif

!==========Density and capacity from materials=======================================


ALLOCATE(dens(N))
ALLOCATE(Density(Model % MaxElementNodes))

do i=1,Solver % NumberOfActiveElements

   Element => GetActiveElement(i)
   nb = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes
   Material => GetMaterial(Element)
   Density(1:nb) = ListGetReal( Material,'Enthalpy Density', nb, NodeIndexes)
   dens(NodeIndexes)=Density(1:nb)

enddo


!======Get enthalpy and depth================================================

Hs => VariableGet( Model % Variables, 'Phase Change Enthalpy' )
Enth => VariableGet( Model % Variables, 'enthalpy_h')
Depth => VariableGet( Model % Variables, 'Depth')

!===================================================================================
!========Get surface melting=========================================================
ALLOCATE(fonte_loc(Model % MaxElementNodes))
ALLOCATE(fonte(N))

fonte=0.0
fonte_loc=0.0

DO i=1,Model % NumberOfBoundaryElements

   Element => GetBoundaryElement(i)
   nb = GetElementNOFNodes(Element)
   NodeIndexes => Element % NodeIndexes

   BC => GetBC(Element)

   fonte_loc(1:nb) = ListGetReal( BC,'Surf_melt', nb, NodeIndexes,GotIt)
   IF (.NOT.GotIt) THEN
	!fonte(NodeIndexes)=0.0
   ELSE
        fonte(NodeIndexes)=fonte_loc(1:nb)/(1.0)!======Correction as a function of timestep=======
   END IF

ENDDO

!===========================================================================================
!=========Loop over surface nodes=========================================================
DO i=1,Nsurf
!Surface melting en kg m-2:
cont=fonte(N+1-i)*rho_w
j=0


!=========Loop over verticaly aligned nodes ===============================================
DO WHILE (cont.ne.0.0)

!Curent node number:
k = N+1-i-j*Nsurf

!layer thickness:
dz = Depth % values (Depth % perm (k-Nsurf))-Depth % values (Depth % perm (k))

IF (dens(k)<800.0) THEN

!===========Residual water content========================================================

! kg m^{-3}
wres=Sr*(1.0-dens(k)/rho_ice)*rho_w

H_lim = Hs % values (Hs % perm (k)) + wres/dens(k)*L_heat

!================percolation ============================================================


Enth % values (Enth % perm (k))= Enth % values (Enth % perm (k)) + L_heat*cont/dz/dens(k)

if (Enth % values (Enth % perm (k))> H_lim) then
	cont=(Enth % values (Enth % perm (k))-H_lim)*dz*dens(k)/L_heat
        Enth % values (Enth % perm (k))=H_lim
else
	cont=0.0
endif

ELSE

cont=0.0

ENDIF

j=j+1

if (j==Model % NumberOfNodes/Nsurf-1) then
	cont=0.0
endif



ENDDO
ENDDO


END SUBROUTINE percol_1D_solver
