!/******************************************************************************
! *
! *  Module containing a functions for friction heat based
! *       on residuals from (Navier-)Stokes solver
! *  This function should be preferably used to compute the heat production
! *  at the base, as it utilizes the natural way to couple the flow
! *  and temperature solution
! *
! *  Andrew Nolan Edits: 
! *    Add suport for the 'Friction Loads Scaling Factor' keyword in the 
! *    boundary condition section, so that ensure output is complient with 
! *    of thermodynamic soltion, irrespective of what the units are of the 
! *    flow solution. 
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Juha Ruokolainen, Martina Schäfer, Olivier Gagliardini, 
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 08 Jun 1997
! *  Current date:  28 August 2014 (Martina/Thomas)
! *
! *****************************************************************************/

FUNCTION getFrictionLoads(  Model, Node, DummyInput )RESULT(frictionLoad)
  
    USE DefUtils
  
    IMPLICIT NONE
    
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: Node
    REAL(KIND=dp) :: DummyInput, frictionLoad
    !----------------------------------------------------------------------------
  
    INTEGER :: DIM, i, other_body_id
    REAL(KIND=dp), POINTER :: FlowValues(:),FlowLoadValues(:),NormalValues(:),MaskValues(:)
    REAL(KIND=dp) :: normal(3), velo(3), normalvelocity, flowload(3), tangvelocity(3), ScaleFactor
    INTEGER, POINTER :: FlowPerm(:),FlowLoadPerm(:), NormalPerm(:), MaskPerm(:)
    LOGICAL :: FirstTime=.TRUE., GotIt,UnFoundFatal,UseMask = .FALSE., Warned=.FALSE.
    TYPE(Variable_t), POINTER :: FlowVar,FlowLoadVar, NormalVar, MaskVar
    TYPE(ValueList_t), POINTER :: Equation, BC
    CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolutionName, FlowLoadsName, MaskName
    CHARACTER(LEN=MAX_NAME_LEN), PARAMETER :: FunctionName='USF_GetFrictionHeating(getFrictionLoads)'
    TYPE(Element_t), POINTER ::  BoundaryElement, ParentElement
    
    SAVE FirstTime,  DIM, UseMask
  
    
    IF (FirstTime) THEN
      !WRITE(FunctionName,'(A)') 'USF_GetFrictionHeating(getFrictionLoads)'    
  
      DIM = CoordinateSystemDimension()
    END IF
     
    ! Get variable names from Equation section
    !-----------------------------------------
    BoundaryElement => Model % CurrentElement

    ! Get pointer to list corresponding to boundry condition 
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) THEN
      CALL FATAL(FunctionName,'No boundary condition found')
    END IF

    other_body_id = BoundaryElement % BoundaryInfo % outbody
    IF (other_body_id < 1) THEN ! only one body in calculation
      ParentElement => BoundaryElement % BoundaryInfo % Right
      IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
    ELSE ! we are dealing with a body-body boundary and assume that the normal is pointing outwards
      ParentElement => BoundaryElement % BoundaryInfo % Right
      IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
    END IF
    Equation => GetEquation(Element=ParentElement,Found=GotIt)
    IF (.NOT.ASSOCIATED(Equation) .OR. .NOT.GotIt) THEN
      IF (FirstTime) THEN
        WRITE (Message,'(A,I3)') 'No "Equation" found. Using default values for variables'
        CALL WARN(FunctionName,Message)    
        Warned = .TRUE.
      END IF
      WRITE(FlowSolutionName,'(A)') 'Flow Solution'
      WRITE(FlowLoadsName,'(A)') TRIM(FlowSolutionName)//' Loads'
      UseMask = .FALSE.
    ELSE
      FlowSolutionName = GetString( Equation , 'Flow Solution Name', GotIt )
      IF (.NOT.GotIt) THEN
        WRITE(FlowSolutionName,'(A)') 'Flow Solution'
        IF (FirstTime) THEN
          WRITE(Message,'(A,A)') 'Using default name for flow solution: ', &
               FlowSolutionName
          CALL WARN(FunctionName,Message)
          Warned = .TRUE.
        END IF
      END IF
      FlowLoadsName = GetString( Equation , 'Flow Loads Name', GotIt )  
      IF (.NOT. GotIt) THEN
        WRITE(FlowLoadsName,'(A)') TRIM(FlowSolutionName)//' Loads'
        IF (FirstTime) THEN
          WRITE(Message,'(A,A)') 'Using default name for flow solution loads: ', &
               FlowLoadsName
          CALL WARN(FunctionName,Message)
          Warned = .TRUE.
        END IF
      END IF
      MaskName = GetString( Equation , 'Friction Load Mask', UseMask )    
      IF (UseMask) THEN
        WRITE (Message, '(A,A)') '>Friction Load Mask< found and set to ', &
             TRIM(MaskName)
        CALL INFO(FunctionName, Message, Level=1)
      END IF
    END IF
    IF (Warned .AND. FirstTime) &
         CALL WARN(FunctionName,"All Warnings will be further omitted")
    FirstTime = .FALSE.
    
    ! ----------------------------------------------
    ! get parameters from boundary condition section
    ! ----------------------------------------------
    ScaleFactor = GetConstReal( BC,'Friction Load Scaling Factor',GotIt)
    IF (.NOT. GotIt ) ScaleFactor = 1.0_dp
        
    ! Get the variable velocity
    !---------------------------
    FlowVar => VariableGet( Model % Variables, TRIM(FlowSolutionName),UnFoundFatal=UnFoundFatal)
    FlowPerm    => FlowVar % Perm
    FlowValues  => FlowVar % Values
  
    
    ! Get the Stokes loads
    !---------------------------
    FlowLoadVar => VariableGet( Model % Variables, TRIM(FlowLoadsName),UnFoundFatal=UnFoundFatal)
    FlowLoadPerm    => FlowLoadVar % Perm
    FlowLoadValues  => FlowLoadVar % Values
    
  
    ! Get the variable for normal vector
    !-----------------------------------
    NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
    NormalPerm => NormalVar % Perm
    NormalValues => NormalVar % Values
  
    IF (UseMask) THEN
      MaskVar => VariableGet( Model % Variables, TRIM(MaskName),UnFoundFatal=UnFoundFatal)
      MaskPerm    => MaskVar % Perm
      MaskValues  => MaskVar % Values
    END IF
  
    IF (UseMask) THEN
      IF ( MaskValues(MaskPerm(Node)).LE.0.0_dp ) THEN
        frictionLoad = 0.0_dp
      END IF
    ELSE
      DO i=1, DIM
        normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
        velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
        flowload(i) = FlowLoadValues( (DIM+1)*(FlowLoadPerm(Node)-1) + i )
      END DO
   
      normalvelocity    = SUM( velo(1:DIM) * normal(1:DIM) )
  
      DO i=1, DIM
        tangvelocity(i) = velo(i) - normalvelocity * normal(i)
      END DO
      
      frictionLoad = &
           MAX( (-1.0_dp * (SUM(tangvelocity(1:DIM) * flowLoad(1:DIM)))), 0.0_dp) & 
           * ScaleFactor
    END IF
  
  END FUNCTION getFrictionLoads