!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION quadratic(  Model, Node, z) RESULT(accum)
  ! provides you with most Elmer functionality
  USE DefUtils
  ! saves you from stupid errors
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model         ! the access point to everything about the model
  INTEGER       :: Node          ! the current Node number
  REAL(KIND=dp) :: z             ! nodal surface elevation [m a.s.l.]
  REAL(KIND=dp) :: accum         ! the result
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  LOGICAL       :: FirstTime=.TRUE., GotIt

  REAL(KIND=dp) :: Delta_mb   ! mass balance offset     [m i.e.q. yr^{-1}]

  ! Quadratic coefficents from numpy regression (hard_coded for now)
  REAL(KIND=dp), parameter :: p0 = -1.6877583854732304e-06, &
                              p1 = 0.012837495164557717,    &
                              p2 = -20.253792745949053

  ! Variables that only need to be read in once, are saved for future uses
  SAVE FirstTime, Delta_mb

  IF (FirstTime) THEN
    FirstTime=.FALSE.

    Delta_mb  = ListGetConstReal( Model % Constants, 'Mass Balance Offset', GotIt )
    IF (.NOT. GotIt) THEN
      CALL WARN('getAccumulation','Keyword >Mass Balance Offset< not found in >Constant< section')
      CALL WARN('getAccumulation','Taking default value >Mass Balance Offset< of 0.0 (m a^{-1})')
      Delta_mb = 0.0_dp
    END IF
  END IF

  accum = (p0*z**2 + p1*z + p2) + Delta_mb

  RETURN

END FUNCTION quadratic


! FUNCTION spline(  Model, Node, z) RESULT(accum)
!   ! provides you with most Elmer functionality
!   USE DefUtils
!   ! saves you from stupid errors
!   IMPLICIT NONE
!   ! the external variables
!   !----------------------------------------------------------------------------
!   TYPE(Model_t) :: Model         ! the access point to everything about the model
!   INTEGER       :: Node          ! the current Node number
!   REAL(KIND=dp) :: z             ! nodal surface elevation [m a.s.l.]
!   REAL(KIND=dp) :: accum         ! the result
!   !----------------------------------------------------------------------------
!   ! internal variables
!   !----------------------------------------------------------------------------
!   LOGICAL       :: FirstTime=.TRUE., GotIt
!
!   REAL(KIND=dp) :: Delta_mb   ! mass balance offset     [m i.e.q. yr^{-1}]
!
!   ! Quadratic coefficents from numpy regression (hard_coded for now)
!   REAL(KIND=dp), parameter :: p0 = -1.6877583854732304e-06, &
!                               p1 = 0.012837495164557717,    &
!                               p2 = -20.253792745949053
!
!   ! Variables that only need to be read in once, are saved for future uses
!   SAVE FirstTime, Delta_mb
!
!   IF (FirstTime) THEN
!     FirstTime=.FALSE.
!
!     Delta_mb  = ListGetConstReal( Model % Constants, 'Mass Balance Offset', GotIt )
!     IF (.NOT. GotIt) THEN
!       CALL WARN('getAccumulation','Keyword >Mass Balance Offset< not found in >Constant< section')
!       CALL WARN('getAccumulation','Taking default value >Mass Balance Offset< of 0.0 (m a^{-1})')
!       Delta_mb = 0.0_dp
!     END IF
!   END IF
!
!   accum = (p0*z**2 + p1*z + p2) + Delta_mb
!
!   RETURN
!
! END FUNCTION spline
