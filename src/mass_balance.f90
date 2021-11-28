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


FUNCTION cubic_spline(  Model, Node, z) RESULT(accum)
  ! provides you with most Elmer functionality
  USE DefUtils
  USE fitpack_interface
  ! saves you from stupid errors
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model            ! the access point to everything about the model
  INTEGER       :: Node             ! the current Node number
  REAL(KIND=dp) :: accum, &         ! the result
                   z                ! nodal surface elevation [m a.s.l.]
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  integer, parameter :: nn=9     ! number of knots in spline
  character(len=*), parameter ::   &
  knots_fp="../../input_data/mass_balance/cubic_spline_knots_s_1500_weighted.dat",   &
  coefs_fp="../../input_data/mass_balance/cubic_spline_coefs_s_1500_weighted.dat"

  LOGICAL :: FirstTime=.TRUE., GotIt
  REAL*8  :: Delta_mb                    ! mass balance offset [m i.e.q. yr^{-1}]
  REAL*8  :: x(1),                     & ! where to evaluate spline (x vec)
             y(1)                        ! evaluated spline y=s(x)

  REAL*8, dimension(nn) :: knots,        & ! vector of knots
                           coefs           ! vector of coefficients
  integer       ::         spl_err,      & ! spline error
                           alloc           ! allocaction error

  ! Variables that only need to be read in once, are saved for future uses
  SAVE FirstTime, Delta_mb, knots, coefs

  IF (FirstTime) THEN
    FirstTime=.FALSE.

    ! Read the mass balance offset from the constants section
    Delta_mb  = ListGetConstReal( Model % Constants, 'Mass Balance Offset', GotIt )
    IF (.NOT. GotIt) THEN
      CALL WARN('cubic_spline','Keyword >Mass Balance Offset< not found in >Constant< section')
      CALL WARN('cubic_spline','Taking default value >Mass Balance Offset< of 0.0 (m a^{-1})')
      Delta_mb = 0.0_dp
    END IF
    ! Load the knots
    call read_vector(knots, nn, trim(knots_fp))
    ! Load the coefs
    call read_vector(coefs, nn, trim(coefs_fp))

    !write(*,*) knots
  END IF

  ! dump our nodal surface elevation input a vector
  x(1) = z

  ! evaluate the spline
  call splev(knots,nn,coefs,3,x,y,1,0,spl_err)
  ! make sure the spline evaluation worked
  if (spl_err/=0) then
    CALL FATAL('cubic_spline', 'splev raised error')
  end if

  !write(*,*) x(1), y(1)

  ! dump the y vector to the returned scalar and add the mass balance offset
  accum = y(1) + Delta_mb
  RETURN

contains
  !-----------------------------------------------------------------------------
  subroutine read_vector(vec, len, fn)
    implicit none
    integer, intent(in)                 :: len ! length of vector to be read
    real*8, dimension(len), intent(inout) :: vec ! vector to be read from disk
    character(len=*), intent(in)        :: fn  ! filename

    integer            :: funit,   &  ! iounit number
                          stat,    &  ! status of io read
                          i           ! loop counter

    open(funit, file=fn, iostat=stat)
    if (stat /= 0) then
       write (*,*) 'Error, could not open ' // trim(fn)
       stop
    end if

    ! read data into inner regions
    do i = 1, len
       read(funit, *) vec(i)
    end do

    close(funit)
  end subroutine read_vector
  !-----------------------------------------------------------------------------
END FUNCTION cubic_spline
