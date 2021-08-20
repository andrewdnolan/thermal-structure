module utils
  use netcdf
  use iso_fortran_env
  implicit none
  ! private
  ! public :: unique_sort, dp, ReadRecur
  integer, parameter :: dp = real64

contains

! read a (logical) kine from FORTRAN unit device. Inspied by "ReadAndTrim"
! function from elmer source code (https://github.com/ElmerCSC/elmerfem/blob/6dfb482454dba8245cf35d0e1591927156b6e1ec/fem/src/GeneralUtils.F90#L842)
!-------------------------------------------------------------------------------
  recursive function ReadRecur( Unit, str ) result(l)
!-------------------------------------------------------------------------------
    integer :: Unit                       !< Fortran unit number to read from
    character(len=:), allocatable :: str  !< The string read from the file
    LOGICAL :: l                          !< Success of the read operation
!-------------------------------------------------------------------------------

    integer :: outlen ! length of output character array
    integer, parameter :: maxlen = 16384
    character(len=maxlen) :: readstr = ' '

    l = .TRUE.

    read( Unit,'(A)',end=10,err=10 ) readstr

    outlen = len(trim(readstr))
    if(.not.allocated(str)) allocate(character(512)::str)
    str(1:outlen) = trim(readstr)
    str(outlen+1:)  = ''
    return
    ! Return false (ending do while loop) when we read the end of file,
    ! or if there is an error
10  CONTINUE
    l = .FALSE.
  end function
!-------------------------------------------------------------------------------

! Helper Function for dealing with netcdf calls
!-------------------------------------------------------------------------------
  subroutine nc_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine nc_check
!-------------------------------------------------------------------------------

  function unique_sort(in_array) result(out_array)
    implicit none
    integer :: len, i=0,  n_unique, j, alloc_stat
    real(kind=dp) :: min_val, max_val
    real(kind=dp) :: in_array(:)
    real(kind=dp), dimension(:), allocatable :: out_array, unique

    if (allocated(unique)) deallocate(unique)

    allocate(unique(size(in_array)))

    min_val = minval(in_array) - 1.0
    max_val = maxval(in_array)

    i = 0
    do while (min_val<max_val)
      i = i+1
      min_val = minval(in_array, mask=in_array>min_val)
      unique(i) = min_val
    end do

    if (allocated(out_array)) deallocate(out_array)
    allocate(out_array(i), source=unique(1:i))
    if (alloc_stat /= 0) call abort()

  end function unique_sort




!-------------------------------------------------------------------------------
! ELMER FUNCTIONS DIRECTLY FROM ELMER
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Compare equality of start of s1 to (in most uses string literal) s2.
!> From elmer GeneralUtils.90
!------------------------------------------------------------------------------
  PURE FUNCTION SEQL(s1,s2) RESULT(L)
!------------------------------------------------------------------------------
    LOGICAL :: L
    CHARACTER(LEN=*), INTENT(IN) :: s1,s2
!------------------------------------------------------------------------------
    INTEGER :: n
!------------------------------------------------------------------------------
    L = .FALSE.
    n = LEN(s2)
    IF(LEN(s1) < n) RETURN
    IF (s1(1:n)==s2) L=.TRUE.
!------------------------------------------------------------------------------
  END FUNCTION SEQL
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Sort an real array, and change the order of an index array accordingly.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortD( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: b(:)
     REAL(KIND=dp), INTENT(inout) :: a(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,rb
     REAL(KIND=dp) :: ra
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN

      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortD
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine print( Caller, String, Verbose)
!------------------------------------------------------------------------------
    implicit none
    character(len=*)  :: Caller  !< The function where print is called from
    character(len=*)  :: String  !< The message to print
    logical, optional :: Verbose !< Wether the message should be printed or not
!------------------------------------------------------------------------------


  if (.not. present(Verbose)) Verbose=.FALSE.

  if ( Verbose ) then
    write(*,"(A)") trim(Caller) // ": " // trim(String)
  else
    return
  end if

  end subroutine print
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine fatal( Caller, String)
!------------------------------------------------------------------------------
    implicit none
    character(len=*)  :: Caller  !< The function where print is called from
    character(len=*)  :: String  !< The message to print
!------------------------------------------------------------------------------


  write(*,"(A)") trim(Caller) // ": " // trim(String)
  stop

  end subroutine fatal
!------------------------------------------------------------------------------
end module utils
