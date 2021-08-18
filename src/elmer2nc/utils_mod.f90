module utils
  use iso_fortran_env
  implicit none
  ! private
  ! public :: unique_sort, dp, ReadRecur
  integer, parameter :: dp = real64

contains


!------------------------------------------------------------------------------
!> Compare equality of start of s1 to (in most uses string literal) s2.
!> From elmer GeneralUtils.90
!------------------------------------------------------------------------------
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

  function unique_sort(in_array, len) result(n_unique)
    implicit none
    integer :: len, i=0,  n_unique, j, alloc_stat
    real(kind=dp) :: min_val, max_val
    real(kind=dp), dimension(len) :: in_array, unique
    real(kind=dp), dimension(:), allocatable :: out_array

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

    n_unique =  size(out_array)

  end function unique_sort
end module utils
