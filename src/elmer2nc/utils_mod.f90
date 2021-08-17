module utils
  use iso_fortran_env
  implicit none
  private
  public :: unique_sort, dp
  integer, parameter :: dp = real64

contains

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
