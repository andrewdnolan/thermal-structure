program test_fitpack
  ! gfortran -c fitpack_interface.f90 -shared fpbspl.f splev.f
  ! gfortran test.f90 -o test fitpack_interface.o  fpbspl.o  splev.o
  use iso_fortran_env
  use fitpack_interface

  implicit none

  integer, parameter :: dp = real64
  integer            :: ier         ! error flag
  integer, parameter :: n = 9,   &  ! number of knots
                        k = 3,   &  ! spline order
                        e = 0,   &  ! extrapolate
                        m = 1000    ! number of points in x vec

  real*8, dimension(n) :: t,     &  ! vector of knots
                          c         ! vector of coefficients
  real*8, dimension(m) :: x,     &  ! where to evaluate spline (x vec)
                                    ! evaluated spline y=s(x)

  ! Load the knots
  call read_vector(t, n, 'data/knots.dat')
  ! Load the coefs
  call read_vector(c, n, 'data/coefs.dat')
  ! Load the points to evaluate spline at
  call read_vector(x, m, 'data/x_vec.dat')

  ! evaluate the spline
  call splev(t,n,c,k,x,y,m,e,ier)

  ! write the evaluated spline to disc
  call write_field(y, m, 'data/y_vec.dat')
contains


  !-----------------------------------------------------------------------------
  subroutine read_vector(vec, len, fn)
    implicit none
    integer, intent(in)                   :: len ! length of vector to be read
    real*8, dimension(len), intent(inout) :: vec ! vector to be read from disk
    character(len=*), intent(in)          :: fn  ! filename

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


  !-----------------------------------------------------------------------------
  subroutine write_field(vec, len, fn)

    implicit none
    integer, intent(in)                 :: len ! length of vector to be written
    real*8, dimension(len), intent(in)    :: vec ! vector to be written to disk
    character(len=*), intent(in)        :: fn  ! filename

    integer            :: funit,   &  ! iounit number
                          stat,    &  ! status of io read
                          i           ! loop counter


    open(unit=funit, file=fn, iostat=stat)
    if (stat /= 0) call abort()
    do i = 1, len
      write(funit, *) vec(i)
    end do
    close(funit)
  end subroutine write_field
  !-----------------------------------------------------------------------------

end program
