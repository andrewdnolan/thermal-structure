! ******************************************************************************
! *
! *  Authors: Andrew Nolan
! *  Email:   anolan@sfu.ca
! *  Github:  andrewdnolan
! *
! *  Date Written:
! *   2021/11/06
! *
! * f90 interface to f77 code for evaluating a spline using fitpack (dierckx)
! * library. fitpack (dierckx) is what's used under the hood for spline interpolation
! * by scipy.interpolate. Spline should be fit within python, but tck tuple
! * retuned by splrep should be written to disk inorder to be evaluated in fortran
!
! ******************************************************************************
module fitpack_interface

  implicit none

  interface

    subroutine splev(t,n,c,k,x,y,m,ier)
      ! input parameters:
      !   t    : array,length n, which contains the position of the knots.
      !   n    : integer, giving the total number of knots of s(x).
      !   c    : array,length n, which contains the b-spline coefficients.
      !   k    : integer, giving the degree of s(x).
      !   x    : array,length m, which contains the points where s(x) must
      !          be evaluated.
      !   m    : integer, giving the number of points where s(x) must be
      !          evaluated.
      !
      ! output parameter:
      !   y    : array,length m, giving the value of s(x) at the different
      !          points.
      !   ier  : error flag
      !     ier = 0 : normal return
      !     ier =10 : invalid input data (see restrictions)
      real, dimension(n), intent(in)  :: t
      integer, intent(in)  :: n
      real, dimension(n), intent(in)  :: c
      integer, intent(in)  :: k
      real, dimension(n), intent(in)  :: x
      integer, intent(in)  :: m
      real, dimension(m), intent(out) :: y
      integer, intent(out) :: ier
    end subroutine splev

  end interface


end module fitpack_interface
