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
