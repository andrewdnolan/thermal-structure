module special_funcs
  ! Use elmer type definitions to avoid numeric errors in B.C. code
  use DefUtils
  implicit none

  private            ! All entities are now module-private
  public erfc        ! Explicitly export public entities

  interface erfc
    MODULE PROCEDURE erfc_s, erfc_v
  end interface

  contains

  function erfc_v(x)
    ! Vecotrized complementary error function [Press et al. 1996]
    ! https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
    real(kind=dp), dimension(:), intent(in) :: x
    real(kind=dp), dimension(size(x))       :: erfc_v


    erfc_v = tau_v(x)
    where (x .lt. 0.0) erfc_v = 2 - erfc_v

    contains

    function tau_v(x)
      real(kind=dp), dimension(:), intent(in) :: x
      real(kind=dp), dimension(size(x))       :: t, tau_v, z

      z = abs(x)
      t = 1 / (1 + 0.5 * z)

      tau_v = t * exp(-x**2   &
            - 1.26551223      &
            + 1.00002368*t    &
            + 0.37409196*t**2 &
            + 0.09678418*t**3 &
            - 0.18628806*t**4 &
            + 0.27886807*t**5 &
            - 1.13520398*t**6 &
            + 1.48851587*t**7 &
            - 0.82215223*t**8 &
            + 0.17087277*t**9 )
    end function tau_v
  end function erfc_v

  function erfc_s(x)
    ! Scalar complementary error function [Press et al. 1996]
    ! https://en.wikipedia.org/wiki/Error_function#Numerical_approximations
    real(kind=dp), intent(in) :: x
    real(kind=dp)             :: erfc_s


    erfc_s = tau_s(x)
    if (x .lt. 0.0) erfc_s = 2 - erfc_s

    contains

    function tau_s(x)
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: t, tau_s, z

      z = abs(x)
      t = 1 / (1 + 0.5 * z)

      tau_s = t * exp(-x**2   &
            - 1.26551223      &
            + 1.00002368*t    &
            + 0.37409196*t**2 &
            + 0.09678418*t**3 &
            - 0.18628806*t**4 &
            + 0.27886807*t**5 &
            - 1.13520398*t**6 &
            + 1.48851587*t**7 &
            - 0.82215223*t**8 &
            + 0.17087277*t**9 )
    end function tau_s
  end function erfc_s

end module special_funcs
