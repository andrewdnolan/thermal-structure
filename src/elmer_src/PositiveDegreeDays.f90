!-------------------------------------------------------------------------------
! Computation of the positive degree days (PDD) with statistical temperature
! fluctuations; based on semi-analytical solution by Calov and Greve (2005).
!<------------------------------------------------------------------------------
module PositiveDegreeDays
  use defutils
  use special_funcs
  implicit none

  interface Cavlo_Greve_PDDs
    module procedure CV_PDDs_vstd, CV_PDDs_sstd
  end interface Cavlo_Greve_PDDs

contains

  function CV_PDDs_vstd(T, sigma)

    real(dp), dimension(:),       intent(in) :: T
    real(dp), dimension(size(T)), intent(in) :: sigma
    real(dp), dimension(size(T))             :: CV_PDDs_vstd, T_norm

    where (sigma .eq. 0.0)
    ! Handle division by zero without needing use compiler dependent ieee infinity
      T_norm       = sign(1.0_dp, T) * huge(T)
      CV_PDDs_vstd = 0.0_dp &
                   + T/2.0_dp*erfc(-T_norm)
    elsewhere
      T_norm       = T / (sqrt(2.0_dp)*sigma)
      CV_PDDs_vstd = sigma / sqrt(2.0_dp*pi) * exp(-T_norm**2.0_dp) &
                   + T/2.0_dp*erfc(-T_norm)
    end where
  end function CV_PDDs_vstd


  function CV_PDDs_sstd(T, sigma)

    real(dp), dimension(:), intent(in) :: T
    real(dp), intent(in)               :: sigma
    real(dp), dimension(size(T))       :: CV_PDDs_sstd, T_norm

    real(dp) :: inv_sqrt2pi, inv_s_stat, inv_sqrt2

    if (sigma .eq. 0.0) then
    ! Handle division by zero without needing use compiler dependent ieee infinity
      T_norm       = sign(1.0_dp, T) * huge(T)
      CV_PDDs_sstd = 0.0_dp &
                   + T/2.0_dp*erfc(-T_norm)
    else
      T_norm       = T / (sqrt(2.0_dp)*sigma)
      CV_PDDs_sstd = sigma / sqrt(2.0_dp*pi) * exp(-T_norm**2.0_dp) &
                   + T/2.0_dp*erfc(-T_norm)
    endif

  end function CV_PDDs_sstd
end module PositiveDegreeDays
