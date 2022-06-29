module SurfaceTemperature

  implicit none

contains

  subroutine SurfTemp(z, T, alpha, grad_T, ref_z, T_mean, T_peak, coefs, seed)
    !
    ! calculate surface temp with random (but consistent) variability
    !

    use DefUtils

    implicit none

    ! input params
    Integer, intent(in)       :: T_peak     ! DOY of annual temp peak [DOY]
    Integer, intent(in), &
             optional         :: seed       ! seed for random num. generator
    real(kind=dp), intent(in) :: z,       & ! nodal surface elevation [m]
                                 alpha,   & ! Anual air temp. amp.    [C]
                                 grad_T,  & ! Air temp lapse rate     [K m^-1]
                                 ref_z,   & ! Reference surface elev. [m a.s.l.]
                                 T_mean,  & ! Mean annual T @ ref_z   [C]
                                 coefs(3)   ! Coefs for T_std(doy)    [?]
    ! resulting vector
    real(kind=dp), intent(out) :: T(365)    ! T(DOY) @ nodal z          [m]

    ! internal params
    integer :: i,       &                   ! index counter for DOY
              seed__                        ! internal value of seed
    real(kind=dp) :: std                    ! standard dev. for DOY     [?]

    ! Check if seed was passed
    if (present(seed)) then
      seed__ = seed
    else
      seed__ = 123456789
    end if

    ! seed the random number generator
    call set_seed(seed)

    ! Itterate over the julian calendar days
    DO i=1,365

        ! Calculate the daily standard deviation
        call temp_std(std, i, coefs)

        ! Find surface temp for DOY(i)
        T(i) = T_mean                                       & ! mean signal
             + alpha*cos(2.0_dp*PI*float(i - T_peak)/365.0) & ! seasonal cycle
             + grad_T*(ref_z-z)                             & ! elevation dependence
             + norm_rand(0.0_dp, std)                         ! random variability

    ENDDO

  contains

    subroutine temp_std(std, d, coefs)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Time dependent standard deviation in air temperature, which facilitates
      ! greater varaince in winter than summer, but is also neccisary for
      ! recreating melt curve
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      integer,       intent(in)  :: d            ! Day of year            [DOY]
      real(kind=dp), intent(in)  :: coefs(3)     ! Coefs for T_std(doy)   [?]
      real(kind=dp), intent(out) :: std          ! standard deviation for DOY

      ! evaluate the quadaratic function for the daily std
      std = coefs(1)**2 * d + coefs(2) * d + coefs(3)
    end subroutine temp_std

    subroutine set_seed(seed)
      !https://masuday.github.io/fortran_tutorial/random.html
      implicit none
      integer :: n                      ! length of seed array
      integer, intent(in) :: seed       ! interger to seed all elements with
      integer, allocatable :: seed_a(:) ! seed vector

      ! find the number of elements in seed array
      call random_seed(size=n)

      ! allocate seed array
      allocate(seed_a(n))
      ! put seed passed to subroutine in at all elements
      seed_a = seed
      ! seed random number generator w/ our seed array
      call random_seed(put=seed_a)
      ! deallocate seed array
      deallocate(seed_a)
    end subroutine set_seed

    function norm_rand(mean, std_dev)
        ! Draw random samples from a normal (Gaussian) distribution with passed
        ! mean and standard devitaion.
        !
        !
        implicit none

        real(KIND=dp) :: norm_rand
        real(KIND=dp), intent(in) :: mean, std_dev
        real(KIND=dp) :: v1, v2, r
        real(KIND=dp), save :: spare
        logical, save :: has_spare


        if (has_spare) then                    ! We have an extra deviate handy,
          norm_rand = mean + (std_dev * spare) ! so return it,
          has_spare = .FALSE.                  ! and unset the flag
          return
        else                                   ! We do not have extra deviate handy,
            r = 1.0
            do while ( r >= 1.0 )
                ! pick two random number between 0 and 1
                CALL RANDOM_NUMBER(v1)
                CALL RANDOM_NUMBER(v2)
                ! remap random number to square extending from -1 to +1
                v1 = (v1 * 2.0) - 1.0
                v2 = (v2 * 2.0) - 1.0
                ! see if variates within unit circle (while loop)
                r = v1**2 + v2**2
            end do

            ! make Box-Muller transformation to get two normal deviates
            r = sqrt((-2.0 * log(r)) / r)

            ! scale the deviate by the specified mean and standard deviation
            norm_rand = mean + (std_dev * v1 * r)
            ! store the transformed (but unscaled) deviate for next time
            spare = v2 * r
            has_spare = .TRUE.
            return
        end if
    end function norm_rand

  end subroutine SurfTemp

end module SurfaceTemperature
