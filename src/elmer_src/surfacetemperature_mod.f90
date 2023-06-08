module SurfaceTemperature

  implicit none

contains

  subroutine SurfTemp(T, z, time_i, time_ip1, alpha, grad_T, ref_z, T_mean, T_peak)
    !
    ! calculate surface temp. 
    !
    ! returns an array containing the surface air temperature for 100 linearly spaced 
    ! points in time between `time_i` and `time_ip1`

    use DefUtils

    implicit none

    ! input params
    integer, parameter        :: n=100       ! fixed length of the time vector
    real(kind=dp), intent(in) :: z,        & ! nodal surface elevation [m]
                                 time_i,   & ! time at beginging of the current timestep [a]
                                 time_ip1, & ! time at    end    of the current timestep [a]
                                 alpha,    & ! Anual air temp. amp.    [C]
                                 grad_T,   & ! Air temp lapse rate     [K m^-1]
                                 ref_z,    & ! Reference surface elev. [m a.s.l.]
                                 T_mean,   & ! Mean annual T @ ref_z   [C]
                                 T_peak      ! DOY of annual temp peak [DOY]
                                 ! coefs(3)   ! Coefs for T_std(doy)    [?]
  
    ! resulting vector
    real(kind=dp), intent(inout) :: T(n)    ! T(DOY) @ nodal z          [m]

    ! internal params
    integer :: i                ! index counter 
    real(kind=dp) :: time(n), & ! time array [a]
                     h          ! step between `time_i` and `time_ip1`

    ! mimicing `np.linspace` over the current timestep range 
    ! ref: https://math.unm.edu/~motamed/Teaching/OLD/Fall20/HPSC/fortran.html#more-on-arrays
    ! get the increment [a] betweent the two timesteps
    h = (time_ip1 - time_i)/real(n-1, dp)
    ! populate time array, in years
    ! NOTE: Need to itterate b/w `0` and `n-1` for linspace like calc to work properly. 
    time = time_i + h * (/(real(i, dp), i = 0, n-1)/)

    DO i=1,n
        ! Find surface temp for DOY(i+1), index is `i+1` (instead of just `i`) to 
        ! componsate for the zero indexing in the `DO` loop
        T(i) = T_mean                                               & ! mean signal
               + alpha*cos(2.0_dp*PI*(time(i) - T_peak/365.0_dp))   & ! seasonal cycle
               + grad_T*(ref_z-z)                                     ! elevation dependence
               ! + norm_rand(0.0_dp, std)                             ! random variability
    ENDDO

  contains
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
