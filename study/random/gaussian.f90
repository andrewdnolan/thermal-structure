program test_normal
  !https://rosettacode.org/wiki/Statistics/Normal_distribution#Fortran
  use iso_fortran_env

  integer, parameter :: dp  = selected_real_kind(15)
  integer, parameter :: i64 = selected_int_kind(18)


  integer:: i, ind
  real(dp) :: mu = 0.0, sd = 1.0
  integer, allocatable :: seed(:)

  integer(i64), parameter :: samples = 1000000_i64
  real(dp) :: mean, stddev
  real(dp) :: sumn = 0, sumnsq = 0
  integer(i64) :: n = 0
  integer(i64) :: bin(-50:50) = 0


  ! set seed so same "random" numbers are generated everytime
  ! call set_seed(123456789)

  n = 0
  do while(n <= samples)
    rand = norm_rand(mu, sd)
    ! print *, rand

    ind = floor(5.0*rand)
    bin(ind) = bin(ind) + 1_i64
    sumn = sumn + rand
    sumnsq = sumnsq + rand*rand

    n = n + 1_i64
  end do

  mean = sumn / n
  stddev = sqrt(sumnsq/n - mean*mean)

  write(*, "(a, i0)") "sample size = ", samples
  write(*, "(a, f17.15)") "Mean :   ", mean
  write(*, "(a, f17.15)") "Stddev : ", stddev


  do i = -15, 15
    write(*, "(f4.1, a, a)") real(i)/5.0, ": ", repeat("=", int(bin(i)*500/samples))
  end do

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

end program
