module  parse_variable

  use utils
  implicit none

!-------------------------------------------------------------------------------
  type :: variable_t
!-------------------------------------------------------------------------------
    integer :: &
    nfield,    &      ! size of variable field
    nperm,     &      ! size of associated permutation
    dofs,      &      ! Variable degrees of freedom, should always be 1
    nc_varid          ! variable id, will be assigned by the netcdf module

    real(kind=dp), allocatable :: values(:) ! array to store variable data
    integer,       allocatable :: perm(:)   ! variable permutation table

    character(len=:), allocatable :: name,  & ! name of the variable
                                     solver   ! name of solver variable is from
  contains
    procedure :: allocate
    procedure :: deallocate
  end type variable_t

contains

!-------------------------------------------------------------------------------
! Allocate the values and perm array for this instance of the variable
!-------------------------------------------------------------------------------
  subroutine allocate(self)
!-------------------------------------------------------------------------------
    implicit none
    integer :: err
    class(variable_t), intent(inout) :: self

    allocate(self%values(self%nfield), stat=err)
    if (err /=0 ) then
      call fatal('variable_t', 'problem allocating values array for'//trim(self%name))
    end if

    allocate(self%perm(self%nfield), stat=err)
    if (err /=0 ) then
      call fatal('variable_t', 'problem allocating perm array for'//trim(self%name))
    end if

  end subroutine allocate

!-------------------------------------------------------------------------------
! Allocate the values and perm array for this instance of the variable
!-------------------------------------------------------------------------------
  subroutine deallocate(self)
!-------------------------------------------------------------------------------
    implicit none
    integer :: err
    class(variable_t), intent(inout) :: self

    if (allocated(self%values)) deallocate(self%values, stat=err)
    if (err /=0 ) then
      call fatal('variable_t', 'problem deallocating values array for'//trim(self%name))
    end if

    if (allocated(self%perm)) deallocate(self%perm, stat=err)
    if (err /=0 ) then
      call fatal('variable_t', 'problem deallocating perm array for'//trim(self%name))
    end if

  end subroutine deallocate

end module parse_variable
