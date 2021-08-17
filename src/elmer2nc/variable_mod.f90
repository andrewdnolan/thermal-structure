module  parse_variable

  use utils
  implicit none

  type :: variable_t
    integer :: nfield, dofs, nperm
    real(kind=dp), dimension(:,:), allocatable:: data
    character(len=:), allocatable :: name, solver
  contains
    procedure :: init_data
  end type variable_t

contains
  subroutine init_data(self, nt)
    implicit none
    integer, intent(in) :: nt
    class(variable_t), intent(inout) :: self

    allocate(self%data(self%nfield, nt))
  end subroutine init_data
end module parse_variable
