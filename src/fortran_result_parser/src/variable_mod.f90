module  parse_variable

  implicit none

  type :: variable_t
    integer :: nfield, dofs, perm
    character(len=:), allocatable :: name, solver
  end type variable_t
end module parse_variable
