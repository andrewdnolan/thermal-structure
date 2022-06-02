module nodes
  use utils
  implicit none
  private
  public node_file_t, parse_nodes

  type :: node_file_t
     integer :: nnodes
     integer,       dimension(:), allocatable :: nn, p
     real(kind=dp), dimension(:), allocatable :: x, y, z
  contains
    procedure ::   allocate_arrays
    procedure :: deallocate_arrays
  end type node_file_t

contains
  ! procedures (i.e methods) for node_file type
  subroutine allocate_arrays(self, nnodes)
    integer, intent(in) :: nnodes
    class(node_file_t), intent(inout) :: self

    allocate(self%nn(nnodes))
    allocate(self %p(nnodes))
    allocate(self %x(nnodes))
    allocate(self %y(nnodes))
    allocate(self %z(nnodes))

  end subroutine allocate_arrays

  subroutine deallocate_arrays(self)
    class(node_file_t) :: self

    deallocate(self%nn)
    deallocate(self %p)
    deallocate(self %x)
    deallocate(self %y)
    deallocate(self %z)

  end subroutine deallocate_arrays
  !------------------------------------------

  !-----------------------------------------------------------------------------
  ! Public function to parse the mesh.header and mesh.nodes files
  !
  subroutine parse_nodes(mesh_db, parsed)

    implicit none
    integer  :: nnodes
    character(len=2000), intent(in) :: mesh_db
    type(node_file_t), intent(inout) :: parsed

    ! Read the header file to find the total number of nodes
    call parse_header_file(mesh_db, nnodes)

    !Allocate the arrays of the instance of our node_file type
    call parsed%allocate_arrays(nnodes)

    ! parse mesh.nodes file and store the data with our derived type
    call parse_node_file(mesh_db, nnodes, parsed %nn, parsed %p, parsed %x, parsed %y, parsed %z)

    parsed%nnodes=nnodes
  end subroutine parse_nodes


  subroutine parse_node_file(mesh_db, nnodes, nn, p, x, y, z)

    implicit none

    integer :: nnodes, stat, i
    integer, parameter :: uid=11
    integer, dimension(nnodes) :: nn, p
    character(len=2000), intent(in)  :: mesh_db
    real(kind=dp), dimension(nnodes), intent(out) :: x, y, z

    open(uid, file=trim(mesh_db)//"mesh.nodes", status='old', iostat=stat)
    if (stat /= 0) then
      write(*,'(a)') "Error, could not open"//mesh_db//"mesh.nodes"
      stop
    end if
    do i = 1, nnodes
      read(uid,*) nn(i), p(i), x(i), y(i), z(i)
    end do
    close(uid)
  end subroutine parse_node_file

  subroutine parse_header_file(mesh_db, nnodes)
    !---------------------------------------------------------------------------
    ! parse the mesh.header and find the number of nodes
    !---------------------------------------------------------------------------

    implicit none

    integer, parameter :: uid=10
    integer, intent(out) :: nnodes
    character(len=2000), intent(in) :: mesh_db
    integer :: elements, boundary_elements, stat

    open(uid, file = trim(mesh_db)//"mesh.header", status='old', iostat=stat)
    if (stat /= 0) then
      write(*,'(a)') "Error, could not open "//trim(mesh_db)//"mesh.header"
      stop
    end if
    read(uid, *) nnodes, elements, boundary_elements
    close(uid)
  end subroutine parse_header_file


end module nodes
