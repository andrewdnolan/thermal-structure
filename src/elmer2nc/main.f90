program main
  use utils
  use nodes
  use result_parser
  implicit none

  integer :: Nx, Ny, i
  type(node_file_t) :: parsed
  integer, dimension(11,139) :: test3
  real(kind=dp), dimension(11,139) :: test1, test2
  character(len= *), parameter :: mesh_db  = "./results/mesh/"

  ! Parse the mesh.nodes file, which will return the x, y, z values and
  ! node indexes to be used with variable permutation tables
  !--------------------------------------------------------------------
  call parse_nodes(mesh_db, parsed)


  !n_unique =  unique_sort(parsed%x, parsed%nnodes)

  Nx = unique_sort(parsed%x, parsed%nnodes)
  Ny = unique_sort(parsed%y, parsed%nnodes)

  test1 = (reshape(parsed%x,  shape=(/Ny, Nx/), order=(/2,1/)))
  test2 = (reshape(parsed%y,  shape=(/Ny, Nx/), order=(/2,1/)))
  test3 = (reshape(parsed%nn, shape=(/Ny, Nx/), order=(/2,1/)))


  

  ! Deallocate the arrays of the instance of our node_file type
  call parsed%deallocate_arrays()

  !
  call parser_results(mesh_db)
end program main
