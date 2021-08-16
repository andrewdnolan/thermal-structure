program main
  use nodes
  use result_parser
  implicit none

  type(node_file) :: parsed
  character(len= *), parameter :: mesh_db  = "./results/mesh/"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parse the mesh.nodes file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call parse(mesh_db, parsed)

  ! write(*,*) parsed%x(100)
  ! write(*,*)
  ! write(*,*) shape(parsed%x)
  ! write(*,*)

  ! Deallocate the arrays of the instance of our node_file type
  call parsed%deallocate_arrays()

  !
  call parser_results(mesh_db)
end program main
