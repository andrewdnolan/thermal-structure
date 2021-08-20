program main

  use utils
  use nodes
  use netcdf
  use result_parser
  implicit none

  integer :: Nx,         &   ! Number of horizontal nodes
             Ny,         &   ! Number of vertical   nodes
             TotalDOFs,  &   ! Number of variables to be saved
             j,          &   ! counter to loop over nodes with
             n,          &   ! number of nodes to loop over
             k,          &   ! permutation index
             nt,         &   ! local timestep count
             iostat,     &   ! status of io read
             err,        &   ! allocation status
             var,        &   ! variable index
             SavedCount, &   ! Saved timestep index
             Timestep,   &   ! Total timestep index
             ncid,       &   ! nc dataset id
             x_dimid,    &   ! x dim nc id
             x_varid,    &   ! x var nc id
             y_dimid,    &   ! y dim nc id
             y_varid         ! y var nc id

  integer, allocatable :: &
            perm(:),      &  ! array to store current perm table
            pidx(:)          ! index array for sorting arrays

  ! ! Not Sure I really Want or need this here in the Main loop
  ! integer, parameter :: maxlen=16384 ! Maximum length of character string
  ! character(len=:), allocatable :: line
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  logical  :: GotPerm         ! whether permutation read was success
  real(dp) :: Val,        &   ! Value read from .result file
              Time            ! Current time in years

  real(dp), allocatable :: X(:), & ! unique x values
                           Y(:), & ! unique y values
                           Z(:)    ! unique z values

  type(node_file_t) :: parsed
  type(variable_t), dimension(:), allocatable :: variable_list

  integer, dimension(11,139) :: test3
  real(kind=dp), dimension(11,139) :: test1, test2
  character(len= *), parameter :: mesh_db  = "./results/mesh/"



  ! Parse the mesh.nodes file, which will return the x, y, z values and
  ! node indexes to be used with variable permutation tables
  !--------------------------------------------------------------------
  call parse_nodes(mesh_db, parsed)

  ! Allocate dummy array while I figure out what the hell is up with the
  ! unique / sort subroutines
  allocate(X(139), source=parsed%x(1:139), stat=err)
  if ( err /= 0) print *, "variable_list: Allocation request denied"

  allocate(Y(11), source=parsed%y(1:1529:139), stat=err)
  if ( err /= 0) print *, "variable_list: Allocation request denied"

  Nx = size(X)
  Ny = size(Y)
  write(*,*) size(X)
  write(*,*) size(Y)

  ! Open the input file from which we will parse all our data
  open(10,file=mesh_db//"Accumulation_Flux.result", status='old')

  ! Find the number of variables to be parsed
  call ReadTotalDOFs(10, TotalDOFs, iostat)

  ! Allocate an array of our derived type to store the variable info
  allocate(variable_list(TotalDOFs), stat=err)
  if ( err /= 0) print *, "variable_list: Allocation request denied"

  ! Parse the variale info the the header of the .result file and store in
  ! the array of our derived types
  call ReadResultHeader(10, variable_list, TotalDOFs)

  ! Now that all the header info has been parsed, lets create our NetCDF file
  ! where we will store all out data.
  !--------------------------------------------------------------------------

  ! Create the NetCDF file.
  call nc_check( nf90_create("results/test.nc", nf90_clobber, ncid) )

  ! Define the dimensions.
  call nc_check( nf90_def_dim(ncid, "coord_1", Nx, x_dimid) )
  call nc_check( nf90_def_dim(ncid, "coord_2", Ny, y_dimid) )

  ! Define the coordinate variables.
  call nc_check( nf90_def_var(ncid, "coord_1", NF90_REAL, x_dimid, x_varid) )
  call nc_check( nf90_def_var(ncid, "coord_2", NF90_REAL, y_dimid, y_varid) )

  ! add for loop to itterate over all the variables + x and z, and add them
  ! as NetCDF variables

  ! End define mode.
  call nc_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data.
  call nc_check( nf90_put_var(ncid, x_varid, X) )
  call nc_check( nf90_put_var(ncid, y_varid, Y) )

  ! Close the file.
  call nc_check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee!
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"
  ! Loop over the timesteps
  do while ( Timestep < 10 )

    ! Read the current timestep info
    call ReadTime( 10, SavedCount, Timestep, Time, iostat )

    !write(*,*) SavedCount,Timestep,Time

    ! For a given timestep loop over each of the variable to be read
    do var=1,TotalDOFs
      ! if ( Timestep < 2 ) then
      !   write(*,*) variable_list(var)%nc_varid
      ! end if

      call ReadVariableName( 10, line, iostat )

      write(*,*) trim(variable_list(var)%name) == trim(line)

      call ReadPermuation( 10, perm, GotPerm )

      n = variable_list(var)%nfield

      do j = 1, n

        call ReadValue( 10, perm, j, k, Val)

        variable_list(var)%perm(j) = k
        variable_list(var)%values(k) = Val

        !
      end do
    end do
  end do

  close(10)
  !
  ! ! Lets write the perm tables of zs and depth to dedug this perutation table problem
  open(12,file="results/zs_perm.dat")
  do j=1, size(variable_list(14)%perm)
    write(12,*) variable_list(14)%perm(j)
  end do
  close(12)

  open(13,file="results/zs_vals.dat")
  do j=1, size(variable_list(14)%values)
    write(13,*) (variable_list(14)%values(j))
  end do
  close(13)
  ! ! write(*,*)  trim(variable_list(var-1)%name)
  ! ! do j = 1, 1512
  ! !   write(*,*) perm(j), variable_list(var-1)%perm(j)
  ! ! end do
  ! test1 = (reshape(variable_list(11)%perm, shape=(/Ny, Nx/),   order=(/2,1/)))
  ! test2 = (reshape(variable_list(11)%values, shape=(/Ny, Nx/), order=(/2,1/)))
  !
    do var=1, 1
      do j =1, 139
      write(*,*) variable_list(11)%perm(var), variable_list(11)%values(var)
    end do
    !write(*,*) variable_list(15)%perm(j), variable_list(15)%values(variable_list(15)%perm(j))
  end do

  ! Deallocate the arrays of the instance of our node_file type
  call parsed%deallocate_arrays()

  ! Deallocate the arrays in out list of variable
  do var = 1, TotalDOFs
    call variable_list(var)%deallocate()
  end do

end program main
