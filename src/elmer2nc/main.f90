program main

  use utils
  use nodes
  use netcdf
  use result_parser
  implicit none

  integer :: Nx,         &   ! Number of horizontal nodes
             Nz,         &   ! Number of vertical   nodes
             TotalDOFs,  &   ! Number of variables to be saved
             j,          &   ! counter to loop over nodes with
             NT,         &   ! number of timesteps to loop over
             k,          &   ! permutation index
             i,          &   ! counter to loop over variables
             iostat,     &   ! status of io read
             err,        &   ! allocation status
             SavedCount, &   ! Saved timestep index
             Timestep,   &   ! Total timestep index
             ndims,      &   ! number of dimensions of the data
             ncid,       &   ! nc dataset id
             x_dimid,    &   ! x dim nc id
             x_varid,    &   ! x var nc id
             z_dimid,    &   ! z dim nc id
             z_varid,    &   ! z var nc id
             t_dimid,    &   ! t dim nc id
             t_varid,    &   ! t var nc id
             XX_varid,   &   ! X-grid var id
             ZZ_varid,   &   ! Z-grid var id
             nn_varid        ! node number var id

  integer, allocatable :: &
            perm(:),      &  ! array to store current perm table
            dimids(:),    &  ! array of netcdf ids for each dimension
            pidx(:),      &  ! index array for sorting finding unique coords
            nodes_num(:,:)   ! gridded array node number

  logical  :: GotPerm,            &   ! whether permutation read was success
              transient,          &   ! whether data is time dependent
              bool_2D=.FALSE.,    &   ! whether data 2D
              bool_3D=.FALSE.         ! whether data 3D

  real(dp) :: Val,        &   ! Value read from .result file
              Time            ! Current time in years

  real(dp), allocatable :: X(:),        & ! unique x values
                           Z(:),        & ! unique z values
                           temp(:),     & ! temporary array to fill in with values from lower dim. fields
                           vals(:,:),   & ! gridded variable values to write to NetCDF
                           coord_1(:,:),& ! girdded x-coordinate
                           coord_2(:,:)   ! girdded Z-coordinate

  type(node_file_t) :: parsed
  type(variable_t), pointer :: var
  type(variable_t), allocatable, target :: variable_list(:)

  character(len= *), parameter :: mesh_db  = "./results/mesh/"
  character(len= *), parameter :: caller = "main.f90"
  character(len= *), parameter :: SOLVER = "Solver" ! The solver the variable was written by


  ! TODO: Manually setting transient to true, but should be parsed from cmd line
  transient = .TRUE.
  NT = 250


  ! Parse the mesh.nodes file, which will return the x, y, z values and
  ! node indexes to be used with variable permutation tables
  !--------------------------------------------------------------------
  call parse_nodes(mesh_db, parsed)

  ! Allocate array of unique indexes, only need one will be written over each call
  allocate(pidx(size(parsed%x(:))), source= (/(i,i=1,size(parsed%x(:)))/))

  ! Call unique rank from ORDERPACK 2.0 for x coordinate
  call D_unirnk(parsed%x(:), pidx, Nx)
  allocate(X(Nx), source=(/(parsed%x(pidx(i)),i=1,Nx)/), stat=err)
  if (err /= 0) print *, "X(Nx): Allocation request denied"

  ! Call unique rank from ORDERPACK 2.0 for z coordinate
  call D_unirnk(parsed%y(:), pidx, Nz)
  allocate(Z(Nz), source=(/(parsed%y(pidx(i)),i=1,Nz)/), stat=err)
  if (err /= 0) print *, "Z(Nz): Allocation request denied"

  ! Allocate the data/coordinate grids
  allocate(coord_1(Nz,Nx),   stat=err)
  if (err /= 0) print *, "coord_1(Nz, Nx): Allocation request denied"
  allocate(coord_2(Nz,Nx),   stat=err)
  if (err /= 0) print *, "coord_2(Nz,Nx): Allocation request denied"
  allocate(nodes_num(Nz,Nx), stat=err)
  if (err /= 0) print *, "nodes_num(Nz,Nx): Allocation request denied"

  ! Reshape the x,y corrdinates and node number to the grid
  coord_1   = (reshape(parsed%x,  shape=(/Nz, Nx/), order=(/2,1/)))
  coord_2   = (reshape(parsed%y,  shape=(/Nz, Nx/), order=(/2,1/)))
  nodes_num = (reshape(parsed%nn, shape=(/Nz, Nx/), order=(/2,1/)))

  ! 2D Case. Primary use case
  !------------------------------------------------
  if (count(parsed%z - 0.0 >= 1e-6) == 0) then
    bool_2D = .TRUE.

    ! Set the number of dimensions
    if (transient) then; ndims=3; else; ndims=2; endif

    ! Allocate the dimids array
    allocate(dimids(ndims), stat=err)
    if ( err /= 0) print *, "ndims: Allocation request denied"

  ! 3D Case. Experimentally functional
  !------------------------------------------------
  else
    bool_3D = .TRUE.

    ! Set the number of dimensions
    if (transient) then; ndims=4; else; ndims=3; endif

    ! Allocate the dimids array
    allocate(dimids(ndims), stat=err)
    if ( err /= 0) print *, "ndims: Allocation request denied"

  endif

  ! Only going to support 2D for now, haven't tested 3D enough to use it
  !---------------------------------------------------------------------
  if (.not.bool_2D) then
    call fatal(caller, "Only 2D or 3D data is currently supported")
  end if

  ! test3 = (reshape(parsed%nn, shape=(/Nz, Nx/), order=(/2,1/)))
  ! test4 = (reshape(parsed%p,  shape=(/Nz, Nx/), order=(/2,1/)))

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
  ! TODO: Take variable whih is the outfile name/path
  call nc_check( nf90_create("results/test.nc", nf90_clobber, ncid) )

  ! Define the dimensions.
  call nc_check( nf90_def_dim(ncid, "coord_1", Nx, x_dimid) )
  call nc_check( nf90_def_dim(ncid, "coord_2", Nz, z_dimid) )
  ! if the data is time dependent, then add (unlimited) time dimension
  if (transient) then
    call nc_check( nf90_def_dim(ncid, "t", NF90_UNLIMITED , t_dimid) )
    dimids=(/z_dimid, x_dimid, t_dimid/)
  else
    dimids=(/z_dimid, x_dimid/)
  end if

  ! Define the coordinate variables.
  call nc_check( nf90_def_var(ncid, "coord_1", NF90_REAL, x_dimid, x_varid) )
  call nc_check( nf90_def_var(ncid, "coord_2", NF90_REAL, z_dimid, z_varid) )

  ! Define X, Z, and node number (NN) variables
  call nc_check( nf90_def_var(ncid, "X", NF90_REAL,  (/ z_dimid, x_dimid /), XX_varid) )
  call nc_check( nf90_def_var(ncid, "Z", NF90_REAL,  (/ z_dimid, x_dimid /), ZZ_varid) )
  call nc_check( nf90_def_var(ncid, "NN", NF90_INT , (/ z_dimid, x_dimid /), nn_varid) )
  ! Define the time variable
  if (transient) then
    call nc_check( nf90_def_var(ncid, "t", NF90_REAL, t_dimid, t_varid))
  end if

  ! Initialize each variable and it's attributes from the header within the ouput file
  do i=1,TotalDOFs
    ! Initialize pointer to variable_t instance
    var => variable_list(i)
    ! Create the variable
    call nc_check( nf90_def_var(ncid, trim(var%name), NF90_REAL, dimids, var%nc_varid) )
    ! Add Solver attribute to each varibale
    call nc_check( nf90_put_att(ncid, var%nc_varid, SOLVER, trim(var%solver) ) )
  end do

  ! End define mode.
  call nc_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data.
  call nc_check( nf90_put_var(ncid, x_varid, X) )
  call nc_check( nf90_put_var(ncid, z_varid, Z) )

  ! Write the gridded coordinate data
  call nc_check( nf90_put_var(ncid, XX_varid, coord_1) )
  call nc_check( nf90_put_var(ncid, ZZ_varid, coord_2) )
  call nc_check( nf90_put_var(ncid, nn_varid, nodes_num) )

  ! Loop over the timesteps
  do while ( Timestep < NT   )

    ! Read the current timestep info
    call ReadTime( 10, SavedCount, Timestep, Time, iostat )

    !write(*,*) SavedCount,Timestep,Time

    ! For a given timestep loop over each of the variable to be read
    do i=1,TotalDOFs

      ! Initialize pointer to variable_t instance
      var => variable_list(i)

      ! Read the variable name
      call ReadVariableName( 10, line, iostat )

      ! Read the permutation for the variable
      call ReadPermuation( 10, perm, GotPerm )

      ! Store the compressed permuation table within our derived type.
      ! We only keep the non-zero node indexes, this is neccessary when lower dim
      ! bodies are present (e.g. the free surface.)
      var%perm = pack(perm, perm/=0)

      ! Loop over the variable values
      do j = 1, var%nfield

        ! if permuation read was succesfull use indexes from perutation table
        if (GotPerm) k = Perm(j)

        ! actually read the variable value
        call ReadValue( 10, j, Val)

        !-----------------------------------------------------------------------
        ! Deal with permutation table error/exceptions and actually write the
        ! data to the allocated array in out derived type
        !-----------------------------------------------------------------------
        if (.not.GotPerm) then

          call fatal("ReadPermuation", &
          "Error somethings gone wrong reading permutation table for "//trim(var%name) )

        ! If perm is non-zero, use indexes from perm to asign variable values
        ! (This condition should always be true is steps above worked properly)
        else if ( var%perm(j) > 0 ) then
          var%values(var%perm(j)) = Val

        ! If perm is full of zeros, then just index from loop to assign values
        ! TODO: this will likely cause problems and should be looked into more
        else

          print *, "Warning, while the permutation table was read successully \n &
                    &the values were all zero, using loop indexes to assign values"

          var%perm(j) = k
          var%values(k) = Val
        end if
      end do

      ! things most likely are already sorted, but just incase lets sort nodes indexes
      ! in increasing order
      ! ALLOCATE (ind(Ns))
      ! ind = (/(i,i=1,Ns)/)
      ! CALL SortD(Ns,xs,ind)

      !-------------------------------------------------------------------------
      ! Deal with lower dimmmentional data (e.g. free surface) and ensure it's
      ! compatiable to be written to the open NetCDF file
      !-------------------------------------------------------------------------
      if (size(var%values) < size(parsed%x)) then

        ! Allocate temp array, or wipe and fill with zeros if neccessary
        if (allocated(temp)) then
          temp = 0.0
        else
          allocate(temp(size(parsed%x)), stat=err)
          if (err /= 0)  print *, "temp: Allocation request denied"
        end if

        ! use the permutation indexes to fill the temp array
        do j=1,size(var%values)
          temp(var%perm(j)) = var%values(var%perm(j))
        end do

        ! reshape the data to the 2D grid
        vals = (reshape(temp, shape=(/Nz, Nx/), order=(/2,1/)))
      else
        ! reshape the data to the 2D grid
        vals = (reshape(var%values, shape=(/Nz, Nx/), order=(/2,1/)))
      end if

      !-------------------------------------------------------------------------
      ! Actually do NetCDF writing with the parsed and sorted data
      !-------------------------------------------------------------------------
      ! write the variable data
      call nc_check( nf90_put_var(ncid, var%nc_varid, vals,    &
                                  start=(/1, 1, SavedCount /), &
                                  count=(/Nz, Nx, 1/) ) )
      ! write the timestep info
      call nc_check( nf90_put_var(ncid, t_varid, Time, start=(/SavedCount/) ) )

    end do
  end do

  ! Close the file.
  call nc_check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee!
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"

  close(10)

  ! Deallocate the arrays of the instance of our node_file type
  call parsed%deallocate_arrays()

  ! Deallocate the arrays in out list of variable
  do i = 1, TotalDOFs
    ! Initialize pointer to variable_t instance
    var => variable_list(i)
    ! Deallocate the variable values and perm arrays
    call var%deallocate()
  end do

end program main
