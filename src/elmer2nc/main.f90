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
             i,         &    ! DUMMY VARIABLE
             iostat,     &   ! status of io read
             err,        &   ! allocation status
             var,        &   ! variable index
             SavedCount, &   ! Saved timestep index
             Timestep,   &   ! Total timestep index
             ncid,       &   ! nc dataset id
             x_dimid,    &   ! x dim nc id
             x_varid,    &   ! x var nc id
             y_dimid,    &   ! y dim nc id
             y_varid,    &   ! y var nc id
             XX_varid,   &   ! X-grid var id
             YY_varid,   &   ! Y-grid var id
             pp_varid,   &   ! node number var id
             nn_varid        ! node number var id

  integer, allocatable :: &
            perm(:)!,      &  ! array to store current perm table
            !pidx(:)          ! index array for sorting arrays

  ! ! Not Sure I really Want or need this here in the Main loop
  ! integer, parameter :: maxlen=16384 ! Maximum length of character string
  ! character(len=:), allocatable :: line
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  logical  :: GotPerm         ! whether permutation read was success
  real(dp) :: Val,        &   ! Value read from .result file
              Time            ! Current time in years

  real(dp), allocatable :: X(:), & ! unique x values
                           Y(:)!, & ! unique y values
                           !Z(:)    ! unique z values

  type(node_file_t) :: parsed
  type(variable_t), dimension(:), allocatable :: variable_list

  integer, dimension(11,139) :: test3, test4
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


  test1 = (reshape(parsed%x, shape=(/Ny, Nx/), order=(/2,1/)))
  test2 = (reshape(parsed%y, shape=(/Ny, Nx/), order=(/2,1/)))
  test3 = (reshape(parsed%nn, shape=(/Ny, Nx/), order=(/2,1/)))
  test4 = (reshape(parsed%p, shape=(/Ny, Nx/), order=(/2,1/)))

  ! write(*,*) size(X)
  ! write(*,*) size(Y)

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

  ! Define Dummy X and Z variables for de-bugging
  call nc_check( nf90_def_var(ncid, "X", NF90_REAL,  (/ y_dimid, x_dimid /), XX_varid) )
  call nc_check( nf90_def_var(ncid, "Y", NF90_REAL,  (/ y_dimid, x_dimid /), YY_varid) )
  call nc_check( nf90_def_var(ncid, "NN", NF90_INT , (/ y_dimid, x_dimid /), nn_varid) )
  !call nc_check( nf90_def_var(ncid, "PP", NF90_INT , (/ y_dimid, x_dimid /), pp_varid) )

  ! add for loop to itterate over all the variables + x and z, and add them
  ! as NetCDF variables
  ! Also add the solver the variable is from as and attribute for each
  do var=5,5
    call nc_check( nf90_def_var(ncid, trim(variable_list(var)%name), NF90_REAL, &
                               (/ y_dimid, x_dimid /), variable_list(var)%nc_varid) )

     call nc_check( nf90_def_var(ncid, trim(variable_list(var)%name)//"_perm", NF90_INT, &
                               (/ y_dimid, x_dimid /), pp_varid) )
  end do
  ! End define mode.
  call nc_check( nf90_enddef(ncid) )

  ! Write the coordinate variable data.
  call nc_check( nf90_put_var(ncid, x_varid, X) )
  call nc_check( nf90_put_var(ncid, y_varid, Y) )

  ! Write the gridded coordinate data
  call nc_check( nf90_put_var(ncid, XX_varid, test1) )
  call nc_check( nf90_put_var(ncid, YY_varid, test2) )
  call nc_check( nf90_put_var(ncid, nn_varid, test3) )
  !call nc_check( nf90_put_var(ncid, pp_varid, test4) )

  ! Loop over the timesteps
  do while ( Timestep < 1 )

    ! Read the current timestep info
    call ReadTime( 10, SavedCount, Timestep, Time, iostat )

    write(*,*) SavedCount,Timestep,Time

    ! For a given timestep loop over each of the variable to be read
    do var=1,TotalDOFs

      call ReadVariableName( 10, line, iostat )

      !write(*,*) "Current Variable: ", trim(variable_list(var)%name)

      call ReadPermuation( 10, perm, GotPerm )

      variable_list(var)%perm = pack(perm, perm/=0)

      n = variable_list(var)%nfield

      !write(*,*) perm
      !
      ! write(*,*) "Length of Variable    Vector: ", variable_list(var)%nfield
      ! write(*,*) "Length of Permutation Vector: ", variable_list(var)%nperm
      ! write(*,*) "Length of Crrent Permutation Vector: ", size(perm)
      !
      ! write(*,*) "Fist Value from Perm", perm(0)

      do j = 1, variable_list(var)%nfield

        i =  variable_list(var)%nperm - variable_list(var)%nfield
        k = Perm(j)
        ! write(*,*) k
        ! if ( k == 0 ) then
        !   write(*,*) "Cycling over "
        !   cycle
        ! end if

        call ReadValue( 10, j, Val)


        if (var == 5) then
          write(*,*) j, i+j, perm(i+j), variable_list(var)%perm(j), k, Val
        end if
        ! if ( GotPerm .and. j > variable_list(var)%nfield) then
        !   write(*,*) "Cycling over ", j
        !    cycle
        ! end if
        !
        ! if (.not.GotPerm) then
        !   variable_list(var)%values(j) = Val
        ! else if ( variable_list(var)%perm(j) > 0 ) then
        !   variable_list(var)%values(variable_list(var)%perm(j)) = Val
        ! end if
      end do

    !  write(*,*) "Non-zero count:", count(perm/=0)

      if ( var == 5 ) then
        write(*,*) size(variable_list(var)%values)
        test1 = (reshape(variable_list(var)%values, shape=(/Ny, Nx/), order=(/2,1/)))
        test2 = (reshape(variable_list(var)%perm,   shape=(/Ny, Nx/), order=(/2,1/)))
        call nc_check( nf90_put_var(ncid, variable_list(var)%nc_varid, test1) )
        call nc_check( nf90_put_var(ncid, pp_varid, test2) )

        ! do k = 1, 1592
        !   write(*,*) variable_list(var)%perm(k)
        ! end do
      end if

      ! write(*,*)
      ! write(*,*) '======================='
      ! write(*,*)
    end do
  end do

  ! Close the file.
  call nc_check( nf90_close(ncid) )

  ! If we got this far, everything worked as expected. Yipee!
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"

  close(10)
  !
  ! ! Lets write the perm tables of zs and depth to dedug this perutation table problem
  ! open(12,file="results/zs_perm.dat")
  ! do j=1, size(variable_list(14)%perm)
  !   write(12,*) variable_list(14)%perm(j)
  ! end do
  ! close(12)
  !
  ! open(13,file="results/zs_vals.dat")
  ! do j=1, size(variable_list(14)%values)
  !   write(13,*) (variable_list(14)%values(j))
  ! end do
  ! close(13)
  ! ! write(*,*)  trim(variable_list(var-1)%name)
  ! ! do j = 1, 1512
  ! !   write(*,*) perm(j), variable_list(var-1)%perm(j)
  ! ! end do
  ! test1 = (reshape(variable_list(11)%perm, shape=(/Ny, Nx/),   order=(/2,1/)))
  ! test2 = (reshape(variable_list(11)%values, shape=(/Ny, Nx/), order=(/2,1/)))
  !
  !   do var=1, 1
  !     do j =1, 139
  !     write(*,*) variable_list(11)%perm(var), variable_list(11)%values(var)
  !   end do
  !   !write(*,*) variable_list(15)%perm(j), variable_list(15)%values(variable_list(15)%perm(j))
  ! end do

  ! Deallocate the arrays of the instance of our node_file type
  call parsed%deallocate_arrays()

  ! Deallocate the arrays in out list of variable
  do var = 1, TotalDOFs
    call variable_list(var)%deallocate()
  end do

end program main
