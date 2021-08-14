module result_parser
  use parse_variable
  
  implicit none
! no neeed to re-invent the wheel, lets look at what Elmer crowd is doing:

!https://github.com/ElmerCSC/elmerfem/blob/b8b0be7c8a25d8c39b3dc03bc5a57b539b9d272f/fem/src/ModelDescription.F90#L3739
contains

  subroutine parser_results(mesh_db)
    integer, parameter :: maxlen=16384
    character(len=maxlen) :: line   = ' '
    character(len=256)    :: solver = ' ', Variable_name = ' '
    integer ierr, s, k, TotalDOFs, nlen, j, nfield, nperm, ndof
    character(len=15), intent(in) :: mesh_db


    open(10,file=mesh_db//"Accumulation_Flux.result", status='old')

    ! end do

    !------------------------------------------------------------
    ! Parse Total DOFs from header
    !------------------------------------------------------------
    ! index of string match, index intrinsic returns 0 by default
    k = 0
    ! Read through file until match found
    do while (k == 0)
      !
      read(10, '(a)', iostat=ierr) Line

      k = index( Line,'Total DOFs:',.TRUE.)

      if( k /= 0 ) then
        read(line(k+11:maxlen), *, iostat=ierr) TotalDOFs
        ! rewind the file back to the begining for more parsing
        rewind(10)
        exit
      end if
    end do


    !------------------------------------------------------------
    ! Parse Variables (and their info) from header
    !------------------------------------------------------------

    ! Pass over the first three lines, since we know the variable
    ! info starts on line four
    read(10, '(a)', iostat=ierr) Line
    read(10, '(a)', iostat=ierr) Line
    read(10, '(a)', iostat=ierr) Line

    k = 0
    do while ( k == 0 )
      read(10, '(a)', iostat=ierr) Line

      ! get the length of the trimed line
      nlen = len_trim(Line)
      ! Abort when we have reached the end of the variable list
      k = INDEX( Line(1:nlen), 'Total DOFs:',.TRUE.)
      !write(*,*) trim(Line)
      IF( k /= 0 ) EXIT

      ! Find the name of the variable
      j = index( Line(1:nlen),'[')
      if ( j /=0 ) then
        read(Line(0:j-1), *, iostat=ierr) Variable_name
      else
        k =  index( Line(1:nlen), ':')
        read(Line(1:k-1), '(A)', iostat=ierr) Variable_name
      end if

      ! read the size of filed, size of perm table, and number of dofs per node
      k = index( Line(1:nlen),']',.TRUE.)
      k = k + index( Line(k:nlen),':')
      s = k + index( Line(k:nlen),':', .TRUE.)
      read(Line(k:s-2), *, iostat=ierr) nfield, nperm, ndof
      write(*,*) nfield, nperm, ndof

      ! Find the solver the variable is from
      k = index( Line(1:nlen),']',.TRUE.)
      k = k + index( Line(k:nlen),':', .TRUE.)
      read(Line(k:nlen), *, iostat=ierr) Solver

      ! set k to zero so loop continues until the case above
      k=0
    end do

    close(10)
  end subroutine parser_results
end module result_parser
