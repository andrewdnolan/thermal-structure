module result_parser
  use utils
  use parse_variable

  implicit none
  private
  public parser_results

  integer, parameter :: maxlen=16384
  character(len=maxlen) :: line   = ' '
! no neeed to re-invent the wheel, lets look at what Elmer crowd is doing:

!https://github.com/ElmerCSC/elmerfem/blob/b8b0be7c8a25d8c39b3dc03bc5a57b539b9d272f/fem/src/ModelDescription.F90#L3739
contains


  subroutine parser_results(mesh_db)
    integer :: TotalDOFs, err=0, idx0=23, k, n_line, psize, np, j, &
    NT, & ! local timestep count
    var   ! variable index

    integer, dimension(250) :: t_idx
    integer, dimension(:), allocatable :: p_idx, perm
    character(len=15), intent(in)   :: mesh_db
    type(variable_t), dimension(:), allocatable :: variable_list

    call parse_TotalDOFs(mesh_db, TotalDOFs)

    allocate(variable_list(TotalDOFs), stat=err)
    if ( err /= 0) print *, "variable_list: Allocation request denied"

    call parse_result_header(mesh_db, variable_list, TotalDOFs)

    call parse_result_NT(mesh_db, 250, t_idx)

    open(10,file=mesh_db//"Accumulation_Flux.result", status='old')

    k = 0
    n_line = 0
    nt = 0

    do while(k == 0)
      n_line = n_line + 1
      read(10, '(a)', iostat=err) line
      if ( any(n_line == t_idx) ) then
        nt = nt + 1 ! timestep counter
        n_line = n_line + 1
        read(10, '(a)', iostat=err) line
        ! itterate over the written variables
        do var = 1, TotalDOFs
          n_line = n_line + 1

          !write(*,*) trim(variable_list(var)%name), variable_list(var)%nfield

          j = index(line, variable_list(var)%name)
          if (j /= 0) then
            n_line = n_line + 1
            read(10, '(a)', iostat=err) line
            j = index(line, ':')

            read(line(j+1:300), *) psize, np

            allocate(perm( np))
            allocate(p_idx(np))

            do j = 1, np
              n_line = n_line + 1
              read(10, *, iostat=err) perm(j), p_idx(j)
              if (err /= 0) call abort()
            end do

            do j = 1, np
              n_line = n_line + 1
              read(10, *, iostat=err) variable_list(var)%data(perm(j), nt)
              if (err /= 0) call abort()
            end do

            do j = 1, 139
              write(*,*) variable_list(var)%data(j, nt)
            end do
          end if
        end do
        ! n_line = n_line + 1
        ! read(10, '(a)', iostat=err) line
        ! !write(*,*) trim(line)
        ! n_line = n_line + 1
        !
        ! allocate(perm( np))
        ! allocate(p_idx(np))
        !
        ! do j = 1, psize
        !   n_line = n_line + 1
        !   read(10, *) perm(j), p_idx(j)
        ! end do
        !
        ! write(*,*) n_line
        !
        ! do j = 1, psize
        !   write(*,*) perm(j), p_idx(j)
        ! end do
        ! !write(*,*) parse_perm(line, j)
        ! k = 1
        if ( nt >= 1 ) then
          exit
        end if
      end if
    end do
    close(10)
  end subroutine parser_results

  subroutine parse_TotalDOFs(mesh_db, TotalDOFs)
    implicit none
    integer :: k=0, ierr
    integer, intent(out) :: TotalDOFs
    character(len=15), intent(in) :: mesh_db

    open(10,file=mesh_db//"Accumulation_Flux.result", status='old')
    ! Read through file until match found
    do while (k == 0)
      !
      read(10, '(a)', iostat=ierr) line

      k = index( Line,'Total DOFs:',.TRUE.)

      if( k /= 0 ) then
        read(line(k+11:maxlen), *, iostat=ierr) TotalDOFs
        exit
      end if
    end do
    close(10)
  end subroutine parse_TotalDOFs


  subroutine parse_result_header(mesh_db, variable_list, TotalDOFs)
    integer :: i, ierr=0
    integer, intent(in) :: TotalDOFs
    character(len=15), intent(in)   :: mesh_db
    type(variable_t), dimension(TotalDOFs), intent(inout) :: variable_list

    character(len=256)    :: solver = ' ', Variable_name = ' '
    integer s, k, nlen, j, nfield, nperm, ndof, m

    open(10,file=mesh_db//"Accumulation_Flux.result", status='old')


    ! Pass over the first three lines, since variable info starts on line four
    do i = 1, 3
      read(10, '(a)', iostat=ierr) Line
      if (ierr /= 0) then
        write(*,*) "parse_result_header could not open file"
        stop
      end if
    end do

    k = 0
    m = 1
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

      ! Find the solver the variable is from
      k = index( Line(1:nlen),']',.TRUE.)
      k = k + index( Line(k:nlen),':', .TRUE.)
      read(Line(k:nlen), '(A)', iostat=ierr) Solver

      if ( ndof == 1 ) then
        variable_list(m)%name   = trim(Variable_name)
        variable_list(m)%solver = trim(Solver)
        variable_list(m)%dofs   = ndof
        variable_list(m)%nperm  = nperm
        variable_list(m)%nfield = nfield

        call variable_list(m)%init_data(1)
        m = m + 1
      end if

      ! set k to zero so loop continues until the case above
      k=0
    end do

    close(10)
  end subroutine parse_result_header


  subroutine parse_result_NT(mesh_db, TNT, t_idx)
    implicit none
    integer :: k
    real(kind=dp), dimension(TNT) :: t
    integer, dimension(TNT) :: idx, ns, nt
    integer, dimension(TNT), intent(out) :: t_idx
    integer, intent(in) :: TNT !total number of timesteps


    character(len=15), intent(in) :: mesh_db
    open(20,file=mesh_db//"timesteps.dat", status='old')


    do k = 1, 250
      read(20, *) idx(k), ns(k), nt(k), t(k)
    end do

    t_idx(:) = idx(:)
    close(20)
  end subroutine parse_result_NT
end module result_parser
