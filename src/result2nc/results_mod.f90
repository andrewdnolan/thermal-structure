!-------------------------------------------------------------------------------
! results_mod.f90
!
! A collection of subroutines to parse the various parts of the .result file.
! These functions were heavily inspired/influenced by the elmer source code, which
! can be found here:
! https://github.com/ElmerCSC/elmerfem/blob/devel/fem/src/ModelDescription.F90
!-------------------------------------------------------------------------------
module result_parser
  use utils
  use parse_variable

  implicit none

  integer, parameter :: maxlen=16384
  character(len=:), allocatable :: line

contains

!-------------------------------------------------------------------------------
! Read the TotalDOFs (i.e. number of variables with DOFs of 1) from the end
! of the .result header. Then rewind open file for parsing of reading individual
! varibale info.
!-------------------------------------------------------------------------------
  subroutine ReadTotalDOFs(RestartUnit, TotalDOFs, Stat)
    implicit none
    integer,  intent(in) :: RestartUnit   !< Fortran unit number to read from
    integer, intent(out) :: TotalDOFs     !< Total # of variables in .result file
    integer, intent(out) :: Stat          !< Status of io read
!-------------------------------------------------------------------------------
    integer :: k

    do while (  ReadRecur(RestartUnit, line) )
      k = index( Line,'Total DOFs:',.TRUE.)

      if( k /= 0 ) then
        read(line(k+11:), *, iostat=Stat) TotalDOFs
        if ( Stat /= 0) then
          write(*,'(A)') "Error: "//"Error in ReadDOFs"
          stop
        end if
        rewind(RestartUnit)
        exit
      end if
    end do
  end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Inspired by ReadTime from ModelDescription.f90 in elmer source
! Read SavedCount, Timestep and Time.  Stat is set to 0 for success, < 0 for
! end-of-file, and > 0 for error.
!-------------------------------------------------------------------------------
  subroutine ReadTime( RestartUnit,SavedCount,Timestep,Time, Stat)
    implicit none
    integer, intent(in)  :: RestartUnit      !< Fortran unit number to read from
    integer, intent(out) :: Stat             !< Status of io read
    integer, intent(out) :: SavedCount       !< Saved timestep
    integer, intent(out) :: Timestep         !< Total timestep
    real(kind=dp), intent(out):: Time        !< Physical Time
!-------------------------------------------------------------------------------
    integer :: iostat

    do while (  ReadRecur(RestartUnit, line) )
      if (SEQL(line, 'Time:')) then
        read( line(7:), *, iostat=iostat) SavedCount, Timestep, Time
        if ( iostat /= 0 ) then
          write(*,'(A)') "Error: Error in ReadTime!"
          STOP
        end if
        stat = 0
        return
      end if
    end do
    stat = -1
  end subroutine ReadTime
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read Variable name from the open .result file
!-------------------------------------------------------------------------------
  subroutine ReadVariableName( RestartUnit, VarName, iostat )
    implicit none
    integer,       intent(in) :: RestartUnit !< Fortran unit number to read from
    integer,      intent(out) :: iostat      !< Status of io read
    character(*), intent(out) :: VarName     !< Variable Name
!-------------------------------------------------------------------------------
    character(*), parameter :: Caller="ReadVariableName"

    read(RestartUnit, '(a)', iostat=iostat) VarName
    if (iostat /=0 ) then
      call fatal(Caller, "Error reading Variable Name ")
    end if
  end subroutine ReadVariableName
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read permutation for each variable. Subroutine handles "use previous" scenarios
! with no problem, but does not handle "NULL" scenario correctly. In practice, I've
! never come across this, but could be an avenue for future development.
!
! Currently returns the permutation INDEXES not the node numbers, as in the elmer
! src code. Indexes seem to work, and I was unable to figure out how to use node
! number succesfully.
!
! This function is high on the list for more error catching / handling.
!-------------------------------------------------------------------------------
  subroutine ReadPermuation( RestartUnit, Perm, GotPerm)
    implicit none
    integer, intent(in)  :: RestartUnit      !< Fortran unit number to read from
    logical, intent(out) :: GotPerm          !< true is succesfully read perm
    integer, allocatable :: Perm(:)          !< Where to store the perm table
!-------------------------------------------------------------------------------
    integer :: &
    nPerm,     & !< Length of permutation table
    nPositive, & !< Number of postive nodes ?
    i,         & !< Counter in read Perm loop
    j, k,      & !< index and node number (resectively) from perm table
    iostat       !< Status of the variable read

    character(*), parameter :: Caller="ReadPermuation"

    ! Read the line defining the permutation info
    read(RestartUnit, '(A)') line
    if ( line(7:10) == "NULL" ) then
      nPerm = 0
    else if (line(7:18) == "use previous") then
      nPerm = -1
    else
      read(line(7:), *, iostat=iostat) nPerm, nPositive
      if ( iostat /=0 ) then

        call fatal(Caller, "Error reading sizes in ReadPermuation:  "//trim(line))
      end if
    end if

    ! Special cases where permutation is not directly provided
    if (nPerm < 0) then
      !TO DO: Need to pass Verbose as variable not hard code it
      call print(Caller, 'Using pervious permutation table', .FALSE.)
      GotPerm=.TRUE.
      return
    else if (nPerm == 0) then
      ! TO DO: Need to give some kind of warning about whats going on here
      write(*, "(A)") "Fuck"
      return
    end if

    ! If the Perm vector is already allocated make sure its the correct length
    if (allocated(Perm)) then
      if ( size(Perm) < nPerm ) then
        call print(Caller, 'Permutation vector too small??', .TRUE.)
        deallocate(Perm)
      else if ( size(Perm) > nPerm ) then
        call print(Caller, 'Permutation vector too big??', .TRUE.)
        deallocate(Perm)
      end if
    end if

    ! If Perm vector is not allocated, then allocate
    if( .not.allocated(Perm) ) allocate(Perm(nPerm))

    ! Initialize the Perm vector, or wipe all pervious data
    Perm = 0

    !Actually read the permutation table
    do i = 1, nPositive
      read(RestartUnit, *, iostat=iostat) j, k
      if ( iostat /= 0 ) then
        call fatal(Caller, 'Error reading values in ReadPermuation')
      end if
      !Perm(j) = k
      Perm(i) = j
    end do

    GotPerm = .TRUE.
  end subroutine ReadPermuation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read numeric value line by line for a given variable.
!-------------------------------------------------------------------------------
  subroutine ReadValue( RestartUnit, iNode, Val) !, Perm, iPerm,
    implicit none
    integer, intent(in)   :: RestartUnit      !< Fortran unit number to read from
    integer, intent(in)   :: iNode            !< Node index
    real(dp), intent(out) :: Val              !< Value read from the .result file
!-------------------------------------------------------------------------------
    integer :: iostat
    character(*), parameter :: Caller="ReadValue"

    read(RestartUnit, *, iostat=iostat) Val
    if (iostat /= 0) then
      print '("Error at Line ", I0)', iNode
      call fatal(Caller, "Error in ReadValue for ")
    end if

  end subroutine ReadValue
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Read/Parse the header of the .result file. Here we parse the variable name,
! the size of the field, size of the permutation table, the DOFs per variable,
! and the solver the variable is from.
!-------------------------------------------------------------------------------
  subroutine ReadResultHeader(RestartUnit, variable_list, TotalDOFs)
!-------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: RestartUnit      !< Fortran unit number to read from
    integer, intent(in) :: TotalDOFs        !< Total number of variable to be read
    !< An array of our variable_t class to store variable info in
    type(variable_t), dimension(TotalDOFs), intent(inout) :: variable_list
!-------------------------------------------------------------------------------
    integer :: ierr=0
    character(len=256)    :: solver = ' ', Variable_name = ' '
    integer s, k, nlen, j, nfield, nperm, ndof, m

    m = 1
    do while ( ReadRecur(RestartUnit, line) )

      ! get the length of the trimed line
      nlen = len_trim(Line)
      ! Abort when we have reached the end of the variable list
      k = INDEX( Line(1:nlen), 'Total DOFs:',.TRUE.)
      IF( k /= 0 ) EXIT

      ! Find the name of the variable
      !------------------------------
      j = index( Line(1:nlen),'[')
      if ( j /=0 ) then
        read(Line(0:j-1), *, iostat=ierr) Variable_name
      else
        k =  index( Line(1:nlen), ':')
        read(Line(1:k-1), '(A)', iostat=ierr) Variable_name
      end if

      ! read the size of filed, size of perm table, and number of dofs per node
      !------------------------------------------------------------------------
      k = index( Line(1:nlen),']',.TRUE.)
      k = k + index( Line(k:nlen),':')
      s = k + index( Line(k:nlen),':', .TRUE.)
      read(Line(k:s-2), *, iostat=ierr) nfield, nperm, ndof

      ! Find the solver the variable is from
      !-------------------------------------
      k = index( Line(1:nlen),']',.TRUE.)
      k = k + index( Line(k:nlen),':', .TRUE.)
      read(Line(k:nlen), '(A)', iostat=ierr) Solver

      ! Initialize the mth instance of our class in the inout array
      !------------------------------------------------------------
      if ( ndof == 1 ) then
        !allocate(variable_list(m), mold=nfield)

        variable_list(m)%name   = trim(Variable_name)
        variable_list(m)%solver = trim(Solver)
        variable_list(m)%dofs   = ndof
        variable_list(m)%nperm  = nperm
        variable_list(m)%nfield = nfield

        ! Allocate the values and perm arrays
        call variable_list(m)%allocate()

        m = m + 1
      end if
    end do

  end subroutine ReadResultHeader
!-------------------------------------------------------------------------------

end module result_parser
