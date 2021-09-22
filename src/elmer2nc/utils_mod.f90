module utils
  use netcdf
  use iso_fortran_env
  implicit none
  ! private
  ! public :: unique_sort, dp, ReadRecur
  integer, parameter :: dp = real64
  Integer, Parameter :: kdp = selected_real_kind(15)

contains
!-------------------------------------------------------------------------------
  subroutine argparse( in_path, mesh_db, out_path, NT, Transient )
!-------------------------------------------------------------------------------
    integer, intent(out) :: NT
    logical, intent(out) :: Transient
    character(len=*), intent(inout) :: in_path, mesh_db, out_path
!-------------------------------------------------------------------------------
    character(len=1000)   ::  argument

    NT = -1; !in_path  = ' '; out_path = ' '; argument = ' '

    select case(command_argument_count())
    case(0:3)
      write(*,'(A)') "Insufficent command line arguments!"
      stop
    case(4)
      call get_command_argument(1, in_path)
      call get_command_argument(2, mesh_db)
      write(*,'(a)') trim(mesh_db)//"mesh.header"
      call get_command_argument(3, out_path)
      call get_command_argument(4, argument)
      read(argument, *) NT
      if (NT > 1) then
        transient=.TRUE.
      elseif (NT == 1) then
        transient=.FALSE.
      elseif ( NT < 0 ) then
        write(*,"(A)") "Invalid number of timesteps parsed"
        stop
      end if
    case default
      write(*,'(A)') "Too many commad line arguments passed, unable to parse."
      stop
    end select
  end subroutine argparse
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read a (logical) kine from FORTRAN unit device. Inspied by "ReadAndTrim"
! function from elmer source code (https://github.com/ElmerCSC/elmerfem/blob/6dfb482454dba8245cf35d0e1591927156b6e1ec/fem/src/GeneralUtils.F90#L842)
!-------------------------------------------------------------------------------
  recursive function ReadRecur( Unit, str ) result(l)
!-------------------------------------------------------------------------------
    integer :: Unit                       !< Fortran unit number to read from
    character(len=:), allocatable :: str  !< The string read from the file
    LOGICAL :: l                          !< Success of the read operation
!-------------------------------------------------------------------------------

    integer :: outlen ! length of output character array
    integer, parameter :: maxlen = 16384
    character(len=maxlen) :: readstr = ' '

    l = .TRUE.

    read( Unit,'(A)',end=10,err=10 ) readstr

    outlen = len(trim(readstr))
    if(.not.allocated(str)) allocate(character(512)::str)
    str(1:outlen) = trim(readstr)
    str(outlen+1:)  = ''
    return
    ! Return false (ending do while loop) when we read the end of file,
    ! or if there is an error
10  CONTINUE
    l = .FALSE.
  end function
!-------------------------------------------------------------------------------

! Helper Function for dealing with netcdf calls
!-------------------------------------------------------------------------------
  subroutine nc_check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine nc_check
!-------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  subroutine print( Caller, String, Verbose)
!------------------------------------------------------------------------------
    implicit none
    character(len=*)  :: Caller  !< The function where print is called from
    character(len=*)  :: String  !< The message to print
    logical, optional :: Verbose !< Wether the message should be printed or not
!------------------------------------------------------------------------------


  if (.not. present(Verbose)) Verbose=.FALSE.

  if ( Verbose ) then
    write(*,"(A)") trim(Caller) // ": " // trim(String)
  else
    return
  end if

  end subroutine print
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine fatal( Caller, String)
!------------------------------------------------------------------------------
    implicit none
    character(len=*)  :: Caller  !< The function where print is called from
    character(len=*)  :: String  !< The message to print
!------------------------------------------------------------------------------


  write(*,"(A)") trim(Caller) // ": " // trim(String)
  stop

  end subroutine fatal
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! ELMER FUNCTIONS DIRECTLY FROM ELMER
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Compare equality of start of s1 to (in most uses string literal) s2.
!> From elmer GeneralUtils.90
!------------------------------------------------------------------------------
  PURE FUNCTION SEQL(s1,s2) RESULT(L)
!------------------------------------------------------------------------------
    LOGICAL :: L
    CHARACTER(LEN=*), INTENT(IN) :: s1,s2
!------------------------------------------------------------------------------
    INTEGER :: n
!------------------------------------------------------------------------------
    L = .FALSE.
    n = LEN(s2)
    IF(LEN(s1) < n) RETURN
    IF (s1(1:n)==s2) L=.TRUE.
!------------------------------------------------------------------------------
  END FUNCTION SEQL
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Sort an real array, and change the order of an index array accordingly.
!------------------------------------------------------------------------------
   PURE SUBROUTINE SortD( n,a,b )
!------------------------------------------------------------------------------
     INTEGER, INTENT(in) :: n
     INTEGER, INTENT(inout) :: b(:)
     REAL(KIND=dp), INTENT(inout) :: a(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir,rb
     REAL(KIND=dp) :: ra
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN

      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

!------------------------------------------------------------------------------
   END SUBROUTINE SortD
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Code from ORDERPACK: http://www.fortran-2000.com/rank/
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
Subroutine D_unirnk (XVALT, IRNGT, NUNI)
! __________________________________________________________
!   UNIRNK = Merge-sort ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it discards indices that correspond to
!   duplicate entries.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (Kind=kdp), Dimension (:), Intent (In) :: XVALT
      Integer, Dimension (:), Intent (Out) :: IRNGT
      Integer, Intent (Out) :: NUNI
! __________________________________________________________
      Integer, Dimension (SIZE(IRNGT)) :: JWRKT
      Integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2
      Integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
      Real (Kind=kdp) :: XTST, XVALA, XVALB

!
!
      NVAL = Min (SIZE(XVALT), SIZE(IRNGT))
      NUNI = NVAL
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XVALT(IIND-1) < XVALT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XVALT(IRNGT(IWRKD+2)) <= XVALT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XVALT(IRNGT(IWRKD+1)) <= XVALT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XVALT(IRNG1) <= XVALT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XVALT(IRNG2) <= XVALT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XVALA = XVALT (JWRKT(IINDA))
            XVALB = XVALT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XVALT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XVALT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = NEARLESS (Min(XVALT(JWRKT(1)), XVALT(IRNGT(IINDB))))
      Else
        XTST = NEARLESS (XVALT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XVALT(JWRKT(IINDA)) > XVALT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XVALT(IRNG) > XTST) Then
            XTST = XVALT (IRNG)
            NUNI = NUNI + 1
            IRNGT (NUNI) = IRNG
         End If
!
      End Do
!
      Return
!
End Subroutine D_unirnk

Function nearless (XVAL) result (D_nl)
!  Nearest value less than given value
! __________________________________________________________
      Real (kind=kdp), Intent (In) :: XVAL
      Real (kind=kdp) :: D_nl
! __________________________________________________________
      D_nl = nearest (XVAL, -1.0_kdp)
      return
!
End Function nearless

end module utils
