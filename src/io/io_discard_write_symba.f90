!**********************************************************************************************************************************
!
!  Unit Name   : io_discard_write_symba
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write out information about discarded and merged massive bodies and test particles in SyMBA
!
!  Input
!    Arguments : t             : time
!                mtiny         : smallest self-gravitating mass
!                npl           : number of massive bodies
!                nsppl         : number of spilled massive bodies
!                nsptp         : number of spilled test particles
!                nmergeadd     : number of merged massive bodies to add
!                nmergesub     : number of merged massive bodies to subtract
!                symba_pl1P    : pointer to head of SyMBA massive body structure linked-list
!                symba_pld1P   : pointer to head of discard SyMBA massive body structure linked-list
!                symba_tpd1P   : pointer to head of discard SyMBA test particle structure linked-list
!                mergeadd_list : array of structures of merged massive bodies to add
!                mergesub_list : array of structures of merged massive bodies to subtract
!                fname         : name of file to write
!                lbig_discard  : logical flag indicating whether to dump massive body data with discards
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : discard data to discard file
!
!  Invocation  : CALL io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergeadd, nmergesub, symba_pl1P, symba_pld1P,
!                                            symba_tpd1P, mergeadd_list, mergesub_list, fname, lbig_discard)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_discard_mass.f and io_discard_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, symba_plA, & 
     discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)

! Modules
     USE swiftest, EXCEPT_THIS_ONE => io_discard_write_symba
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                       :: lbig_discard
     INTEGER(I4B), INTENT(IN)                       :: npl, ntp, nsppl, nsptp, nmergeadd
     REAL(DP), INTENT(IN)                           :: t, mtiny
     CHARACTER(*), INTENT(IN)                       :: fname
     TYPE(symba_pl), INTENT(INOUT)                  :: symba_plA
     TYPE(swiftest_tp), INTENT(INOUT)               :: discard_tpA
     TYPE(swiftest_pl), INTENT(INOUT)               :: discard_plA
     TYPE(symba_merger), INTENT(INOUT)              :: mergeadd_list, mergesub_list

! Internals
     INTEGER(I4B), PARAMETER   :: LUN = 40
     INTEGER(I4B)              :: i, index, ncomp, ierr, nplm

! Executable code
     CALL io_open(LUN, fname, "APPEND", "FORMATTED", ierr)
     IF (ierr /= 0) THEN
          CALL io_open(LUN, fname, "NEW", "FORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTEST Error:"
               WRITE(*, *) "   Unable to open discard output file, ", fname
               CALL util_exit(FAILURE)
          END IF
     END IF
     WRITE(LUN, 100) t, nsppl + nsptp + 2*nmergeadd, lbig_discard 
 100 FORMAT(E23.16, 1X, I8, 1X, L1)
     index = 0
     DO i = 1, nmergeadd
          WRITE(LUN, 200) ADD, mergeadd_list%name(i), mergeadd_list%status(i)
 200      FORMAT(A, 2(1X, I8))
          WRITE(LUN, 300) mergeadd_list%xh(:,i)
 300      FORMAT(3(E23.16, 1X))
          WRITE(LUN, 300) mergeadd_list%vh(:,i)
          ncomp = mergeadd_list%ncomp(i)
          DO index = 1, ncomp !j=0 -> index=1
               !index = index + 1
               WRITE(LUN, 200) SUB, mergesub_list%name(index), mergesub_list%status(index)
               WRITE(LUN, 300) mergesub_list%xh(:,index)
               WRITE(LUN, 300) mergesub_list%vh(:,index)
               WRITE(LUN, 500) mergesub_list%name(index), mergesub_list%mass(index), mergesub_list%radius(index)
          END DO
     END DO
     DO i = 1, nsppl
          IF (discard_plA%status(i) /= MERGED) THEN
               WRITE(LUN, 200) SUB, discard_plA%name(i), discard_plA%status(i)
               WRITE(LUN, 300) discard_plA%xh(1,i),discard_plA%xh(2,i), discard_plA%xh(3,i)
               WRITE(LUN, 300) discard_plA%vh(1,i),discard_plA%vh(2,i), discard_plA%vh(3,i)
               WRITE(LUN, 500) discard_plA%name(i),discard_plA%mass(i), discard_plA%radius(i)
          END IF
     END DO
     DO i = 1, nsptp
          WRITE(LUN, 200) SUB, discard_tpA%name(i), discard_tpA%status(i)
          WRITE(LUN, 300) discard_tpA%xh(1,i),discard_tpA%xh(2,i), discard_tpA%xh(3,i)
          WRITE(LUN, 300) discard_tpA%vh(1,i),discard_tpA%vh(2,i), discard_tpA%vh(3,i)
     END DO
     IF (lbig_discard) THEN
          nplm = 0
          DO i = 1, npl
               IF (symba_plA%mass(i) < mtiny) EXIT
               nplm = nplm + 1
          END DO
          IF (nplm > 1) THEN
               WRITE(LUN, 400) nplm
 400           FORMAT(I8)
               DO i = 2, nplm
                    WRITE(LUN, 500) symba_plA%name(i), symba_plA%mass(i),& 
                     symba_plA%radius(i)
 500                FORMAT(I8, 2(1X, E23.16))
                    WRITE(LUN, 300) symba_plA%xh(:,i)
                    WRITE(LUN, 300) symba_plA%vh(:,i)
               END DO
          END IF
     END IF
     CLOSE(LUN)

     RETURN

END SUBROUTINE io_discard_write_symba
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
