!**********************************************************************************************************************************
!
!  Unit Name   : io_discard_write_symba
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write out information about discarded and merged planets and test particles in SyMBA
!
!  Input
!    Arguments : t             : time
!                mtiny         : smallest self-gravitating mass
!                npl           : number of planets
!                nsppl         : number of spilled planets
!                nsptp         : number of spilled test particles
!                nmergeadd     : number of merged planets to add
!                nmergesub     : number of merged planets to subtract
!                symba_pl1P    : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P   : pointer to head of discard SyMBA planet structure linked-list
!                symba_tpd1P   : pointer to head of discard SyMBA test particle structure linked-list
!                mergeadd_list : array of structures of merged planets to add
!                mergesub_list : array of structures of merged planets to subtract
!                fname         : name of file to write
!                lbig_discard  : logical flag indicating whether to dump planet data with discards
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
SUBROUTINE io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergeadd, nmergesub, symba_pl1P, symba_pld1P, symba_tpd1P,        &
     mergeadd_list, mergesub_list, fname, lbig_discard)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => io_discard_write_symba
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                     :: lbig_discard
     INTEGER(I4B), INTENT(IN)                     :: npl, nsppl, nsptp, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                         :: t, mtiny
     CHARACTER(*), INTENT(IN)                     :: fname
     TYPE(symba_pl), POINTER                      :: symba_pl1P, symba_pld1P
     TYPE(symba_tp), POINTER                      :: symba_tpd1P
     TYPE(symba_merger), DIMENSION(:), INTENT(IN) :: mergeadd_list, mergesub_list

! Internals
     INTEGER(I4B), PARAMETER   :: LUN = 40
     INTEGER(I4B)              :: i, index, j, ncomp, ierr, nplm
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(symba_pl), POINTER   :: symba_plP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     CALL io_open(LUN, fname, "APPEND", "FORMATTED", ierr)
     IF (ierr /= 0) THEN
          CALL io_open(LUN, fname, "NEW", "FORMATTED", ierr)
          IF (ierr /= 0) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Unable to open discard output file, ", fname
               CALL util_exit(FAILURE)
          END IF
     END IF
     WRITE(LUN, 100) t, nsppl + 2*nmergeadd + nsptp, lbig_discard
 100 FORMAT(E23.16, 1X, I8, 1X, L1)
     index = 0
     DO i = 1, nmergeadd
          WRITE(LUN, 200) ADD, mergeadd_list(i)%id, mergeadd_list(i)%status
 200      FORMAT(A, 2(1X, I8))
          WRITE(LUN, 300) mergeadd_list(i)%xh(:)
 300      FORMAT(3(E23.16, 1X))
          WRITE(LUN, 300) mergeadd_list(i)%vh(:)
          ncomp = mergeadd_list(i)%ncomp
          DO j = 1, ncomp
               index = index + 1
               WRITE(LUN, 200) SUB, mergesub_list(index)%id, mergesub_list(index)%status
               WRITE(LUN, 300) mergesub_list(index)%xh(:)
               WRITE(LUN, 300) mergesub_list(index)%vh(:)
          END DO
     END DO
     symba_plP => symba_pld1P
     DO i = 1, nsppl
          swifter_plP => symba_plP%helio%swifter
          IF (swifter_plP%status /= MERGED) THEN
               WRITE(LUN, 200) SUB, swifter_plP%id, swifter_plP%status
               WRITE(LUN, 300) swifter_plP%xh(:)
               WRITE(LUN, 300) swifter_plP%vh(:)
          END IF
          symba_plP => symba_plP%nextP
     END DO
     symba_tpP => symba_tpd1P
     DO i = 1, nsptp
          swifter_tpP => symba_tpP%helio%swifter
          WRITE(LUN, 200) SUB, swifter_tpP%id, swifter_tpP%status
          WRITE(LUN, 300) swifter_tpP%xh(:)
          WRITE(LUN, 300) swifter_tpP%vh(:)
          symba_tpP => symba_tpP%nextP
     END DO
     IF (lbig_discard) THEN
          nplm = 0
          swifter_plP => symba_pl1P%helio%swifter
          DO i = 1, npl
               IF (swifter_plP%mass < mtiny) EXIT
               nplm = nplm + 1
               swifter_plP => swifter_plP%nextP
          END DO
          IF (nplm > 1) THEN
               WRITE(LUN, 400) nplm
 400           FORMAT(I8)
               swifter_plP => symba_pl1P%helio%swifter
               DO i = 2, nplm
                    swifter_plP => swifter_plP%nextP
                    WRITE(LUN, 500) swifter_plP%id, swifter_plP%mass, swifter_plP%radius
 500                FORMAT(I8, 2(1X, E23.16))
                    WRITE(LUN, 300) swifter_plP%xh(:)
                    WRITE(LUN, 300) swifter_plP%vh(:)
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
