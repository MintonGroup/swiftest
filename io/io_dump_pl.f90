!**********************************************************************************************************************************
!
!  Unit Name   : io_dump_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Dump planet data to file
!
!  Input
!    Arguments : npl            : number of planets
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                lclose         : logical flag indicating whether to check for planet-test particle encounters
!                lrhill_present : logical flag indicating whether Hill's sphere radii are present in planet data
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : to dump file
!                npl            : number of planets
!                id             : planet identifier     (from planet structure, for each planet)
!                mass           : mass                  (from planet structure, for each planet)
!                rhill          : Hill's sphere radius  (from planet structure, for each planet except the Sun, if lrhill_present)
!                radius         : planet radius         (from planet structure, for each planet except the Sun, if lclose)
!                xh             : heliocentric position (from planet structure, for each planet)
!                vh             : heliocentric velocity (from planet structure, for each planet)
!
!  Invocation  : CALL io_dump_pl(npl, swifter_pl1P, lclose, lrhill_present)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_dump_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE io_dump_pl(npl, swiftest_plA, lclose, lrhill_present)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_dump_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)         :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)         :: npl
     TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA

! Internals
     INTEGER(I4B)                     :: i, iu, ierr
     INTEGER(I4B), SAVE               :: idx = 1

! Executable code
     CALL io_open_fxdr(DUMP_PL_FILE(idx), "W", .TRUE., iu, ierr)
     IF (ierr /= 0) THEN
          WRITE(*, *) "SWIFTEST Error:"
          WRITE(*, *) "   Unable to open binary dump file ", TRIM(DUMP_PL_FILE(idx))
          CALL util_exit(FAILURE)
     END IF
     ierr = ixdrint(iu, npl)
     ierr = ixdrint(iu, swiftest_plA%id(1))
     ierr = ixdrdouble(iu, swiftest_plA%mass(1))
     ierr = ixdrdmat(iu, NDIM, swiftest_plA%xh(:,1))
     ierr = ixdrdmat(iu, NDIM, swiftest_plA%vh(:,1))
     DO i = 2, npl
          ierr = ixdrint(iu, swiftest_plA%id(i))
          ierr = ixdrdouble(iu, swiftest_plA%mass(i))
          IF (lrhill_present) ierr = ixdrdouble(iu, swifter_plA%rhill(i))
          IF (lclose) ierr = ixdrdouble(iu, swifter_plA%radius(i))
          ierr = ixdrdmat(iu, NDIM, swifter_plA%xh(:,i))
          ierr = ixdrdmat(iu, NDIM, swifter_plA%vh(:,i))
     END DO
     ierr = ixdrclose(iu)
     idx = idx + 1
     IF (idx > 2) idx = 1

     RETURN

END SUBROUTINE io_dump_pl
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
