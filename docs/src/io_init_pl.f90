!**********************************************************************************************************************************
!
!  Unit Name   : io_init_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in planet data
!
!  Input
!    Arguments : inplfile       : name of input file for planets
!                in_type        : format of input data file
!                lclose         : logical flag indicating whether planets' physical radii are present in input file
!                lrhill_present : logical flag indicating whether planets' Hill's sphere radii are present in input file
!                npl            : number of planets
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : id             : planet identifier               (all planets)
!                mass           : mass                            (all planets)
!                rhill          : Hill's sphere radius (optional) (if present, all planets except the Sun)
!                radius         : physical radius (optional)      (if present, all planets except the Sun)
!                xh             : heliocentric position           (all planets)
!                vh             : heliocentric velocity           (all planets)
!
!  Output
!    Arguments : swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_init_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)

! Modules
     USE swiftest
     USE module_symba
     USE module_helio
     USE module_swiftest
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_init_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)         :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)         :: npl
     CHARACTER(*), INTENT(IN)         :: inplfile, in_type
     TYPE(symba_pl), INTENT(INOUT)    :: symba_plA

! Internals
     INTEGER(I4B), PARAMETER          :: LUN = 7
     INTEGER(I4B)                     :: i, iu, ierr, inpl
    
! Executable code
     IF (in_type == "ASCII") THEN
          CALL io_open(LUN, inplfile, "OLD", "FORMATTED", ierr)
          READ(LUN, *) inpl
          READ(LUN, *) symba_plA%helio%swiftest%name(1), symba_plA%helio%swiftest%mass(1)
          symba_plA%helio%swiftest%rhill(1) = 0.0_DP
          symba_plA%helio%swiftest%radius(1) = 0.0_DP
          READ(LUN, *) symba_plA%helio%swiftest%xh(:,1)
          READ(LUN, *) symba_plA%helio%swiftest%vh(:,1)
          DO i = 1, NDIM
               IF ((symba_plA%helio%swiftest%xh(i,1) /= 0.0_DP) .OR. (symba_plA%helio%swiftest%vh(i,1) /= 0.0_DP)) THEN
                    WRITE(*, *) "SWIFTEST Error:"
                    WRITE(*, *) " Input MUST be in heliocentric coordinates."
                    WRITE(*, *) " Position/velocity components of Body 1 are"
                    WRITE(*, *) symba_plA%helio%swiftest%xh(:,1)
                    WRITE(*, *) symba_plA%helio%swiftest%vh(:,1)
                    CALL util_exit(FAILURE)
               END IF
          END DO
          symba_plA%helio%swiftest%status(1) = ACTIVE
          DO i = 2, npl
               IF (lrhill_present) THEN
                    READ(LUN, *) symba_plA%helio%swiftest%name(i), symba_plA%helio%swiftest%mass(i), &
                    symba_plA%helio%swiftest%rhill(i)
               ELSE
                    READ(LUN, *) symba_plA%helio%swiftest%name(i), symba_plA%helio%swiftest%mass(i)
                    symba_plA%helio%swiftest%rhill(i) = 0.0_DP
               END IF
               IF (lclose) THEN
                    READ(LUN, *) symba_plA%helio%swiftest%radius(i)
               ELSE
                    symba_plA%helio%swiftest%radius(i) = 0.0_DP
               END IF
               READ(LUN, *) symba_plA%helio%swiftest%xh(:,i)
               READ(LUN, *) symba_plA%helio%swiftest%vh(:,i)
               symba_plA%helio%swiftest%status(i) = ACTIVE
          END DO
          CLOSE(UNIT = LUN)
     ELSE
          CALL io_open_fxdr(inplfile, "R", .TRUE., iu, ierr)
          ierr = ixdrint(iu, inpl)
          ierr = ixdrint(iu, symba_plA%helio%swiftest%name(1))
          ierr = ixdrdouble(iu, symba_plA%helio%swiftest%mass(1))
          symba_plA%helio%swiftest%rhill(1) = 0.0_DP
          symba_plA%helio%swiftest%radius(1) = 0.0_DP
          ierr = ixdrdmat(iu, NDIM, symba_plA%helio%swiftest%xh(:,1))
          ierr = ixdrdmat(iu, NDIM, symba_plA%helio%swiftest%vh(:,1))
          DO i = 1, NDIM
               IF ((symba_plA%helio%swiftest%xh(i,1) /= 0.0_DP) .OR. &
                    (symba_plA%helio%swiftest%vh(i,1) /= 0.0_DP)) THEN
                    WRITE(*, *) "SWIFTEST Error:"
                    WRITE(*, *) " Input MUST be in heliocentric coordinates."
                    WRITE(*, *) " Position/velocity components of Body 1 are"
                    WRITE(*, *) symba_plA%helio%swiftest%xh(:,1)
                    WRITE(*, *) symba_plA%helio%swiftest%vh(:,1)
                    CALL util_exit(FAILURE)
               END IF
          END DO
          symba_plA%helio%swiftest%status(1) = ACTIVE
          DO i = 2, npl
               ierr = ixdrint(iu, symba_plA%helio%swiftest%name(i))
               ierr = ixdrdouble(iu, symba_plA%helio%swiftest%mass(i))
               IF (lrhill_present) THEN
                    ierr = ixdrdouble(iu, symba_plA%helio%swiftest%rhill(i))
               ELSE
                    symba_plA%helio%swiftest%rhill(i) = 0.0_DP
               END IF
               IF (lclose) THEN
                    ierr = ixdrdouble(iu, symba_plA%helio%swiftest%radius(i))
               ELSE
                    symba_plA%helio%swiftest%radius(i) = 0.0_DP
               END IF
               ierr = ixdrdmat(iu, NDIM, symba_plA%helio%swiftest%xh(:,i))
               ierr = ixdrdmat(iu, NDIM, symba_plA%helio%swiftest%vh(:,i))
               symba_plA%helio%swiftest%status(i) = ACTIVE
          END DO
          ierr = ixdrclose(iu)
     END IF

     RETURN

END SUBROUTINE io_init_pl
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann (Checked by Jennifer Pouplin & Carlisle Wishard)
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
