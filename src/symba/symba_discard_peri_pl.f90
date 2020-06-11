!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_peri_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their pericenter distances
!
!  Input
!    Arguments : t          : time
!                npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!                msys       : total system mass
!                qmin       : minimum allowed pericenter distance
!                qmin_alo   : minimum semimajor axis for qmin
!                qmin_ahi   : maximum semimajor axis for qmin
!                qmin_coord : coordinate frame for qmin
!                ldiscards  : logical flag indicating whether any planets are discarded
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!                ldiscards  : logical flag indicating whether any planets are discarded
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL symba_discard_peri_pl(t, npl, symba_pl1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_peri.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)

! Modules
     USE swiftest
     USE helio
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_peri_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(INOUT)   :: ldiscards
     INTEGER(I4B), INTENT(IN)      :: npl
     REAL(DP), INTENT(IN)          :: t, msys, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)      :: qmin_coord
     TYPE(symba_pl), INTENT(INOUT) :: symba_plA

! Internals
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B)              :: i, j, ih
     REAL(DP)                  :: r2
     REAL(DP), DIMENSION(NDIM) :: dx


! Executable code
     IF (lfirst) THEN
          CALL symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
          lfirst = .FALSE.
     ELSE
          CALL symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
          DO i = 2, npl
               IF (symba_plA%helio%swiftest%status(i) == ACTIVE) THEN
                    IF ((symba_plA%isperi(i) == 0) .AND. (symba_plA%nplenc(i)== 0)) THEN
                         IF ((symba_plA%atp(i) >= qmin_alo) .AND. (symba_plA%atp(i) <= qmin_ahi) &
                          .AND. (symba_plA%peri(i) <= qmin)) THEN
                              ldiscards = .TRUE.
                              symba_plA%helio%swiftest%status(i) = DISCARDED_PERI
                              WRITE(*, *) "Particle ", symba_plA%helio%swiftest%name(i), &
                               " perihelion distance too small at t = ", t
                         END IF
                    END IF
               END IF
          END DO
     END IF

     RETURN

END SUBROUTINE symba_discard_peri_pl
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
