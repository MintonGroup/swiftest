!**********************************************************************************************************************************
!
!  Unit Name   : discard_peri
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : discard
!  Language    : Fortran 90/95
!
!  Description : Check to see if a test particle should be discarded because its perihelion distance becomes too small
!
!  Input
!    Arguments : t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!                msys           : total system mass
!                qmin           : minimum pericenter distance for test particle
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                qmin_coord     : coordinate frame to use for qmin
!                lrhill_present : logical flag indicating whether Hill sphere radii for planets are present
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P   : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL discard_peri(t, npl, ntp, swifter_pl1P, swifter_tp1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,
!                                  lrhill_present)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_peri.f
!
!**********************************************************************************************************************************
SUBROUTINE discard_peri(t, npl, ntp, swifter_pl1P, swifter_tp1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, lrhill_present)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => discard_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lrhill_present
     INTEGER(I4B), INTENT(IN)  :: npl, ntp
     REAL(DP), INTENT(IN)      :: t, msys, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)  :: qmin_coord
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P

! Internals
     LOGICAL(LGT), SAVE        :: lfirst = .TRUE.
     INTEGER(I4B)              :: i, j, ih
     REAL(DP)                  :: r2
     REAL(DP), DIMENSION(NDIM) :: dx
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     IF (lfirst) THEN
          IF (.NOT. lrhill_present) CALL util_hills(npl, swifter_pl1P)
          CALL util_peri(lfirst, ntp, swifter_tp1P, swifter_pl1P%mass, msys, qmin_coord)
          lfirst = .FALSE.
     ELSE
          CALL util_peri(lfirst, ntp, swifter_tp1P, swifter_pl1P%mass, msys, qmin_coord)
          swifter_tpP => swifter_tp1P
          DO i = 1, ntp
               IF (swifter_tpP%status == ACTIVE) THEN
                    IF (swifter_tpP%isperi == 0) THEN
                         ih = 1
                         swifter_plP => swifter_pl1P
                         DO j = 2, npl
                              swifter_plP => swifter_plP%nextP
                              dx(:) = swifter_tpP%xh(:) - swifter_plP%xh(:)
                              r2 = DOT_PRODUCT(dx(:), dx(:))
                              IF (r2 <= swifter_plP%rhill*swifter_plP%rhill) ih = 0
                         END DO
                         IF (ih == 1) THEN
                              IF ((swifter_tpP%atp >= qmin_alo) .AND. (swifter_tpP%atp <= qmin_ahi) .AND.                         &
                                  (swifter_tpP%peri <= qmin)) THEN
                                   swifter_tpP%status = DISCARDED_PERI
                                   WRITE(*, *) "Particle ", swifter_tpP%id, " perihelion distance too small at t = ", t
                              END IF
                         END IF
                    END IF
               END IF
               swifter_tpP => swifter_tpP%nextP
          END DO
     END IF

     RETURN

END SUBROUTINE discard_peri
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
