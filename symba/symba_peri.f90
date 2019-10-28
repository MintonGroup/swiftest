!**********************************************************************************************************************************
!
!  Unit Name   : symba_peri
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Determine system pericenter passages for planets in SyMBA
!
!  Input
!    Arguments : lfirst     : logical flag indicating whether current invocation is the first
!                npl        : number of planets
!                symba_pl1P : pointer to head of SyMBA planet structure linked-list
!                msys       : total system mass
!                qmin_coord : coordinate frame for qmin
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P : pointer to head of SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_peri(lfirst, npl, symba_pl1P, msys, qmin_coord)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_mass_peri.f
!
!                If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
!                planet structures are up-to-date and are not recomputed
!
!**********************************************************************************************************************************
SUBROUTINE symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lfirst
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: msys
     CHARACTER(*), INTENT(IN) :: qmin_coord
     TYPE(symba_pl), INTENT(INOUT)  :: symba_plA

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: vdotr, e, mu, msun

! Executable code
     msun = symba_plA%helio%swifter%mass(1)
     IF (lfirst) THEN
          IF (qmin_coord == "HELIO") THEN
               DO i = 2, npl
                    IF (symba_plA%helio%swiftest%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%vh(:,i))
                         IF (vdotr > 0.0_DP) THEN
                              symba_plP%isperi(i) = 1
                         ELSE
                              symba_plP%isperi(i) = -1
                         END IF
                    END IF
               END DO
          ELSE
               DO i = 2, npl
                    IF (symba_plA%helio%swiftest%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(symba_plA%helio%swiftest%xb(:,i), symba_plA%helio%swiftest%vb(:,i))
                         IF (vdotr > 0.0_DP) THEN
                              symba_plP%isperi(i) = 1
                         ELSE
                              symba_plP%isperi(i) = -1
                         END IF
                    END IF
               END DO
          END IF
     ELSE
          IF (qmin_coord == "HELIO") THEN
               DO i = 2, npl
                    IF (symba_plA%helio%swiftest%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%vh(:,i))
                         IF (symba_plP%isperi(i) == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   symba_plP%isperi(i) = 0
                                   mu = msun + symba_plA%helio%swiftest%mass(i)
                                   CALL orbel_xv2aeq(symba_plA%helio%swiftest%xh(:,i), symba_plA%helio%swiftest%vh(:,i), mu, symba_plP%atp(i), e, symba_plP%peri(i))
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   symba_plP%isperi(i) = 1
                              ELSE
                                   symba_plP%isperi(i) = -1
                              END IF
                         END IF
                    END IF
               END DO
          ELSE
               DO i = 2, npl
                    IF (symba_plA%helio%swiftest%status(i) == ACTIVE) THEN
                         vdotr = DOT_PRODUCT(symba_plA%helio%swiftest%xb(:,i), symba_plA%helio%swiftest%vb(:,i))
                         IF (symba_plP%isperi(i) == -1) THEN
                              IF (vdotr >= 0.0_DP) THEN
                                   symba_plP%isperi(i) = 0
                                   CALL orbel_xv2aeq(symba_plA%helio%swiftest%xb(:,i), symba_plA%helio%swiftest%vb(:.i), msys, symba_plP%atp(i), e, symba_plP%peri(i))
                              END IF
                         ELSE
                              IF (vdotr > 0.0_DP) THEN
                                   symba_plP%isperi(i) = 1
                              ELSE
                                   symba_plP%isperi(i) = -1
                              END IF
                         END IF
                    END IF
               END DO
          END IF
     END IF

     RETURN

END SUBROUTINE symba_peri
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
