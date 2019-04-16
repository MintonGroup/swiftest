!**********************************************************************************************************************************
!
!  Unit Name   : symba_kick
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of planets and active test particles within SyMBA recursion
!
!  Input
!    Arguments : irec         : input recursion level
!                nplplenc     : number of planet-planet encounters
!                npltpenc     : number of planet-test particle encounters
!                plplenc_list : array of planet-planet encounter structures
!                pltpenc_list : array of planet-test particle encounter structures
!                dt           : time step
!                sgn          : sign to be applied to acceleration
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : plplenc_list : array of planet-planet encounter structures
!                pltpenc_list : array of planet-test particle encounter structures
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_kick.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_kick
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                      :: irec, nplplenc, npltpenc
     REAL(DP), INTENT(IN)                          :: dt, sgn
     TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(:), INTENT(IN) :: pltpenc_list

! Internals
     INTEGER(I4B)              :: i, j, irm1, irecl
     REAL(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
     REAL(DP), DIMENSION(NDIM) :: dx
     TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(helio_pl), POINTER   :: helio_pliP, helio_pljP
     TYPE(helio_tp), POINTER   :: helio_tpP
     TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     irm1 = irec - 1
     IF (sgn < 0.0_DP) THEN
          irecl = irec - 1
     ELSE
          irecl = irec
     END IF
     DO i = 1, nplplenc
          helio_pliP => plplenc_list(i)%pl1P%helio
          helio_pljP => plplenc_list(i)%pl2P%helio
          helio_pliP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          helio_pljP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, npltpenc
          helio_tpP => pltpenc_list(i)%tpP%helio
          helio_tpP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, nplplenc
          IF (plplenc_list(i)%status == ACTIVE) THEN
               symba_pliP => plplenc_list(i)%pl1P
               symba_pljP => plplenc_list(i)%pl2P
               IF ((symba_pliP%levelg >= irm1) .AND. (symba_pljP%levelg >= irm1)) THEN
                    helio_pliP => symba_pliP%helio
                    helio_pljP => symba_pljP%helio
                    swifter_pliP => helio_pliP%swifter
                    swifter_pljP => helio_pljP%swifter
                    ri = ((swifter_pliP%rhill + swifter_pljP%rhill)**2)*(RHSCALE**2)*(RSHELL**(2*irecl))
                    rim1 = ri*(RSHELL**2)
                    dx(:) = swifter_pljP%xh(:) - swifter_pliP%xh(:)
                    r2 = DOT_PRODUCT(dx(:), dx(:))
                    IF (r2 < rim1) THEN
                         fac = 0.0_DP
                    ELSE IF (r2 < ri) THEN
                         ris = SQRT(ri)
                         r = SQRT(r2)
                         rr = (ris - r)/(ris*(1.0_DP - RSHELL))
                         fac = (r2**(-1.5_DP))*(1.0_DP - 3.0_DP*(rr**2) + 2.0_DP*(rr**3))
                    ELSE
                         ir3 = 1.0_DP/(r2*SQRT(r2))
                         fac = ir3
                    END IF
                    faci = fac*swifter_pliP%mass
                    facj = fac*swifter_pljP%mass
                    helio_pliP%ah(:) = helio_pliP%ah(:) + facj*dx(:)
                    helio_pljP%ah(:) = helio_pljP%ah(:) - faci*dx(:)
               END IF
          END IF
     END DO
     DO i = 1, npltpenc
          IF (pltpenc_list(i)%status == ACTIVE) THEN
               symba_pliP => pltpenc_list(i)%plP
               symba_tpP => pltpenc_list(i)%tpP
               IF ((symba_pliP%levelg >= irm1) .AND. (symba_tpP%levelg >= irm1)) THEN
                    helio_pliP => symba_pliP%helio
                    helio_tpP => symba_tpP%helio
                    swifter_pliP => helio_pliP%swifter
                    swifter_tpP => helio_tpP%swifter
                    ri = ((swifter_pliP%rhill)**2)*(RHSCALE**2)*(RSHELL**(2*irecl))
                    rim1 = ri*(RSHELL**2)
                    dx(:) = swifter_tpP%xh(:) - swifter_pliP%xh(:)
                    r2 = DOT_PRODUCT(dx(:), dx(:))
                    IF (r2 < rim1) THEN
                         fac = 0.0_DP
                    ELSE IF (r2 < ri) THEN
                         ris = SQRT(ri)
                         r = SQRT(r2)
                         rr = (ris - r)/(ris*(1.0_DP - RSHELL))
                         fac = (r2**(-1.5_DP))*(1.0_DP - 3.0_DP*(rr**2) + 2.0_DP*(rr**3))
                    ELSE
                         ir3 = 1.0_DP/(r2*SQRT(r2))
                         fac = ir3
                    END IF
                    faci = fac*swifter_pliP%mass
                    helio_tpP%ah(:) = helio_tpP%ah(:) - faci*dx(:)
               END IF
          END IF
     END DO
     DO i = 1, nplplenc
          helio_pliP => plplenc_list(i)%pl1P%helio
          helio_pljP => plplenc_list(i)%pl2P%helio
          swifter_pliP => helio_pliP%swifter
          swifter_pljP => helio_pljP%swifter
          swifter_pliP%vb(:) = swifter_pliP%vb(:) + sgn*dt*helio_pliP%ah(:)
          swifter_pljP%vb(:) = swifter_pljP%vb(:) + sgn*dt*helio_pljP%ah(:)
          helio_pliP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          helio_pljP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, npltpenc
          helio_tpP => pltpenc_list(i)%tpP%helio
          swifter_tpP => helio_tpP%swifter
          IF (swifter_tpP%status == ACTIVE) swifter_tpP%vb(:) = swifter_tpP%vb(:) + sgn*dt*helio_tpP%ah(:)
          helio_tpP%ah(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO

     RETURN

END SUBROUTINE symba_kick
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
