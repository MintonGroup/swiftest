!**********************************************************************************************************************************
!
!  Unit Name   : symba_kick
!  Unit Type   : subroutine
!  Project     : Swiftest
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
SUBROUTINE symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_plA, symba_tpA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_kick
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: irec, nplplenc, npltpenc
     REAL(DP), INTENT(IN)                             :: dt, sgn
     TYPE(symba_plplenc), INTENT(IN)                  :: plplenc_list
     TYPE(symba_pltpenc), INTENT(IN)                  :: pltpenc_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA

! Internals
     INTEGER(I4B)              :: i, irm1, irecl, index_i,index_j,index_tp,index_pl
     REAL(DP)                  :: r, rr, ri, ris, rim1, r2, ir3, fac, faci, facj
     REAL(DP), DIMENSION(NDIM) :: dx

! Executable code
     irm1 = irec - 1
     IF (sgn < 0.0_DP) THEN
          irecl = irec - 1
     ELSE
          irecl = irec
     END IF
     DO i = 1, nplplenc
          index_i  = plplenc_list%index1(i)
          index_j  = plplenc_list%index2(i)
          symba_plA%helio%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          symba_plA%helio%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, npltpenc
          index_tp  = pltpenc_list%indextp(i)
          symba_tpA%helio%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, nplplenc
          IF (plplenc_list%status(i) == ACTIVE) THEN
               index_i  = plplenc_list%index1(i) 
               index_j  = plplenc_list%index2(i) 
               IF ((symba_plA%levelg(index_i) >= irm1) .AND. (symba_plA%levelg(index_j) >= irm1)) THEN
                    ri = ((symba_plA%helio%swiftest%rhill(index_i) &
                         + symba_plA%helio%swiftest%rhill(index_j))**2)*(RHSCALE**2)*(RSHELL**(2*irecl))
                    rim1 = ri*(RSHELL**2)
                    dx(:) = symba_plA%helio%swiftest%xh(:,index_j) - symba_plA%helio%swiftest%xh(:,index_i)
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
                    faci = fac*symba_plA%helio%swiftest%mass(index_i)
                    facj = fac*symba_plA%helio%swiftest%mass(index_j)
                    symba_plA%helio%ah(:,index_i) = symba_plA%helio%ah(:,index_i) + facj*dx(:)
                    symba_plA%helio%ah(:,index_j) = symba_plA%helio%ah(:,index_j) - faci*dx(:)
               END IF
          END IF
     END DO
     DO i = 1, npltpenc
          IF (pltpenc_list%status(i) == ACTIVE) THEN
               index_pl  = pltpenc_list%indexpl(i) 
               index_tp  = pltpenc_list%indextp(i) 
               IF ((symba_plA%levelg(index_pl) >= irm1) .AND. (symba_tpA%levelg(index_tp) >= irm1)) THEN
                    ri = ((symba_plA%helio%swiftest%rhill(index_pl))**2)*(RHSCALE**2)*(RSHELL**(2*irecl))
                    rim1 = ri*(RSHELL**2)
                    dx(:) = symba_tpA%helio%swiftest%xh(:,index_tp) - symba_plA%helio%swiftest%xh(:,index_pl)
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
                    faci = fac*symba_plA%helio%swiftest%mass(index_pl)
                    symba_tpA%helio%ah(:,index_tp) = symba_tpA%helio%ah(:,index_tp) - faci*dx(:)
               END IF
          END IF
     END DO
     DO i = 1, nplplenc
          index_i  = plplenc_list%index1(i) 
          index_j  = plplenc_list%index2(i) 
          symba_plA%helio%swiftest%vb(:,index_i) = symba_plA%helio%swiftest%vb(:,index_i) + sgn*dt*symba_plA%helio%ah(:,index_i)
          symba_plA%helio%swiftest%vb(:,index_j) = symba_plA%helio%swiftest%vb(:,index_j) + sgn*dt*symba_plA%helio%ah(:,index_j)
          symba_plA%helio%ah(:,index_i) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          symba_plA%helio%ah(:,index_j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO
     DO i = 1, npltpenc
          index_tp  = pltpenc_list%indextp(i)
          IF (symba_tpA%helio%swiftest%status(index_tp) == ACTIVE)     &
          symba_tpA%helio%swiftest%vb(:,index_tp) = symba_tpA%helio%swiftest%vb(:,index_tp) + sgn*dt*symba_tpA%helio%ah(:,index_tp)
          symba_tpA%helio%ah(:,index_tp) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
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
