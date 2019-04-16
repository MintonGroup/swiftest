!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_chk
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Check to see if there are any encounters between planets and test particles
!
!  Input
!    Arguments : npl        : number of planets
!                ntp        : number of active test particles
!                rmvs_pl1P  : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P  : pointer to head of active RMVS test particle structure linked-list
!                xh         : planet positions at beginning of time step
!                vh         : planet velocities at beginning of time step
!                dt         : time step
!                rts        : fraction of Hill's sphere radius to use as radius of encounter region
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_pl1P  : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P  : pointer to head of active RMVS test particle structure linked-list
!                lencounter : logical flag indicating whether there are any encounters at all
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_chk(npl, ntp, rmvs_pl1P, rmvs_tp1P, xh, vh, dt, rts, lencounter)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_chk.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_chk(npl, ntp, rmvs_pl1P, rmvs_tp1P, xh, vh, dt, rts, lencounter)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_chk
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(OUT)                  :: lencounter
     INTEGER(I4B), INTENT(IN)                   :: npl, ntp
     REAL(DP), INTENT(IN)                       :: dt, rts
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh, vh
     TYPE(rmvs_pl), POINTER                     :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER                     :: rmvs_tp1P

! Internals
     INTEGER(I4B)              :: i, j, k, iflag, nenc
     REAL(DP)                  :: r2crit
     REAL(DP), DIMENSION(NDIM) :: xht, vht, xr, vr
     TYPE(rmvs_pl), POINTER    :: rmvs_plP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP, rmvs_tpencP

! Executable code
     lencounter = .FALSE.
     rmvs_plP => rmvs_pl1P
     DO i = 1, npl
          rmvs_plP%nenc = 0
          NULLIFY(rmvs_plP%tpenc1P)
          rmvs_plP => rmvs_plP%nextP
     END DO
     rmvs_tpP => rmvs_tp1P
     DO i = 1, ntp
          IF (rmvs_tpP%whm%swifter%status == ACTIVE) THEN
               NULLIFY(rmvs_tpP%plencP)
               NULLIFY(rmvs_tpP%tpencP)
               iflag = 0
               xht(:) = rmvs_tpP%whm%swifter%xh(:)
               vht(:) = rmvs_tpP%whm%swifter%vh(:)
               rmvs_plP => rmvs_pl1P
               DO j = 2, npl
                    rmvs_plP => rmvs_plP%nextP
                    r2crit = (rts*rmvs_plP%whm%swifter%rhill)**2
                    xr(:) = xht(:) - xh(:, j)
                    vr(:) = vht(:) - vh(:, j)
                    CALL rmvs_chk_ind(xr(:), vr(:), dt, r2crit, iflag)
                    IF (iflag /= 0) THEN
                         lencounter = .TRUE.
                         rmvs_plP%nenc = rmvs_plP%nenc + 1
                         nenc = rmvs_plP%nenc
                         IF (nenc == 1) THEN
                              rmvs_plP%tpenc1P => rmvs_tpP
                         ELSE
                              rmvs_tpencP => rmvs_plP%tpenc1P
                              DO k = 2, nenc - 1
                                   rmvs_tpencP => rmvs_tpencP%tpencP
                              END DO
                              rmvs_tpencP%tpencP => rmvs_tpP
                         END IF
                         rmvs_tpP%plencP => rmvs_plP
                         EXIT
                    END IF
               END DO
          END IF
          rmvs_tpP => rmvs_tpP%nextP
     END DO

     RETURN

END SUBROUTINE rmvs_chk
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
