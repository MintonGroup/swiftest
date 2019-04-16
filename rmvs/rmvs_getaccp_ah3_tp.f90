!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_getaccp_ah3_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross (third) term planetocentric accelerations of test particles closely encountering a planet
!
!  Input
!    Arguments : index     : inner substep number within current set
!                npl       : number of planets
!                nenc      : number of test particles encountering current planet
!                rmvs_pl1P : pointer to head of RMVS planet structure linked-list
!                rmvs_pleP : pointer to RMVS planet structure of planet being closely encountered
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_pleP : pointer to RMVS planet structure of planet being closely encountered
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_getaccp_ah3_tp(index, npl, nenc, rmvs_pl1P, rmvs_pleP)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_a3_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_getaccp_ah3_tp(index, npl, nenc, rmvs_pl1P, rmvs_pleP)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_getaccp_ah3_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: index, npl, nenc
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P, rmvs_pleP

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: rji2, irij3, fac
     REAL(DP), DIMENSION(NDIM) :: dx, acc, xpct
     TYPE(rmvs_pl), POINTER    :: rmvs_plP
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     rmvs_tpP => rmvs_pleP%tpenc1P
     DO i = 1, nenc
          acc(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          xpct(:) = rmvs_tpP%xpc(:)
          rmvs_plP => rmvs_pl1P
          DO j = 1, npl
               IF (.NOT. ASSOCIATED(rmvs_plP, rmvs_pleP)) THEN
                    dx(:) = xpct(:) - rmvs_plP%xpc(:, index)
                    rji2 = DOT_PRODUCT(dx(:), dx(:))
                    irij3 = 1.0_DP/(rji2*SQRT(rji2))
                    fac = rmvs_plP%whm%swifter%mass*irij3
                    acc(:) = acc(:) - fac*dx(:)
               END IF
               rmvs_plP => rmvs_plP%nextP
          END DO
          rmvs_tpP%apc(:) = acc(:)
          rmvs_tpP => rmvs_tpP%tpencP
     END DO

     RETURN

END SUBROUTINE rmvs_getaccp_ah3_tp
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
