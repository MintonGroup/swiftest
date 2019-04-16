!**********************************************************************************************************************************
!
!  Unit Name   : whm_getacch_ah3_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross (third) term heliocentric accelerations of test particles
!
!  Input
!    Arguments : npl      : number of planets
!                ntp      : number of active test particles
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                whm_tp1P : pointer to head of active WHM test particle structure linked-list
!                xh       : heliocentric planet positions
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_tp1P : pointer to head of active WHM test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_getacch_ah3_tp(npl, ntp, whm_pl1P, whm_tp1P, xh)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_getacch_ah3_tp(npl, ntp, whm_pl1P, whm_tp1P, xh)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_getacch_ah3_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                   :: npl, ntp
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(whm_pl), POINTER                      :: whm_pl1P
     TYPE(whm_tp), POINTER                      :: whm_tp1P

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: rji2, irij3, fac
     REAL(DP), DIMENSION(NDIM) :: dx, acc, xht
     TYPE(whm_pl), POINTER     :: whm_plP
     TYPE(whm_tp), POINTER     :: whm_tpP

! Executable code
     whm_tpP => whm_tp1P
     DO i = 1, ntp
          acc(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          xht(:) = whm_tpP%swifter%xh(:)
          whm_plP => whm_pl1P
          DO j = 2, npl
               whm_plP => whm_plP%nextP
               dx(:) = xht(:) - xh(:, j)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               fac = whm_plP%swifter%mass*irij3
               acc(:) = acc(:) - fac*dx(:)
          END DO
          whm_tpP%ah(:) = acc(:)
          whm_tpP => whm_tpP%nextP
     END DO

     RETURN

END SUBROUTINE whm_getacch_ah3_tp
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
