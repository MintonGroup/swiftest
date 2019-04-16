!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch_int_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross term heliocentric accelerations of test particles
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                helio_tp1P   : pointer to head of active helio test particle structure linked-list
!                xh           : heliocentric planet positions
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_tp1P   : pointer to head of active helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_getacch_int_tp(npl, ntp, swifter_pl1P, helio_tp1P, xh)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah3_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_getacch_int_tp(npl, ntp, swifter_pl1P, helio_tp1P, xh)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch_int_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                   :: npl, ntp
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(swifter_pl), POINTER                  :: swifter_pl1P
     TYPE(helio_tp), POINTER                    :: helio_tp1P

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: r2, fac
     REAL(DP), DIMENSION(NDIM) :: dx
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(helio_tp), POINTER   :: helio_tpP

! Executable code
     helio_tpP => helio_tp1P
     DO i = 1, ntp
          swifter_tpP => helio_tpP%swifter
          IF (swifter_tpP%status == ACTIVE) THEN
               swifter_plP => swifter_pl1P
               DO j = 2, npl
                    swifter_plP => swifter_plP%nextP
                    dx(:) = swifter_tpP%xh(:) - swifter_plP%xh(:)
                    r2 = DOT_PRODUCT(dx(:), dx(:))
                    fac = swifter_plP%mass/(r2*SQRT(r2))
                    helio_tpP%ahi(:) = helio_tpP%ahi(:) - fac*dx(:)
               END DO
          END IF
          helio_tpP => helio_tpP%nextP
     END DO

     RETURN

END SUBROUTINE helio_getacch_int_tp
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
