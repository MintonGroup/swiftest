!**********************************************************************************************************************************
!
!  Unit Name   : bs_derivs
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Compute time derivatives of positions and velocities of planets and test particles for Bulirsch-Stoer method
!
!  Input
!    Arguments : lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_derivs(lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_derivs(lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_derivs
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
     TYPE(bs_pl), POINTER     :: bs_pl1P
     TYPE(bs_tp), POINTER     :: bs_tp1P

! Internals
     INTEGER(I4B)         :: i
     TYPE(bs_pl), POINTER :: bs_plP
     TYPE(bs_tp), POINTER :: bs_tpP

! Executable code
     CALL bs_getaccb(lextra_force, t, npl, nplmax, bs_pl1P, j2rp2, j4rp4)
     IF (ntp > 0) CALL bs_getaccb_tp(lextra_force, t, npl, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
     bs_plP => bs_pl1P
     DO i = 1, npl
          bs_plP%dydx(1:NDIM) = bs_plP%y(NDIM+1:NDIM2)
          bs_plP%dydx(NDIM+1:NDIM2) = bs_plP%ab(:)
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%dydx(1:NDIM) = bs_tpP%y(NDIM+1:NDIM2)
          bs_tpP%dydx(NDIM+1:NDIM2) = bs_tpP%ab(:)
          bs_tpP => bs_tpP%nextP
     END DO

     RETURN

END SUBROUTINE bs_derivs
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
