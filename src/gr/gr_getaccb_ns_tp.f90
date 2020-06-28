!**********************************************************************************************************************************
!
!  Unit Name   : gr_getaccb_ns_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : swifter
!  Language    : Fortran 90/95
!
!  Description : Add relativistic correction acceleration to test particles for
!  non-symplectic integrators
!
!  Input
!    Arguments : t         : time
!                ntp       : number of tpanets
!                swifter_pl1P : pointer to head of RA15 panet structure linked-list
!                swifter_tp1P : pointer to head of RA15 tpanet structure linked-list
!                c2        : inverse speed of light squared
!    Terminal  : TBS as needed by user
!    File      : TBS as needed by user
!
!  Output
!    Arguments : agr       : barycentric acceleration
!    Terminal  : 
!    File      : 
!
!  Invocation  : CALL gr_getaccb_ns_tp(ntp, swifter_pl1P, swifter_tp1P, agr, c2)
!
!  Notes       : Based on Quinn, Tremaine, & Duncan 1991
!
!**********************************************************************************************************************************
SUBROUTINE gr_getaccb_ns_tp(ntp, swifter_pl1P, swifter_tp1P, agr, c2)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => gr_getaccb_ns_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     TYPE(swifter_pl), POINTER   :: swifter_pl1P
     TYPE(swifter_tp), POINTER   :: swifter_tp1P
     REAL(DP), INTENT(IN)     :: c2
     REAL(DP), DIMENSION(ntp, NDIM), INTENT(OUT) :: agr

! Internals
     REAL(DP), DIMENSION(NDIM) :: xh, vh
     TYPE(swifter_tp), POINTER :: swifter_tpP
     REAL(DP)                  :: mu, rmag, rv, vmag2
     INTEGER(I4B)              :: i

! Executable code
     mu = swifter_pl1P%mass
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          swifter_tpP => swifter_tpP%nextP
          xh(:) = swifter_tpP%xb(:) - swifter_pl1P%xb(:)
          vh(:) = swifter_tpP%vb(:) - swifter_pl1P%vb(:)
          rmag = SQRT(DOT_PRODUCT(xh, xh))
          vmag2 = DOT_PRODUCT(vh, vh)
          rv = DOT_PRODUCT(xh, vh)
          agr(:, i) = mu * c2 / rmag**3 *((4*mu /rmag - vmag2)*xh(:) + 4 * rv*vh(:))
     END DO

   

     RETURN

END SUBROUTINE gr_getaccb_ns_tp
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
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
