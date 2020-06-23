!**********************************************************************************************************************************
!
!  Unit Name   : gr_pseudovel2vel
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Converts the relativistic pseudovelocity back into a
!  heliocentric velocity
!
!  Input
!    Arguments : xh, pvh  : heliocentric position and relativistic pseudovelocity velocity
!                c2       : inverse speed of light squared
!                mu       : Gravitational constant G*(Msun+Mp)
!    Terminal  : none
!
!  Output
!    Arguments : vh      : heliocentric velocity
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL gr_pseudovel2vel(xh,pvh,mu,c2,vh)
!
!  Notes       : Based on Saha & Tremaine (1994) Eq. 32
!
!
!**********************************************************************************************************************************
SUBROUTINE gr_pseudovel2vel(xh, pvh, mu, c2, vh)

! Modules
     USE module_parameters
     USE module_interfaces, EXCEPT_THIS_ONE => gr_pseudovel2vel
     IMPLICIT NONE

! Arguments
     REAL(DP), DIMENSION(NDIM), INTENT(IN)  :: xh,pvh
     REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: vh
     REAL(DP), INTENT(IN)                   :: mu, c2

! Internals
    REAL(DP) :: vmag2, rmag, grterm

! Executable code
     vmag2 = DOT_PRODUCT(pvh, pvh)
     rmag  = SQRT(DOT_PRODUCT(xh, xh))
     grterm = 1.0_DP - c2*(0.5_DP*vmag2 + 3*mu/rmag)
     vh(:) = pvh(:) * grterm

     RETURN

END SUBROUTINE gr_pseudovel2vel
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
