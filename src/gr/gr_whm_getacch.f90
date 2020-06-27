!**********************************************************************************************************************************
!
!  Unit Name   : gr_whm_getacch
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Compute relativisitic accelerations of planets
!
!  Input
!    Arguments : npl      : number of planets
!                whm_pl1P : pointer to head of WHM planet structure linked-list
!                ir3j     : inverse cubed Jacobi radii of planets
!                c2       : inverse speed of light squared
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : whm_pl1P : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL gr_whm_getacch(npl, whm_pl1P, c2)
!
!  Notes       : Based on Saha & Tremaine (1994) eq. 28
!
!**********************************************************************************************************************************
SUBROUTINE gr_whm_getacch(npl, whm_pl1P, c2)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => gr_whm_getacch
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)           :: npl
     TYPE(whm_pl), POINTER              :: whm_pl1P
     REAL(DP), INTENT(IN)               :: c2

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: etajm1, etaj, eta, beta, mu, rjmag4, msun
     TYPE(whm_pl), POINTER     :: whm_plP
     REAL(DP), DIMENSION(NDIM) :: suma
     REAL(DP), DIMENSION(npl, NDIM) :: aj

! Executable code
     whm_plP => whm_pl1P
     msun = whm_pl1P%swifter%mass
     etajm1 = msun
     DO i = 2, npl
          whm_plP => whm_plP%nextP
          etaj = etajm1 + whm_pl1P%swifter%mass
          mu = msun * etaj/etajm1
          rjmag4 = (DOT_PRODUCT(whm_plP%xj(:), whm_plP%xj(:)))**2
          beta = - mu**2*c2 
          aj(:, i) = beta * whm_plP%xj(:) / rjmag4
          IF (i == 2) THEN
               suma(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
               whm_plP%ah(:) = whm_plP%ah(:) + aj(:, 2)
          ELSE
               suma(:) = suma(:) + whm_plP%swifter%mass*aj(:, i)/etaj
               whm_plP%ah(:) = whm_plP%ah(:) + aj(:, i) + suma(:)
          END IF
          etajm1 = etaj
     END DO
     RETURN

END SUBROUTINE gr_whm_getacch
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
