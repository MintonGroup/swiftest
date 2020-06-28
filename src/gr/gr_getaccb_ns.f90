!**********************************************************************************************************************************
!
!  Unit Name   : gr_getaccb_ns
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : ra15
!  Language    : Fortran 90/95
!
!  Description : Add relativistic correction acceleration for non-symplectic
!  integrators
!
!  Input
!    Arguments : t         : time
!                npl       : number of planets
!                swifter_pl1P : pointer to head of RA15 planet structure linked-list
!    Terminal  : 
!    File      : 
!
!  Output
!    Arguments : agr       : Barycentric acceleration accelerations of planets due to GR
!    Terminal  : 
!    File      : 
!
!  Invocation  : CALL gr_getaccb_ns(npl, swiftger_pl1P, c2)
!
!  Notes       : Based on Quinn, Tremaine, & Duncan 1991
!
!**********************************************************************************************************************************
SUBROUTINE gr_getaccb_ns(npl, swifter_pl1P, agr, c2)

! Modules
     USE module_parameters
     USE module_ra15
     USE module_interfaces, EXCEPT_THIS_ONE => gr_getaccb_ns
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                    :: npl
     TYPE(swifter_pl), POINTER                   :: swifter_pl1P
     REAL(DP), INTENT(IN)                        :: c2
     REAL(DP), DIMENSION(npl, NDIM), INTENT(OUT) :: agr

! Internals
     REAL(DP), DIMENSION(NDIM) :: xh, vh
     TYPE(swifter_pl), POINTER :: swifter_plP
     REAL(DP)                  :: mu, msun, rmag, rdotv, vmag2
     INTEGER(I4B)              :: i

! Executable code
     swifter_plP => swifter_pl1P
     msun = swifter_plP%mass
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          mu = msun + swifter_plP%mass
          xh(:) = swifter_plP%xb(:) - swifter_pl1P%xb(:)
          vh(:) = swifter_plP%vb(:) - swifter_pl1P%vb(:)
          rmag = SQRT(DOT_PRODUCT(xh, xh))
          vmag2 = DOT_PRODUCT(vh, vh)
          rdotv = DOT_PRODUCT(xh, vh)
          agr(:,i) =  mu * c2 / rmag**3 * ((4 * mu / rmag - vmag2) * xh(:) + 4 * rdotv * vh(:))
     END DO

     agr(:, 1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          agr(:, 1) = agr(:, 1) - swifter_plP%mass*agr(:, i)/msun
     END DO

   

     RETURN

END SUBROUTINE gr_getaccb_ns
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
