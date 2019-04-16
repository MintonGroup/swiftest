!**********************************************************************************************************************************
!
!  Unit Name   : obl_pot
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : obl
!  Language    : Fortran 90/95
!
!  Description : Compute the contribution to the total gravitational potential due solely to the oblateness of the central body
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                xh           : heliocentric positions of the planets
!                irh          : inverse heliocentric radii of the planets
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : oblpot       : total gravitational potential due to J2 and J4 oblateness terms for the central body
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL obl_pot(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, oblpot)
!
!  Notes       : Adapted from Martin Duncan's Swift routine obl_pot.f
!
!                Returned value does not include monopole term or terms higher than J4
!
!                Reference: MacMillan, W. D. 1958. The Theory of the Potential, (Dover Publications), 363.
!
!**********************************************************************************************************************************
SUBROUTINE obl_pot(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, oblpot)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => obl_pot
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                   :: npl
     REAL(DP), INTENT(IN)                       :: j2rp2, j4rp4
     REAL(DP), INTENT(OUT)                      :: oblpot
     REAL(DP), DIMENSION(npl), INTENT(IN)       :: irh
     REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
     TYPE(swifter_pl), POINTER                  :: swifter_pl1P

! Internals
     INTEGER(I4B)              :: i
     REAL(DP)                  :: rinv2, t0, t1, t2, t3, p2, p4, mu
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     oblpot = 0.0_DP
     mu = swifter_pl1P%mass
     swifter_plP => swifter_pl1P
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          rinv2 = irh(i)**2
          t0 = mu*swifter_plP%mass*rinv2*irh(i)
          t1 = j2rp2
          t2 = xh(3, i)*xh(3, i)*rinv2
          t3 = j4rp4*rinv2
          p2 = 0.5_DP*(3.0_DP*t2 - 1.0_DP)
          p4 = 0.125_DP*((35.0_DP*t2 - 30.0_DP)*t2 + 3.0_DP)
          oblpot = oblpot + t0*(t1*p2 + t3*p4)
     END DO

     RETURN

END SUBROUTINE obl_pot
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
