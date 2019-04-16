!**********************************************************************************************************************************
!
!  Unit Name   : helio_discard
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Call discard routine to determine spilled test particles, then remove them from active list
!
!  Input
!    Arguments : t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                nsp            : number of spilled test particles
!                helio_pl1P     : pointer to head of helio planet structure linked-list
!                helio_tp1P     : pointer to head of active helio test particle structure linked-list
!                helio_tpd1P    : pointer to head of discard helio test particle structure linked-list
!                dt             : time step
!                rmin           : minimum allowed heliocentric radius for test particles
!                rmax           : maximum allowed heliocentric radius for test particles
!                rmaxu          : maximum allowed heliocentric radius for unbound test particles
!                qmin           : minimum allowed pericenter distance for test particles
!                qmin_coord     : coordinate frame for qmin
!                qmin_alo       : minimum semimajor axis for qmin
!                qmin_ahi       : maximum semimajor axis for qmin
!                lclose         : logical flag indicating whether to check for close planet-test particle encounters
!                lrhill_present : logical flag indicating whether Hill sphere radii for planets are present
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : ntp            : number of active test particles
!                nsp            : number of spilled test particles
!                helio_tp1P     : pointer to head of active helio test particle structure linked-list
!                helio_tpd1P    : pointer to head of discard helio test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_discard(t, npl, ntp, nsp, helio_pl1P, helio_tp1P, helio_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord,
!                                   qmin_alo, qmin_ahi, lclose, lrhill_present)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE helio_discard(t, npl, ntp, nsp, helio_pl1P, helio_tp1P, helio_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord,          &
     qmin_alo, qmin_ahi, lclose, lrhill_present)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_discard
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)    :: npl
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)    :: qmin_coord
     TYPE(helio_pl), POINTER     :: helio_pl1P
     TYPE(helio_tp), POINTER     :: helio_tp1P, helio_tpd1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P, swifter_tpspP
     TYPE(helio_tp), POINTER   :: helio_tpP, helio_tpspP

! Executable code
     swifter_pl1P => helio_pl1P%swifter
     swifter_tp1P => helio_tp1P%swifter
     CALL discard(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, qmin_coord, lclose,   &
          lrhill_present)
     helio_tpP => helio_tp1P
     DO i = 1, ntp
          helio_tpspP => helio_tpP
          helio_tpP => helio_tpP%nextP
          swifter_tpspP => helio_tpspP%swifter
          IF (swifter_tpspP%status /= ACTIVE) CALL helio_discard_spill(ntp, nsp, helio_tp1P, helio_tpd1P, helio_tpspP)
     END DO

     RETURN

END SUBROUTINE helio_discard
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
