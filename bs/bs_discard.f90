!**********************************************************************************************************************************
!
!  Unit Name   : bs_discard
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Call discard routine to determine spilled test particles, then remove them from active list
!
!  Input
!    Arguments : t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                nsp            : number of spilled test particles
!                bs_pl1P        : pointer to head of BS planet structure linked-list
!                bs_tp1P        : pointer to head of active BS test particle structure linked-list
!                bs_tpd1P       : pointer to head of discard BS test particle structure linked-list
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
!                bs_tp1P        : pointer to head of active BS test particle structure linked-list
!                bs_tpd1P       : pointer to head of discard BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_discard(t, npl, ntp, nsp, bs_pl1P, bs_tp1P, bs_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,
!                                qmin_ahi, lclose, lrhill_present)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE bs_discard(t, npl, ntp, nsp, bs_pl1P, bs_tp1P, bs_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,  &
     lclose, lrhill_present)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_discard
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)    :: npl
     INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
     REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)    :: qmin_coord
     TYPE(bs_pl), POINTER        :: bs_pl1P
     TYPE(bs_tp), POINTER        :: bs_tp1P, bs_tpd1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_pl1P
     TYPE(swifter_tp), POINTER :: swifter_tp1P, swifter_tpspP
     TYPE(bs_tp), POINTER      :: bs_tpP, bs_tpspP

! Executable code
     swifter_pl1P => bs_pl1P%swifter
     swifter_tp1P => bs_tp1P%swifter
     CALL discard(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, qmin_coord, lclose,   &
          lrhill_present)
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpspP => bs_tpP
          bs_tpP => bs_tpP%nextP
          swifter_tpspP => bs_tpspP%swifter
          IF (swifter_tpspP%status /= ACTIVE) CALL bs_discard_spill(ntp, nsp, bs_tp1P, bs_tpd1P, bs_tpspP)
     END DO

     RETURN

END SUBROUTINE bs_discard
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
