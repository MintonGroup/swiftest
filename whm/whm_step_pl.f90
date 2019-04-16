!**********************************************************************************************************************************
!
!  Unit Name   : whm_step_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : whm
!  Language    : Fortran 90/95
!
!  Description : Step planets ahead using kick-drift-kick algorithm
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                whm_pl1P     : pointer to head of WHM planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL whm_step_pl(lfirst, lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine step_kdk_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE whm_step_pl(lfirst, lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4, dt)

! Modules
     USE module_parameters
     USE module_whm
     USE module_interfaces, EXCEPT_THIS_ONE => whm_step_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     TYPE(whm_pl), POINTER       :: whm_pl1P

! Internals
     REAL(DP) :: dth

! Executable code
     dth = 0.5_DP*dt
     IF (lfirst) THEN
          CALL coord_h2j(npl, whm_pl1P)
          CALL whm_getacch(lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4)
          lfirst = .FALSE.
     END IF
     CALL whm_kickvh(npl, whm_pl1P, dth)
     CALL coord_vh2vj(npl, whm_pl1P)
     CALL whm_drift(npl, whm_pl1P, dt)
     CALL coord_j2h(npl, whm_pl1P)
     CALL whm_getacch(lextra_force, t+dt, npl, nplmax, whm_pl1P, j2rp2, j4rp4)
     CALL whm_kickvh(npl, whm_pl1P, dth)

     RETURN

END SUBROUTINE whm_step_pl
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
