!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Call discard routine to determine spilled test particles, then remove them from active list
!
!  Input
!    Arguments : t              : time
!                npl            : number of planets
!                ntp            : number of active test particles
!                nsp            : number of spilled test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                symba_tpd1P    : pointer to head of discard SyMBA test particle structure linked-list
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
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                symba_tpd1P    : pointer to head of discard SyMBA test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_tp(t, npl, ntp, nsp, symba_pl1P, symba_tp1P, symba_tpd1P, dt, rmin, rmax, rmaxu, qmin,
!                                      qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_tp(t, npl, ntp, nsp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord,       &
     qmin_alo, qmin_ahi, lclose, lrhill_present)

! Modules
     USE swiftest
     USE helio
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_tp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)     :: lclose, lrhill_present
     INTEGER(I4B), INTENT(IN)     :: npl
     INTEGER(I4B), INTENT(INOUT)  :: ntp, nsp
     REAL(DP), INTENT(IN)         :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
     CHARACTER(*), INTENT(IN)     :: qmin_coord
     TYPE(symba_pl), INTENT(INOUT):: symba_plA
     TYPE(symba_tp), INTENT(INOUT):: symba_tpA

! Internals
     LOGICAL(LGT)              :: lclosel = .FALSE.
     INTEGER(I4B)              :: i

! Executable code
     CALL discard_tpt, dt, npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, rmin, rmax, rmaxu, qmin, &
          qmin_alo, qmin_ahi, qmin_coord, lclosel,  &
          lrhill_present)
     RETURN

END SUBROUTINE symba_discard_tp
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
