!**********************************************************************************************************************************
!
!  Unit Name   : helio_lindrift_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Perform linear drift of test particles due to barycentric momentum of Sun
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                dt           : time step
!                pt           : negative barycentric velocity of the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_lindrift_tp(ntp, swifter_tp1P, dt, pt)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_lindrift_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_lindrift_tp(ntp, swifter_tp1P, dt, pt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => helio_lindrift_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)              :: ntp
     REAL(DP), INTENT(IN)                  :: dt
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: pt
     TYPE(swifter_tp), POINTER             :: swifter_tp1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_tp), POINTER :: swifter_tpP

! Executable code
     swifter_tpP => swifter_tp1P
     DO i = 1, ntp
          IF (swifter_tpP%status == ACTIVE) swifter_tpP%xh(:) = swifter_tpP%xh(:) + pt(:)*dt
          swifter_tpP => swifter_tpP%nextP
     END DO

     RETURN

END SUBROUTINE helio_lindrift_tp
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
