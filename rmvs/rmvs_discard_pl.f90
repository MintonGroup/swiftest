!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_discard_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Check to see if test particles should be discarded based on pericenter passage distances with respect to
!                planets encountered
!
!  Input
!    Arguments : t         : time
!                ntp       : number of active test particles
!                rmvs_tp1P : pointer to head of active RMVS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_tp1P : pointer to head of active RMVS test particle structure linked-list
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL rmvs_discard_pl(t, ntp, rmvs_tp1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_discard_pl(t, ntp, rmvs_tp1P)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: ntp
     REAL(DP), INTENT(IN)     :: t
     TYPE(rmvs_tp), POINTER   :: rmvs_tp1P

! Internals
     INTEGER(I4B)           :: i
     TYPE(rmvs_tp), POINTER :: rmvs_tpP

! Executable code
     rmvs_tpP => rmvs_tp1P
     DO i = 1, ntp
          IF (rmvs_tpP%whm%swifter%status == ACTIVE) THEN
               IF (rmvs_tpP%lperi) THEN
                    IF (rmvs_tpP%peri < rmvs_tpP%plperP%whm%swifter%radius) THEN
                         rmvs_tpP%whm%swifter%status = DISCARDED_PLQ
                         WRITE(*, *) "Particle ", rmvs_tpP%whm%swifter%id, " q with respect to Planet ",                          &
                              rmvs_tpP%plperP%whm%swifter%id, " is too small at t = ", t
                    END IF
               END IF
          END IF
          rmvs_tpP => rmvs_tpP%nextP
     END DO

     RETURN

END SUBROUTINE rmvs_discard_pl
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
