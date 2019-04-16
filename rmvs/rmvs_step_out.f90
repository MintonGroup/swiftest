!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_step_out
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Step active test particles ahead in the outer encounter region
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                lextra_force   : logical flag indicating whether to include user-supplied accelerations
!                t              : time
!                npl            : number of planets
!                nplmax         : maximum allowed number of planets
!                ntp            : number of active test particles
!                ntpmax         : maximum allowed number of test particles
!                rmvs_pl1P      : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P      : pointer to head of active RMVS test particle structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_pl1P      : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P      : pointer to head of active RMVS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_step_out(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,
!                                   encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_step_out.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_step_out(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,               &
     encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_step_out
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lfirst, lextra_force
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     CHARACTER(*), INTENT(IN)    :: encounter_file, out_type
     TYPE(rmvs_pl), POINTER      :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER      :: rmvs_tp1P

! Internals
     LOGICAL(LGT)           :: lfirsttp
     INTEGER(I4B)           :: i, j, k, nenc
     REAL(DP)               :: dto, time
     TYPE(rmvs_pl), POINTER :: rmvs_pleP
     TYPE(rmvs_tp), POINTER :: rmvs_tpP

! Executable code
     lfirsttp = lfirst
     dto = dt/NTENC
     rmvs_tpP => rmvs_tp1P
     DO i = 1, ntp
          IF (.NOT. ASSOCIATED(rmvs_tpP%plencP)) THEN
               rmvs_tpP%whm%swifter%status = INACTIVE
          ELSE
               rmvs_tpP%lperi = .FALSE.
          END IF
          rmvs_tpP => rmvs_tpP%nextP
     END DO
     DO i = 1, NTENC
          time = t + (i - 1)*dto
          CALL rmvs_step_out2(i, lfirsttp, lextra_force, time, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dto, &
               encounter_file, out_type)
          rmvs_pleP => rmvs_pl1P
          DO j = 2, npl
               rmvs_pleP => rmvs_pleP%nextP
               nenc = rmvs_pleP%nenc
               IF (nenc > 0) THEN
                    rmvs_tpP => rmvs_pleP%tpenc1P
                    DO k = 1, nenc
                         IF (rmvs_tpP%whm%swifter%status == INACTIVE) rmvs_tpP%whm%swifter%status = ACTIVE
                         rmvs_tpP => rmvs_tpP%tpencP
                    END DO
               END IF
          END DO
     END DO

     RETURN

END SUBROUTINE rmvs_step_out
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
