!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_step_out2
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Step active test particles ahead in the outer encounter region, setting up and calling the inner region
!                integration if necessary
!
!  Input
!    Arguments : index          : outer substep number within current set
!                lfirst         : logical flag indicating whether current invocation is the first
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
!  Invocation  : CALL rmvs_step_out2(index, lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4,
!                                    dt, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_step_out2.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_step_out2(index, lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,       &
     encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_whm
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_step_out2
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lfirst, lextra_force
     INTEGER(I4B), INTENT(IN) :: index, npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4, dt
     CHARACTER(*), INTENT(IN) :: encounter_file, out_type
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER   :: rmvs_tp1P

! Internals
     LOGICAL(LGT)                                 :: lfirsttp, lencounter
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: rts
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xbeg, vbeg, xend
     TYPE(whm_pl), POINTER                        :: whm_pl1P
     TYPE(whm_tp), POINTER                        :: whm_tp1P, whm_tpP
     TYPE(rmvs_pl), POINTER                       :: rmvs_plP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xbeg(NDIM, nplmax), vbeg(NDIM, nplmax), xend(NDIM, nplmax))
          lmalloc = .FALSE.
     END IF
     lfirsttp = lfirst
     whm_pl1P => rmvs_pl1P%whm
     whm_tp1P => rmvs_tp1P%whm
     rmvs_plP => rmvs_pl1P
     DO i = 2, npl
          rmvs_plP => rmvs_plP%nextP
          xbeg(:, i) = rmvs_plP%xout(:, index-1)
          vbeg(:, i) = rmvs_plP%vout(:, index-1)
          xend(:, i) = rmvs_plP%xout(:, index)
     END DO
     rts = RHPSCALE
     CALL rmvs_chk(npl, ntp, rmvs_pl1P, rmvs_tp1P, xbeg, vbeg, dt, rts, lencounter)
     IF (lencounter) THEN
          rmvs_plP => rmvs_pl1P
          DO i = 2, npl
               rmvs_plP => rmvs_plP%nextP
               rmvs_plP%xin(:, 0) = xbeg(:, i)
               rmvs_plP%vin(:, 0) = vbeg(:, i)
               rmvs_plP%xin(:, NTPHENC) = xend(:, i)
               rmvs_plP%vin(:, NTPHENC) = rmvs_plP%vout(:, index)
          END DO
          CALL rmvs_interp_in(npl, rmvs_pl1P, dt)
          CALL rmvs_step_in(lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt, encounter_file,    &
               out_type)
          lfirsttp = .TRUE.
          CALL whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2, j4rp4, dt)
     ELSE
          CALL whm_step_tp(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2, j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE rmvs_step_out2
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
