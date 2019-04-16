!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in heliocentric coordinates
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
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                rmvs_pl1P      : pointer to head of RMVS planet structure linked-list
!                rmvs_tp1P      : pointer to head of active RMVS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,
!                               encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_step.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt, encounter_file,   &
     out_type)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)    :: lextra_force
     LOGICAL(LGT), INTENT(INOUT) :: lfirst
     INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
     CHARACTER(*), INTENT(IN)    :: encounter_file, out_type
     TYPE(rmvs_pl), POINTER      :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER      :: rmvs_tp1P

! Internals
     LOGICAL(LGT)                                 :: lfirstpl, lfirsttp, lencounter
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i
     REAL(DP)                                     :: rts
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xbeg, vbeg, xend
     TYPE(swifter_pl), POINTER                    :: swifter_plP
     TYPE(swifter_tp), POINTER                    :: swifter_tpP
     TYPE(whm_pl), POINTER                        :: whm_pl1P, whm_plP
     TYPE(whm_tp), POINTER                        :: whm_tp1P
     TYPE(rmvs_pl), POINTER                       :: rmvs_plP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xbeg(NDIM, nplmax), vbeg(NDIM, nplmax), xend(NDIM, nplmax))
          lmalloc = .FALSE.
     END IF
     whm_pl1P => rmvs_pl1P%whm
     whm_tp1P => rmvs_tp1P%whm
     swifter_plP => whm_pl1P%swifter
     DO i = 2, npl
          swifter_plP => swifter_plP%nextP
          xbeg(:, i) = swifter_plP%xh(:)
          vbeg(:, i) = swifter_plP%vh(:)
     END DO
     rts = RHSCALE
     CALL rmvs_chk(npl, ntp, rmvs_pl1P, rmvs_tp1P, xbeg, vbeg, dt, rts, lencounter)
     IF (lencounter) THEN
          rmvs_plP => rmvs_pl1P
          DO i = 2, npl
               rmvs_plP => rmvs_plP%nextP
               rmvs_plP%xout(:, 0) = xbeg(:, i)
               rmvs_plP%vout(:, 0) = vbeg(:, i)
          END DO
          lfirstpl = lfirst
          CALL whm_step_pl(lfirstpl, lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4, dt)
          rmvs_plP => rmvs_pl1P
          DO i = 2, npl
               rmvs_plP => rmvs_plP%nextP
               swifter_plP => rmvs_plP%whm%swifter
               rmvs_plP%xout(:, NTENC) = swifter_plP%xh(:)
               rmvs_plP%vout(:, NTENC) = swifter_plP%vh(:)
               xend(:, i) = swifter_plP%xh(:)
          END DO
          CALL rmvs_interp_out(npl, rmvs_pl1P, dt)
          CALL rmvs_step_out(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,           &
               encounter_file, out_type)
          swifter_tpP => whm_tp1P%swifter
          DO i = 1, ntp
               IF (swifter_tpP%status == ACTIVE) THEN
                    swifter_tpP%status = INACTIVE
               ELSE IF (swifter_tpP%status == INACTIVE) THEN
                    swifter_tpP%status = ACTIVE
               END IF
               swifter_tpP => swifter_tpP%nextP
          END DO
          lfirsttp = .TRUE.
          CALL whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2, j4rp4, dt)
          swifter_tpP => whm_tp1P%swifter
          DO i = 1, ntp
               IF (swifter_tpP%status == INACTIVE) swifter_tpP%status = ACTIVE
               swifter_tpP => swifter_tpP%nextP
          END DO
     ELSE
          CALL whm_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, j2rp2, j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE rmvs_step
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
