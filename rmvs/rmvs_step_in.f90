!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_step_in
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Step active test particles ahead in the inner encounter region
!
!  Input
!    Arguments : lextra_force   : logical flag indicating whether to include user-supplied accelerations
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
!    Arguments : rmvs_tp1P      : pointer to head of active RMVS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_step_in(lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,
!                                  encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_step_in.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_step_in(lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt, encounter_file,        &
     out_type)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_step_in
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
     REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4, dt
     CHARACTER(*), INTENT(IN) :: encounter_file, out_type
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
     TYPE(rmvs_tp), POINTER   :: rmvs_tp1P

! Internals
     LOGICAL(LGT)                                 :: lfirsttp
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, j, k, nenc
     REAL(DP)                                     :: mu, rhill, dti, time
     REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE    :: irh
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xh, aobl
     TYPE(swifter_pl), POINTER                    :: swifter_pl1P
     TYPE(rmvs_pl), POINTER                       :: rmvs_plP, rmvs_pleP
     TYPE(rmvs_tp), POINTER                       :: rmvs_tpP

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(irh(nplmax), xh(NDIM, nplmax), aobl(NDIM, nplmax))
          lmalloc = .FALSE.
     END IF
     dti = dt/NTPHENC
     IF (j2rp2 /= 0.0_DP) THEN
          swifter_pl1P => rmvs_pl1P%whm%swifter
          DO i = 0, NTPHENC
               rmvs_plP => rmvs_pl1P
               DO j = 2, npl
                    rmvs_plP => rmvs_plP%nextP
                    xh(:, j) = rmvs_plP%xin(:, i)
                    irh(j) = 1.0_DP/SQRT(DOT_PRODUCT(xh(:, j), xh(:, j)))
               END DO
               CALL obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
               rmvs_plP => rmvs_pl1P
               DO j = 1, npl
                    rmvs_plP%aobl(:, i) = aobl(:, j)
                    rmvs_plP => rmvs_plP%nextP
               END DO
          END DO
     END IF
     rmvs_pleP => rmvs_pl1P
     DO i = 2, npl
          rmvs_pleP => rmvs_pleP%nextP
          nenc = rmvs_pleP%nenc
          IF (nenc > 0) THEN
! There are inner encounters with this planet...switch to planetocentric coordinates to proceed
! Determine initial planetocentric positions and velocities for those test particles encountering this planet
               rmvs_tpP => rmvs_pleP%tpenc1P
               DO j = 1, nenc
                    rmvs_tpP%xpc(:) = rmvs_tpP%whm%swifter%xh(:) - rmvs_pleP%xin(:, 0)
                    rmvs_tpP%vpc(:) = rmvs_tpP%whm%swifter%vh(:) - rmvs_pleP%vin(:, 0)
                    rmvs_tpP => rmvs_tpP%tpencP
               END DO
! Determine planetocentric positions for planets and sun at all interpolated points in inner encounter
               rmvs_plP => rmvs_pl1P
               DO j = 1, npl
                    rmvs_plP%xpc(:, :) = rmvs_plP%xin(:, :) - rmvs_pleP%xin(:, :)
                    rmvs_plP => rmvs_plP%nextP
               END DO
               time = t
               mu = rmvs_pleP%whm%swifter%mass
               rhill = rmvs_pleP%whm%swifter%rhill
               CALL rmvs_peri(.TRUE., 0, nenc, rmvs_pleP, rmvs_pleP%tpenc1P, mu, rhill, time, dti, encounter_file, out_type)
! Now step the encountering test particles fully through the inner encounter
               lfirsttp = .TRUE.
               DO j = 1, NTPHENC
                    CALL rmvs_step_in_tp(j, lfirsttp, lextra_force, time, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, &
                         j4rp4, dti)
                    time = t + j*dti
                    CALL rmvs_peri(.FALSE., j, nenc, rmvs_pleP, rmvs_pleP%tpenc1P, mu, rhill, time, dti, encounter_file, out_type)
               END DO
               rmvs_tpP => rmvs_pleP%tpenc1P
               DO j = 1, nenc
                    IF (rmvs_tpP%whm%swifter%status == ACTIVE) rmvs_tpP%whm%swifter%status = INACTIVE
                    rmvs_tpP => rmvs_tpP%tpencP
               END DO
          END IF
     END DO

     RETURN

END SUBROUTINE rmvs_step_in
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
