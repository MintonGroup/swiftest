!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_interp_in
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Interpolate planet positions between two Keplerian orbits in inner encounter region
!
!  Input
!    Arguments : npl       : number of planets
!                rmvs_pl1P : pointer to head of RMVS planet structure linked-list
!                dt        : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_pl1P : pointer to head of RMVS planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL rmvs_interp_in(npl, rmvs_pl1P, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_interp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_interp_in(npl, rmvs_pl1P, dt)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_interp_in
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P

! Internals
     INTEGER(I4B)              :: i, j, iflag
     REAL(DP)                  :: msun, dti, frac, dntphenc
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(rmvs_pl), POINTER    :: rmvs_plP

! Executable code
     dntphenc = REAL(NTPHENC, DP)
     msun = rmvs_pl1P%whm%swifter%mass
     dti = dt/dntphenc
     rmvs_plP => rmvs_pl1P
     DO i = 2, npl
          rmvs_plP => rmvs_plP%nextP
          xtmp(:) = rmvs_plP%xin(:, 0)
          vtmp(:) = rmvs_plP%vin(:, 0)
          DO j = 1, NTPHENC - 1
               CALL drift_one(msun, xtmp(:), vtmp(:), dti, iflag)
               IF (iflag /= 0) THEN
                    WRITE(*, *) " Planet ", rmvs_plP%whm%swifter%id, " is lost!!!!!!!!!!"
                    WRITE(*, *) msun, dti
                    WRITE(*, *) xtmp(:)
                    WRITE(*, *) vtmp(:)
                    WRITE(*, *) " STOPPING "
                    CALL util_exit(FAILURE)
               END IF
               frac = 1.0_DP - j/dntphenc
               rmvs_plP%xin(:, j) = frac*xtmp(:)
               rmvs_plP%vin(:, j) = frac*vtmp(:)
          END DO
          xtmp(:) = rmvs_plP%xin(:, NTPHENC)
          vtmp(:) = rmvs_plP%vin(:, NTPHENC)
          DO j = NTPHENC - 1, 1, -1
               CALL drift_one(msun, xtmp(:), vtmp(:), -dti, iflag)
               IF (iflag /= 0) THEN
                    WRITE(*, *) " Planet ", rmvs_plP%whm%swifter%id, " is lost!!!!!!!!!!"
                    WRITE(*, *) msun, -dti
                    WRITE(*, *) xtmp(:)
                    WRITE(*, *) vtmp(:)
                    WRITE(*, *) " STOPPING "
                    CALL util_exit(FAILURE)
               END IF
               frac = j/dntphenc
               rmvs_plP%xin(:, j) = rmvs_plP%xin(:, j) + frac*xtmp(:)
               rmvs_plP%vin(:, j) = rmvs_plP%vin(:, j) + frac*vtmp(:)
          END DO
     END DO
     DO j = 0, NTPHENC
          rmvs_pl1P%xin(:, j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          rmvs_pl1P%vin(:, j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO

     RETURN

END SUBROUTINE rmvs_interp_in
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
