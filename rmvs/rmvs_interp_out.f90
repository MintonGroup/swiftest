!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_interp_out
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Interpolate planet positions between two Keplerian orbits in outer encounter region
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
!  Invocation  : CALL rmvs_interp_out(npl, rmvs_pl1P, dt)
!
!  Notes       : Adapted from Hal Levison's Swift routine rmvs3_interp.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_interp_out(npl, rmvs_pl1P, dt)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_interp_out
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(rmvs_pl), POINTER   :: rmvs_pl1P

! Internals
     INTEGER(I4B)              :: i, j, iflag
     REAL(DP)                  :: msun, dto, frac, dntenc
     REAL(DP), DIMENSION(NDIM) :: xtmp, vtmp
     TYPE(rmvs_pl), POINTER    :: rmvs_plP

! Executable code
     dntenc = REAL(NTENC, DP)
     msun = rmvs_pl1P%whm%swifter%mass
     dto = dt/dntenc
     rmvs_plP => rmvs_pl1P
     DO i = 2, npl
          rmvs_plP => rmvs_plP%nextP
          xtmp(:) = rmvs_plP%xout(:, 0)
          vtmp(:) = rmvs_plP%vout(:, 0)
          DO j = 1, NTENC - 1
               CALL drift_one(msun, xtmp(:), vtmp(:), dto, iflag)
               IF (iflag /= 0) THEN
                    WRITE(*, *) " Planet ", rmvs_plP%whm%swifter%id, " is lost!!!!!!!!!!"
                    WRITE(*, *) msun, dto
                    WRITE(*, *) xtmp(:)
                    WRITE(*, *) vtmp(:)
                    WRITE(*, *) " STOPPING "
                    CALL util_exit(FAILURE)
               END IF
               frac = 1.0_DP - j/dntenc
               rmvs_plP%xout(:, j) = frac*xtmp(:)
               rmvs_plP%vout(:, j) = frac*vtmp(:)
          END DO
          xtmp(:) = rmvs_plP%xout(:, NTENC)
          vtmp(:) = rmvs_plP%vout(:, NTENC)
          DO j = NTENC - 1, 1, -1
               CALL drift_one(msun, xtmp(:), vtmp(:), -dto, iflag)
               IF (iflag /= 0) THEN
                    WRITE(*, *) " Planet ", rmvs_plP%whm%swifter%id, " is lost!!!!!!!!!!"
                    WRITE(*, *) msun, -dto
                    WRITE(*, *) xtmp(:)
                    WRITE(*, *) vtmp(:)
                    WRITE(*, *) " STOPPING "
                    CALL util_exit(FAILURE)
               END IF
               frac = j/dntenc
               rmvs_plP%xout(:, j) = rmvs_plP%xout(:, j) + frac*xtmp(:)
               rmvs_plP%vout(:, j) = rmvs_plP%vout(:, j) + frac*vtmp(:)
          END DO
     END DO
     DO j = 0, NTENC
          rmvs_pl1P%xout(:, j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
          rmvs_pl1P%vout(:, j) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     END DO

     RETURN

END SUBROUTINE rmvs_interp_out
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
