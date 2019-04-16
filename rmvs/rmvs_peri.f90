!**********************************************************************************************************************************
!
!  Unit Name   : rmvs_peri
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : rmvs
!  Language    : Fortran 90/95
!
!  Description : Determine planetocentric pericenter passages for test particles in close encounters with a planet
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                index          : inner substep number within current set
!                nenc           : number of test particles encountering current planet
!                rmvs_pleP      : pointer to RMVS planet structure of planet being closely encountered
!                rmvs_tpenc1P   : pointer to RMVS test particle structure of first test particle encountering planet
!                mu             : mass of planet being encountered
!                rhill          : Hill sphere radius of planet being encountered
!                t              : time
!                dt             : time step
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : rmvs_tpenc1P   : pointer to RMVS test particle structure of first test particle encountering planet
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL rmvs_peri(lfirst, index, nenc, rmvs_pleP, rmvs_tpenc1P, mu, rhill, t, dt, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine util_peri.f
!
!**********************************************************************************************************************************
SUBROUTINE rmvs_peri(lfirst, index, nenc, rmvs_pleP, rmvs_tpenc1P, mu, rhill, t, dt, encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_rmvs
     USE module_interfaces, EXCEPT_THIS_ONE => rmvs_peri
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lfirst
     INTEGER(I4B), INTENT(IN) :: index, nenc
     REAL(DP), INTENT(IN)     :: mu, rhill, t, dt
     CHARACTER(*), INTENT(IN) :: encounter_file, out_type
     TYPE(rmvs_pl), POINTER   :: rmvs_pleP
     TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P

! Internals
     INTEGER(I4B)              :: i, id1, id2
     REAL(DP)                  :: r2, rhill2, vdotr, a, peri, capm, tperi, rpl
     REAL(DP), DIMENSION(NDIM) :: xh1, xh2, vh1, vh2
     TYPE(rmvs_tp), POINTER    :: rmvs_tpP

! Executable code
     rhill2 = rhill*rhill
     rmvs_tpP => rmvs_tpenc1P
     IF (lfirst) THEN
          DO i = 1, nenc
               IF (rmvs_tpP%whm%swifter%status == ACTIVE) THEN
                    vdotr = DOT_PRODUCT(rmvs_tpP%xpc(:), rmvs_tpP%vpc(:))
                    IF (vdotr > 0.0_DP) THEN
                         rmvs_tpP%isperi = 1
                    ELSE
                         rmvs_tpP%isperi = -1
                    END IF
               END IF
               rmvs_tpP => rmvs_tpP%tpencP
          END DO
     ELSE
          DO i = 1, nenc
               IF (rmvs_tpP%whm%swifter%status == ACTIVE) THEN
                    vdotr = DOT_PRODUCT(rmvs_tpP%xpc(:), rmvs_tpP%vpc(:))
                    IF (rmvs_tpP%isperi == -1) THEN
                         IF (vdotr >= 0.0_DP) THEN
                              rmvs_tpP%isperi = 0
                              CALL orbel_xv2aqt(rmvs_tpP%xpc(:), rmvs_tpP%vpc(:), mu, a, peri, capm, tperi)
                              r2 = DOT_PRODUCT(rmvs_tpP%xpc(:), rmvs_tpP%xpc(:))
                              IF ((ABS(tperi) > FACQDT*dt) .OR. (r2 > rhill2)) peri = SQRT(r2)
                              IF (encounter_file /= "") THEN
                                   id1 = rmvs_pleP%whm%swifter%id
                                   rpl = rmvs_pleP%whm%swifter%radius
                                   xh1(:) = rmvs_pleP%xin(:, index)
                                   vh1(:) = rmvs_pleP%vin(:, index)
                                   id2 = rmvs_tpP%whm%swifter%id
                                   xh2(:) = rmvs_tpP%whm%swifter%xh(:)
                                   vh2(:) = rmvs_tpP%whm%swifter%vh(:)
                                   CALL io_write_encounter(t, id1, id2, mu, 0.0_DP, rpl, 0.0_DP, xh1(:), xh2(:), vh1(:), vh2(:),  &
                                        encounter_file, out_type)
                              END IF
                              IF (rmvs_tpP%lperi) THEN
                                   IF (peri < rmvs_tpP%peri) THEN
                                        rmvs_tpP%peri = peri
                                        rmvs_tpP%plperP => rmvs_pleP
                                   END IF
                              ELSE
                                   rmvs_tpP%lperi = .TRUE.
                                   rmvs_tpP%peri = peri
                                   rmvs_tpP%plperP => rmvs_pleP
                              END IF
                         END IF
                    ELSE
                         IF (vdotr > 0.0_DP) THEN
                              rmvs_tpP%isperi = 1
                         ELSE
                              rmvs_tpP%isperi = -1
                         END IF
                    END IF
               END IF
               rmvs_tpP => rmvs_tpP%tpencP
          END DO                               
     END IF

     RETURN

END SUBROUTINE rmvs_peri
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
