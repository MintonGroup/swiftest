!**********************************************************************************************************************************
!
!  Unit Name   : symba_merge_tp
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check for merger between planet and test particle in SyMBA
!
!  Input
!    Arguments : t              : time
!                dt             : time step
!                index          : index of planet-test particle encounter in array pltpenc_list
!                npltpenc       : number of planet-test particle encounters
!                pltpenc_list   : array of planet-test particle encounter structures
!                vbs            : barycentric velocity of the Sun
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : pltpenc_list   : array of planet-test particle encounter structures
!    Terminal  : status message
!    File      : none
!
!  Invocation  : CALL symba_merge_tp(t, dt, index, npltpenc, pltpenc_list, vbs, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_merge_tp(t, dt, index, npltpenc, pltpenc_list, vbs, encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_merge_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index, npltpenc
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list

! Internals
     LOGICAL(LGT)              :: lmerge
     INTEGER(I4B)              :: id1, id2
     REAL(DP)                  :: r2, rlim, rlim2, vdotr, tcr2, dt2, mu, a, e, q, rad1
     REAL(DP), DIMENSION(NDIM) :: xr, vr, xh1, vh1, xh2, vh2
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(swifter_tp), POINTER :: swifter_tpP
     TYPE(symba_pl), POINTER   :: symba_plP
     TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     lmerge = .FALSE.
     symba_plP => pltpenc_list(index)%plP
     symba_tpP => pltpenc_list(index)%tpP
     swifter_plP => symba_plP%helio%swifter
     swifter_tpP => symba_tpP%helio%swifter
     rlim = swifter_plP%radius
     xr(:) = swifter_tpP%xh(:) - swifter_plP%xh(:)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim*rlim
     IF (rlim2 >= r2) THEN
          lmerge = .TRUE.
     ELSE
          vr(:) = swifter_tpP%vb(:) - swifter_plP%vb(:)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (pltpenc_list(index)%lvdotr .AND. (vdotr > 0.0_DP)) THEN
               mu = swifter_plP%mass
               tcr2 = r2/DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt*dt
               IF (tcr2 <= dt2) THEN
                    CALL orbel_xv2aeq(xr(:), vr(:), mu, a, e, q)
                    IF (q < rlim) lmerge = .TRUE.
               END IF
               IF (.NOT. lmerge) THEN
                    IF (encounter_file /= "") THEN
                         id1 = swifter_plP%id
                         rad1 = swifter_plP%radius
                         xh1(:) = swifter_plP%xh(:)
                         vh1(:) = swifter_plP%vb(:) - vbs(:)
                         id2 = swifter_tpP%id
                         xh2(:) = swifter_tpP%xh(:)
                         vh2(:) = swifter_tpP%vb(:) - vbs(:)
                         CALL io_write_encounter(t, id1, id2, mu, 0.0_DP, rad1, 0.0_DP, xh1(:), xh2(:), vh1(:), vh2(:),           &
                              encounter_file, out_type)
                    END IF
               END IF
          END IF
     END IF
     IF (lmerge) THEN
          pltpenc_list(index)%status = MERGED
          swifter_tpP%status = DISCARDED_PLR
          WRITE(*, *) "Particle ", swifter_tpP%id, " too close to Planet ", swifter_plP%id, " at t = ", t
     END IF

     RETURN

END SUBROUTINE symba_merge_tp
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
