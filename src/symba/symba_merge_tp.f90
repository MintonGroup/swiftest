!**********************************************************************************************************************************
!
!  Unit Name   : symba_merge_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
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
SUBROUTINE symba_merge_tp(t, dt, index_enc, pltpenc_list, vbs, encounter_file, out_type, symba_plA, symba_tpA)

! Modules
     USE swiftest
     USE helio
     USE symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_merge_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA

! Internals
     LOGICAL(LGT)              :: lmerge
     INTEGER(I4B)              :: name1, name2, indexpl, indextp
     REAL(DP)                  :: r2, rlim, rlim2, vdotr, tcr2, dt2, mu, a, e, q, rad1
     REAL(DP), DIMENSION(NDIM) :: xr, vr, xh1, vh1, xh2, vh2

! Executable code
     lmerge = .FALSE.
     
     indexpl = pltpenc_list%indexpl(index_enc)
     indextp = pltpenc_list%indextp(index_enc)

     rlim = symba_plA%helio%swiftest%radius(indexpl)
     xr(:) = symba_tpA%helio%swiftest%xh(:,indextp) - symba_plA%helio%swiftest%xh(:,indexpl)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim*rlim
     IF (rlim2 >= r2) THEN
          lmerge = .TRUE.
     ELSE
          vr(:) = symba_tpA%helio%swiftest%vb(:,indextp) - symba_plA%helio%swiftest%vb(:,indexpl)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (pltpenc_list%lvdotr(index_enc) .AND. (vdotr > 0.0_DP)) THEN
               mu = symba_plA%helio%swiftest%mass(indexpl)
               tcr2 = r2/DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt*dt
               IF (tcr2 <= dt2) THEN
                    CALL orbel_xv2aeq(xr(:), vr(:), mu, a, e, q)
                    IF (q < rlim) lmerge = .TRUE.
               END IF
               IF (.NOT. lmerge) THEN
                    IF (encounter_file /= "") THEN
                         name1 = symba_plA%helio%swiftest%name(indexpl)
                         rad1 = symba_plA%helio%swiftest%radius(indexpl)
                         xh1(:) = symba_plA%helio%swiftest%xh(:,indexpl)
                         vh1(:) = symba_plA%helio%swiftest%vb(:,indexpl) - vbs(:)
                         name2 = symba_tpA%helio%swiftest%name(indextp)
                         xh2(:) = symba_tpA%helio%swiftest%xh(:,indextp)
                         vh2(:) = symba_tpA%helio%swiftest%vb(:,indextp) - vbs(:)
                         CALL io_write_encounter(t, name1, name2, mu, 0.0_DP, rad1, 0.0_DP, &
                              xh1(:), xh2(:), vh1(:), vh2(:), encounter_file, out_type)
                    END IF
               END IF
          END IF
     END IF
     IF (lmerge) THEN
          pltpenc_list%status(index_enc) = MERGED
          symba_tpA%helio%swiftest%status = DISCARDED_PLR
          WRITE(*, *) "Particle ", symba_tpA%helio%swiftest%name, " too close to Planet ", &
          symba_plA%helio%swiftest%name, " at t = ", t
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
