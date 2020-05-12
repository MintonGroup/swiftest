!**********************************************************************************************************************************
!
!  Unit Name   : symba_fragmentation
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_fragmentation(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_fragmentation (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
     nplmax, ntpmax, fragmax)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_fragmentation
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
     REAL(DP), INTENT(INOUT)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA

! Internals
 
     INTEGER(I4B)                   :: model, nres
     REAL(DP), DIMENSION(3)         :: mres, rres
     REAL(DP), DIMENSION(NDIM, 3)   :: pres, vres
     INTEGER(I4B)                   :: regime, collresolve_resolve
     INTEGER(I4B)                   :: index1, index2
     INTEGER(I4B)                   :: name1, name2
     REAL(DP)                       :: r2, rlim, rlim2, vdotr, tcr2, dt2, mtot, a, e, q
     REAL(DP)                       :: rad1, rad2, m1, m2
     REAL(DP), DIMENSION(NDIM)      :: xr, vr, x1, v1, x2, v2
     LOGICAL(LGT)                   :: lfrag_add, lmerge


! Executable code

     lmerge = .FALSE.
     lfrag_add = .FALSE.
     ! Model 2 is the model for collresolve_resolve (LS12)
     model = 2

     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)

     rlim = symba_plA%helio%swiftest%radius(index1) + symba_plA%helio%swiftest%radius(index2)
     xr(:) = symba_plA%helio%swiftest%xh(:,index2) - symba_plA%helio%swiftest%xh(:,index1)
     r2 = DOT_PRODUCT(xr(:), xr(:))
     rlim2 = rlim*rlim
     ! checks if bodies are actively colliding in this time step
     IF (rlim2 >= r2) THEN 
          lfrag_add = .TRUE.
     ! if they are not actively colliding in  this time step, 
     !checks if they are going to collide next time step based on velocities and q
     ELSE 
          vr(:) = symba_plA%helio%swiftest%vb(:,index2) - symba_plA%helio%swiftest%vb(:,index1)
          vdotr = DOT_PRODUCT(xr(:), vr(:))
          IF (plplenc_list%lvdotr(index_enc) .AND. (vdotr > 0.0_DP)) THEN 
               tcr2 = r2/DOT_PRODUCT(vr(:), vr(:))
               dt2 = dt*dt
               IF (tcr2 <= dt2) THEN
                    mtot = symba_plA%helio%swiftest%mass(index1) + symba_plA%helio%swiftest%mass(index2)
                    CALL orbel_xv2aeq(xr(:), vr(:), mtot, a, e, q)
                    IF (q < rlim) lfrag_add = .TRUE.
               END IF
               ! if no collision is going to happen, write as close encounter, not  merger
               IF (.NOT. lfrag_add) THEN
                    IF (encounter_file /= "") THEN
                         name1 = symba_plA%helio%swiftest%name(index1)
                         m1 = symba_plA%helio%swiftest%mass(index1)
                         rad1 = symba_plA%helio%swiftest%radius(index1)
                         x1(:) = symba_plA%helio%swiftest%xh(:,index1)
                         v1(:) = symba_plA%helio%swiftest%vb(:,index1) - vbs(:)
                         name2 = symba_plA%helio%swiftest%name(index2)
                         m2 = symba_plA%helio%swiftest%mass(index2)
                         rad2 = symba_plA%helio%swiftest%radius(index2)
                         x2(:) = symba_plA%helio%swiftest%xh(:,index2)
                         v2(:) = symba_plA%helio%swiftest%vb(:,index2) - vbs(:)

                         CALL io_write_encounter(t, name1, name2, m1, m2, rad1, rad2, x1(:), x2(:), &
                              v1(:), v2(:), encounter_file, out_type)
                    END IF
               END IF
          END IF
     END IF

     IF (lfrag_add) THEN 
          regime = collresolve_resolve(model,m1,m2,rad1,rad2,x1(:),x2(:), v1(:),v2(:),nres, &
               mres,rres,pres,vres)
          !WRITE(*,*) "COLLISION REGIME = ", regime 
          CALL symba_caseresolve(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, &
               npltpenc, pltpenc_list, plplenc_list, regime, nplmax, ntpmax, fragmax, mres, rres)

     END IF 
     RETURN

END SUBROUTINE symba_fragmentation
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
