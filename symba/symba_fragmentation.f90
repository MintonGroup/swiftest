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
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA


! Internals
 
     INTEGER(I4B)                   :: model, nres, i
     REAL(DP), DIMENSION(3)         :: mres, rres
     REAL(DP), DIMENSION(NDIM, 3)   :: pres, vres
     INTEGER(I4B)                   :: regime, collresolve_resolve, regime2
     INTEGER(I4B)                   :: index1, index2, index1_child, index2_child, index1_parent, index2_parent
     INTEGER(I4B)                   :: name1, name2, index_big1, index_big2, stat1, stat2
     REAL(DP)                       :: r2, rlim, rlim2, vdotr, tcr2, dt2, a, e, q
     REAL(DP)                       :: rad1, rad2, m1, m2, GU, den1, den2, denchild
     REAL(DP)                       :: m1_cgs, m2_cgs, rad1_cgs, rad2_cgs, mass1, mass2, mmax, mtmp, mtot
     REAL(DP), DIMENSION(NDIM)      :: xr, vr, x1, v1, x2, v2
     REAL(DP), DIMENSION(NDIM)      :: x1_cgs, x2_cgs, v1_cgs, v2_cgs, x1_au, x2_au, v1_auy, v2_auy
     LOGICAL(LGT)                   :: lfrag_add, lmerge
     INTEGER(I4B), DIMENSION(npl)   :: array_index1_child, array_index2_child
     REAL(DP)                       :: MSUN, K2, m1_msun, m2_msun, rad1_au, rad2_au, AU2CM, year


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

     nres = 2
     IF (lfrag_add) THEN 
          symba_plA%lmerged(index1) = .TRUE.
          symba_plA%lmerged(index2) = .TRUE.
          index1_parent = symba_plA%index_parent(index1)
          m1 = symba_plA%helio%swiftest%mass(index1_parent)
          mass1 = m1 
          rad1 = symba_plA%helio%swiftest%radius(index1_parent)
          x1(:) = m1*symba_plA%helio%swiftest%xh(:,index1_parent)
          v1(:) = m1*symba_plA%helio%swiftest%vb(:,index1_parent)
          mmax = m1
          name1 = symba_plA%helio%swiftest%name(index1_parent)
          index_big1 = index1_parent
          stat1 = symba_plA%helio%swiftest%status(index1_parent)
          array_index1_child(1:npl) = symba_plA%index_child(1:npl,index1_parent)
          
          den1 = (m1**2) / ((4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index1_parent)**3)
          DO i = 1, symba_plA%nchild(index1_parent) ! initialize an array of children
               index1_child = array_index1_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index1_child)
               denchild = (mtmp**2) / ((4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index1_child)**3)
               den1 = den1 + denchild
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name1 = symba_plA%helio%swiftest%name(index1_child)
                    index_big1 = index1_child
                    stat1 = symba_plA%helio%swiftest%status(index1_child)
               END IF
               m1 = m1 + mtmp
               x1(:) = x1(:) + mtmp*symba_plA%helio%swiftest%xh(:,index1_child)
               v1(:) = v1(:) + mtmp*symba_plA%helio%swiftest%vb(:,index1_child)
          END DO
          den1 = den1 / m1
          rad1 = ((3.0_DP * m1) / (den1 * 4.0_DP * PI)) ** (1.0_DP / 3.0_DP)
          x1(:) = x1(:)/m1
          v1(:) = v1(:)/m1

          index2_parent = symba_plA%index_parent(index2)
          m2 = symba_plA%helio%swiftest%mass(index2_parent)
          mass2 = m2
          rad2 = symba_plA%helio%swiftest%radius(index2_parent)
          x2(:) = m2*symba_plA%helio%swiftest%xh(:,index2_parent)
          v2(:) = m2*symba_plA%helio%swiftest%vb(:,index2_parent)
          mmax = m2
          name2 = symba_plA%helio%swiftest%name(index2_parent)
          index_big2 = index2_parent
          stat2 = symba_plA%helio%swiftest%status(index2_parent)
          array_index2_child(1:npl) = symba_plA%index_child(1:npl,index2_parent)

          den2 = (m2**2) / ((4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index2_parent)**3)
          DO i = 1, symba_plA%nchild(index2_parent)
               index2_child = array_index2_child(i)
               mtmp = symba_plA%helio%swiftest%mass(index2_child)
               denchild = (mtmp**2) / ((4.0_DP / 3.0_DP) * PI * symba_plA%helio%swiftest%radius(index2_child)**3)
               den2 = den2 + denchild
               IF (mtmp > mmax) THEN
                    mmax = mtmp
                    name2 = symba_plA%helio%swiftest%name(index2_child)
                    index_big2 = index2_child
                    stat2 = symba_plA%helio%swiftest%status(index2_child)
               END IF
               m2 = m2 + mtmp
               x2(:) = x2(:) + mtmp*symba_plA%helio%swiftest%xh(:,index2_child)
               v2(:) = v2(:) + mtmp*symba_plA%helio%swiftest%vb(:,index2_child)
          END DO
          den2 = den2 / m2
          rad2 = ((3.0_DP * m2) / (den2 * 4.0_DP * PI)) ** (1.0_DP / 3.0_DP)
          x2(:) = x2(:)/m2
          v2(:) = v2(:)/m2

          GU = GC / (DU2CM**3 / (MU2GM * TU2S**2))
          m1_cgs = (m1 / GU) * MU2GM 
          m2_cgs = (m2 / GU) * MU2GM 
          MSUN = 1.989e33 !Msun in cgs
          AU2CM = 1.496e+13 !AU in cgs
          m1_msun = m1_cgs / MSUN
          m2_msun = m2_cgs / MSUN
          K2 = 2.959122082855911e-4
          rad1_cgs = (rad1) * DU2CM
          rad2_cgs = (rad2) * DU2CM
          rad1_au = rad1_cgs / AU2CM
          rad2_au = rad2_cgs / AU2CM
          x1_cgs(:) = x1(:) * DU2CM
          x2_cgs(:) = x2(:) * DU2CM
          x1_au(:) = x1_cgs(:) / AU2CM
          x2_au(:) = x2_cgs(:) / AU2CM

          v1_cgs(:) = v1(:) * DU2CM / TU2S 
          v2_cgs(:) = v2(:) * DU2CM / TU2S
          year = 3.154e7
          v1_auy(:) = v1_cgs(:) / AU2CM * (year)
          v2_auy(:) = v2_cgs(:) / AU2CM * (year)

          mres(:) = 0.0_DP
          rres(:) = 0.0_DP
          pres(:,:) = 0.0_DP
          vres(:,:) = 0.0_DP

          ! PROBLEM
          WRITE(*,*) "model: ", model
          WRITE(*,*) "m1_msun: ", m1_msun
          WRITE(*,*) "m2_msun: ", m2_msun
          WRITE(*,*) "rad1_au: ", rad1_au
          WRITE(*,*) "rad2_au: ", rad2_au
          WRITE(*,*) "x1_au: ", x1_au
          WRITE(*,*) "x2_au: ", x2_au
          WRITE(*,*) "v1_auy: ", v1_auy
          WRITE(*,*) "v2_auy: ", v2_auy
          WRITE(*,*) "nres: ", nres 
          WRITE(*,*) "mres: ", mres(:) !THIS IS THE PROBLEM
          WRITE(*,*) "rres: ", rres(:)
          WRITE(*,*) "pres: ", pres(:,:)
          WRITE(*,*) "vres: ", vres(:,:)

          WRITE(*,*) "Before collresolve_resolve"

          regime = collresolve_resolve(model,m1_msun,m2_msun,rad1_au,rad2_au,x1_au(:),x2_au(:), v1_auy(:),v2_auy(:), &
               nres,mres,rres,pres,vres)

          WRITE(*,*) "After collresolve_resolve"

          WRITE(*,*) "nres: ", nres 
          WRITE(*,*) "mres: ", mres(:) !THIS IS THE PROBLEM
          WRITE(*,*) "rres: ", rres(:)
          WRITE(*,*) "pres: ", pres(:,:)
          WRITE(*,*) "vres: ", vres(:,:)


          !PROBLEM

          mres(:) = mres(:)*GU*MSUN/MU2GM
          rres(:) = rres(:)*AU2CM/DU2CM

          CALL symba_caseresolve(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, &
               npltpenc, pltpenc_list, plplenc_list, regime, nplmax, ntpmax, fragmax, mres, rres, &
               array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

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
