!**********************************************************************************************************************************
!
!  Unit Name   : symba_caseresolve
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
!  Invocation  : CALL symba_caseresolve(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, regime, &
     nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_caseresolve
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
     REAL(DP), INTENT(IN)                             :: t, dt
     REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
     REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     INTEGER(I4B), INTENT(IN)                         :: regime
     INTEGER(I4B), DIMENSION(npl), INTENT(INOUT)      :: array_index1_child, array_index2_child

! Internals
 
     INTEGER(I4B)                                     :: model, nres
     REAL(DP)                                         :: pres, vres

! Executable code
          WRITE(*,*) "COLLISION REGIME = ", regime 
          SELECT CASE (regime)

          CASE (COLLRESOLVE_REGIME_DISRUPTION)
               CALL symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
               CALL symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, &
               plplenc_list, nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, &
               rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
               CALL symba_casemerge (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_HIT_AND_RUN)
               CALL symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)

          CASE (COLLRESOLVE_REGIME_MERGE)
               CALL symba_casemerge (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
               encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          
          CASE DEFAULT 
               WRITE(*,*) "ERROR IN SYMBA_CASERESOLVE, NO REGIME SELECTED"
          END SELECT


RETURN
END SUBROUTINE symba_caseresolve
