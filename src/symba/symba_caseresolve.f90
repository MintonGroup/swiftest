submodule (symba) s_symba_caseresolve 
contains
   module procedure symba_caseresolve 
   !! author: Jennifer L. L. Pouplin and Carlisle A. Wishard
   !!
   !! Resolve which of the collision regimes to apply
use swiftest
implicit none

! executable code

      select case (regime)

      case (COLLRESOLVE_REGIME_DISRUPTION)
         call symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, param%nplmax, param%ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

      case (COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         call symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
         eoffset, vbs, symba_plA, nplplenc, &
         plplenc_list, param%nplmax, param%ntpmax, fragmax, mres, rres, m1, m2, rad1, &
         rad2, x1, x2, v1, v2)

      case (COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
         call symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         npl, symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
         x2, v1, v2)

      case (COLLRESOLVE_REGIME_HIT_AND_RUN)
         call symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, &
         param%nplmax, param%ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

      case (COLLRESOLVE_REGIME_MERGE)
         call symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         npl, symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
         x2, v1, v2)
      
      case default 
         write(*,*) "Error in symba_caseresolve, no regime selected"
      end select


return
   end procedure symba_caseresolve
end submodule s_symba_caseresolve
