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

      case (collresolve_regime_disruption)
         call symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, config%nplmax, config%ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

      case (collresolve_regime_supercatastrophic)
         call symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
         eoffset, vbs, symba_plA, nplplenc, &
         plplenc_list, config%nplmax, config%ntpmax, fragmax, mres, rres, m1, m2, rad1, &
         rad2, x1, x2, v1, v2)

      case (collresolve_regime_graze_and_merge)
         call symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         npl, symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
         x2, v1, v2)

      case (collresolve_regime_hit_and_run)
         call symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, &
         config%nplmax, config%ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

      case (collresolve_regime_merge)
         call symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         npl, symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, &
         x2, v1, v2)
      
      case default 
         write(*,*) "error in symba_caseresolve, no regime selected"
      end select


return
   end procedure symba_caseresolve
end submodule s_symba_caseresolve
