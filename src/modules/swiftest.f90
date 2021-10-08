module swiftest
   !! author: David A. Minton
   !! graph: false
   !!
   !! This module serves to combine all of the Swiftest project modules under a single umbrella so that they can be accessed from individual submodule implementations with a simple "use swiftest" line.
   use swiftest_globals
   use swiftest_operators
   use swiftest_classes
   use whm_classes
   use rmvs_classes
   use helio_classes
   use symba_classes
   use fraggle_classes
   use lambda_function
   use walltime_classes
   use encounter_classes
   !use advisor_annotate
   !$ use omp_lib
   implicit none
   public

end module swiftest
