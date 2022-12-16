!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module swiftest
   !! author: David A. Minton
   !! graph: false
   !!
   !! This module serves to combine all of the Swiftest project modules under a single umbrella so that they can be accessed from individual submodule implementations with a simple "use swiftest" line.
   use swiftest_globals
   use swiftest_operators
   use lambda_function
   use swiftest_classes
   use whm_classes
   use rmvs_classes
   use helio_classes
   use symba_classes
   use encounter_classes
   use collision_classes
   use fraggle_classes
   use walltime_classes
   use io_progress_bar
   !use advisor_annotate
   !$ use omp_lib
   implicit none
   public

end module swiftest
