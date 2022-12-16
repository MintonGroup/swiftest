!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(fraggle_classes) s_fraggle_io
   use swiftest

contains

   module subroutine fraggle_io_log_regime(impactors, fragments)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(collision_impactors),   intent(in) :: impactors !! Fraggle collider system object
      class(fraggle_fragments),   intent(in) :: fragments      !! Fraggle fragment object
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle collisional regime determination results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "True number of impactors : ",impactors%ncoll
      write(LUN, *) "Index list of true impactors  : ",impactors%idx(1:impactors%ncoll)
      select case(impactors%regime) 
      case(COLLRESOLVE_REGIME_MERGE)
         write(LUN, *) "Regime: Merge"
      case(COLLRESOLVE_REGIME_DISRUPTION)
         write(LUN, *) "Regime: Disruption"
      case(COLLRESOLVE_REGIME_SUPERCATASTROPHIC)
         write(LUN, *) "Regime: Supercatastrophic disruption"
      case(COLLRESOLVE_REGIME_GRAZE_AND_MERGE)
         write(LUN, *) "Regime: Graze and merge"
      case(COLLRESOLVE_REGIME_HIT_AND_RUN)
         write(LUN, *) "Regime: Hit and run"
      end select
      write(LUN, *) "Energy loss                  : ", impactors%Qloss
      write(LUN, *) "--------------------------------------------------------------------"
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle regime information to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_regime

end submodule s_fraggle_io