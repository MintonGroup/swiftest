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


   module subroutine fraggle_io_encounter_dump(self, param)
      implicit none
      class(fraggle_storage(*)), intent(inout) :: self   !! Encounter storage object
      class(swiftest_parameters),          intent(inout) :: param  !! Current run configuration parameters 
   end subroutine fraggle_io_encounter_dump

   module subroutine fraggle_io_encounter_initialize_output(self, param)
      implicit none
      class(fraggle_io_encounter_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),             intent(in)    :: param   
   end subroutine fraggle_io_encounter_initialize_output

   module subroutine fraggle_io_encounter_write_frame(self, nc, param)
      implicit none
      class(fraggle_encounter_snapshot), intent(in)    :: self   !! Swiftest encounter structure
      class(encounter_io_parameters),    intent(inout) :: nc     !! Parameters used to identify a particular encounter io NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters
   end subroutine fraggle_io_encounter_write_frame

   module subroutine fraggle_io_log_generate(frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the fragment generation
      implicit none
      ! Arguments
      class(fraggle_fragments),   intent(in) :: frag
      ! Internals
      integer(I4B) :: i
      character(STRMAX) :: errmsg
      character(len=*), parameter :: fmtlabel = "(A14,10(ES11.4,1X,:))"

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle fragment generation results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN,    "(' dL_tot should be very small' )")
      write(LUN,fmtlabel) ' dL_tot      |', (.mag.(frag%Ltot_after(:) - frag%Ltot_before(:))) / (.mag.frag%Ltot_before(:))
      write(LUN,        "(' dE_tot should be negative and equal to Qloss' )")
      write(LUN,fmtlabel) ' dE_tot      |', (frag%Etot_after - frag%Etot_before) / abs(frag%Etot_before)
      write(LUN,fmtlabel) ' Qloss       |', -frag%Qloss / abs(frag%Etot_before)
      write(LUN,fmtlabel) ' dE - Qloss  |', (frag%Etot_after - frag%Etot_before + frag%Qloss) / abs(frag%Etot_before)
      write(LUN,        "(' -------------------------------------------------------------------------------------')")

      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_generate


   module subroutine fraggle_io_log_pl(pl, param)
      !! author: David A. Minton
      !!
      !! Writes a single message to the fraggle log file
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(in) :: pl    !! Swiftest massive body object (only the new bodies generated in a collision)
      class(swiftest_parameters), intent(in) :: param !! Current swiftest run configuration parameters
      ! Internals
      integer(I4B) :: i
      character(STRMAX) :: errmsg

      return
      667 continue
      write(*,*) "Error writing Fraggle message to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_pl


   module subroutine fraggle_io_log_regime(colliders, frag)
      !! author: David A. Minton
      !!
      !! Writes a log of the results of the collisional regime determination
      implicit none
      ! Arguments
      class(fraggle_colliders),   intent(in) :: colliders !! Fraggle collider system object
      class(fraggle_fragments),   intent(in) :: frag      !! Fraggle fragment object
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=FRAGGLE_LOG_OUT, status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *, err = 667, iomsg = errmsg)
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "           Fraggle collisional regime determination results"
      write(LUN, *) "--------------------------------------------------------------------"
      write(LUN, *) "True number of colliders : ",colliders%ncoll
      write(LUN, *) "Index list of true colliders  : ",colliders%idx(1:colliders%ncoll)
      select case(frag%regime) 
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
      write(LUN, *) "Energy loss                  : ", frag%Qloss
      write(LUN, *) "--------------------------------------------------------------------"
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing Fraggle regime information to log file: " // trim(adjustl(errmsg))
   end subroutine fraggle_io_log_regime

end submodule s_fraggle_io