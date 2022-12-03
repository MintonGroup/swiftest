!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_peri
   use swiftest
contains

   module subroutine util_peri_tp(self, system, param) 
      !! author: David A. Minton
      !!
      !! Determine system pericenter passages for test particles
      !! Note:  If the coordinate system used is barycentric, then this routine assumes that the barycentric coordinates in the
      !!        test particle structures are up-to-date and are not recomputed
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_peri.f90
      !! Adapted from Hal Levison's Swift routine util_peri.f
      implicit none
      ! Arguments
      class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: i
      real(DP)     :: e
      real(DP), dimension(:), allocatable :: vdotr

      associate(tp => self, ntp => self%nbody)
         allocate(vdotr(ntp))
         if (param%qmin_coord == "HELIO") then
            do i = 1, ntp
               vdotr(i) = dot_product(tp%rh(:, i), tp%vh(:, i))
               if (tp%isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     tp%isperi(i) = 0
                     call orbel_xv2aeq(tp%mu(i), tp%rh(1,i), tp%rh(2,i), tp%rh(3,i), tp%vh(1,i), tp%vh(2,i), tp%vh(3,i), &
                                       tp%atp(i), e, tp%peri(i))
                  end if
               else
                  if (vdotr(i) > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end do
         else
            do i = 1, ntp
               vdotr(i) = dot_product(tp%xb(:, i), tp%vb(:, i))
               if (tp%isperi(i) == -1) then
                  if (vdotr(i) >= 0.0_DP) then
                     tp%isperi(i) = 0
                     call orbel_xv2aeq(system%Gmtot, tp%xb(1,i), tp%xb(2,i), tp%xb(3,i), tp%vb(1,i), tp%vb(2,i), tp%vb(3,i), &
                                       tp%atp(i), e, tp%peri(i))
                  end if
               else
                  if (vdotr(i) > 0.0_DP) then
                     tp%isperi(i) = 1
                  else
                     tp%isperi(i) = -1
                  end if
               end if
            end do
         end if
      end associate

      return
   end subroutine util_peri_tp

end submodule s_util_peri
