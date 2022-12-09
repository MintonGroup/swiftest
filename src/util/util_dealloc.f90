!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_util_dealloc
   use swiftest
contains

   module subroutine util_dealloc_body(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_body),  intent(inout) :: self

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%lmask)) deallocate(self%lmask)
      if (allocated(self%mu)) deallocate(self%mu)
      if (allocated(self%rh)) deallocate(self%rh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%xb)) deallocate(self%xb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%ah)) deallocate(self%ah)
      if (allocated(self%aobl)) deallocate(self%aobl)
      if (allocated(self%agr)) deallocate(self%agr)
      if (allocated(self%atide)) deallocate(self%atide)
      if (allocated(self%ir3h)) deallocate(self%ir3h)
      if (allocated(self%a)) deallocate(self%a)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%inc)) deallocate(self%inc)
      if (allocated(self%capom)) deallocate(self%capom)
      if (allocated(self%omega)) deallocate(self%omega)
      if (allocated(self%capm)) deallocate(self%capm)

      return
   end subroutine util_dealloc_body


   module subroutine util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest massive body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_pl),  intent(inout) :: self !! Swiftest massive body object

      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%Gmass)) deallocate(self%Gmass)
      if (allocated(self%rhill)) deallocate(self%rhill)
      if (allocated(self%renc)) deallocate(self%renc)
      if (allocated(self%radius)) deallocate(self%radius)
      if (allocated(self%density)) deallocate(self%density)
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%Ip)) deallocate(self%Ip)
      if (allocated(self%k2)) deallocate(self%k2)
      if (allocated(self%Q)) deallocate(self%Q)
      if (allocated(self%tlag)) deallocate(self%tlag)
      if (allocated(self%k_plpl)) deallocate(self%k_plpl)

      call util_dealloc_body(self)

      return
   end subroutine util_dealloc_pl

   module subroutine util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest test particle object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_tp),  intent(inout) :: self !! Swiftest test particle object

      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)
      if (allocated(self%k_pltp)) deallocate(self%k_pltp)

      call util_dealloc_body(self)

      return
   end subroutine util_dealloc_tp

end submodule s_util_dealloc