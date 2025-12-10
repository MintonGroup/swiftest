! Copyright 2025 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module ringmoons
   !! author: Andrew Hesselbrock and David A. Minton
   !!
   !! Definition of classes and methods specific to the Ringmoons integrator
   use swiftest
   use symba
   implicit none
   public

   !> Ringmoons central body particle class
   type, extends(symba_cb) :: ringmoons_cb
   end type ringmoons_cb


   !> Ringmoons massive body class
   type, extends(symba_pl) :: ringmoons_pl
   end type ringmoons_pl


   !> Ringmoons test particle class
   type, extends(symba_tp) :: ringmoons_tp
   contains
   end type ringmoons_tp

   type, extends(base_object) :: ringmoons_ring
   contains
      procedure:: dealloc  => ringmoons_util_dealloc_ring
   end type ringmoons_ring

   type, extends(base_object) :: ringmoons_seeds
   contains
      procedure:: dealloc  => ringmoons_util_dealloc_seeds
   end type ringmoons_seeds

   type, extends(symba_nbody_system) :: ringmoons_nbody_system
      class(ringmoons_ring),         allocatable :: ring
         !! Ringmoons ring object
      class(ringmoons_seeds),        allocatable :: seeds
         !! Ringmoons seeds object
   end type ringmoons_nbody_system

   interface
      module subroutine ringmoons_util_dealloc_ring(self)
         !! author: David A. Minton
         !!
         !! Deallocates all allocatabale arrays
         implicit none
         ! Arguments
         class(ringmoons_ring),  intent(inout) :: self 
            !! Ringmoons ring object
      end subroutine ringmoons_util_dealloc_ring

      module subroutine ringmoons_util_dealloc_seeds(self)
         !! author: David A. Minton
         !!
         !! Deallocates all allocatabale arrays
         implicit none
         ! Arguments
         class(ringmoons_seeds),  intent(inout) :: self 
            !! Ringmoons ring object
      end subroutine ringmoons_util_dealloc_seeds
   end interface

end module ringmoons