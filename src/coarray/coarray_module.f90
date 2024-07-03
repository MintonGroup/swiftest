! Copyright 2024 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module coarray
   !! author: David A. Minton
   !!
   !! Utilities that are used for coarray test particles
   !!
   use globals
   implicit none
   public

   interface coclone
      module subroutine coarray_component_clone_char(var,src_img)
         implicit none
         character(len=*), intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine coarray_component_clone_char

      module subroutine coarray_component_clone_DP_arr1D(var,src_img)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine coarray_component_clone_DP_arr1D

      module subroutine coarray_component_clone_DP_vec2D(var,src_img)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine coarray_component_clone_DP_vec2D

      module subroutine coarray_component_clone_I4B_arr1D(var,src_img)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine coarray_component_clone_I4B_arr1D

      module subroutine coarray_component_clone_lgt_arr1D(var,src_img)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: src_img
      end subroutine coarray_component_clone_lgt_arr1D
   end interface


   interface cocollect
      module subroutine coarray_component_collect_DP_arr1D(var,dest_img)
         implicit none
         real(DP), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: dest_img
      end subroutine coarray_component_collect_DP_arr1D

      module subroutine coarray_component_collect_DP_vec2D(var,dest_img)
         implicit none
         real(DP), dimension(:,:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: dest_img
      end subroutine coarray_component_collect_DP_vec2D

      module subroutine coarray_component_collect_I4B_arr1D(var,dest_img)
         implicit none
         integer(I4B), dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: dest_img
      end subroutine coarray_component_collect_I4B_arr1D

      module subroutine coarray_component_collect_lgt_arr1D(var,dest_img)
         implicit none
         logical, dimension(:), allocatable, intent(inout) :: var
         integer(I4B), intent(in),optional :: dest_img
      end subroutine coarray_component_collect_lgt_arr1D

   end interface

end module coarray