!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(operators) s_operator_mag
   use, intrinsic :: ieee_exceptions
   !! author: David A. Minton
   !!
   !! Contains implementations for the .mag. operator for all defined real types
   !! Computes the magnitude of a vector or array of vectors using norm2
   !! Single vector implementations:  B   = .mag. A(1:3)
   !! Vector list implementations:   B(:) = .mag. A(1:3, :)
contains

   pure module function operator_mag_sp(A) result(B)
      implicit none
      real(SP), dimension(:), intent(in) :: A
      real(SP)                           :: B
      call ieee_set_halting_mode(ieee_underflow, .false.)
      B = norm2(A(:))
      return
   end function operator_mag_sp

   pure module function operator_mag_dp(A) result(B)
      implicit none
      real(DP), dimension(:), intent(in) :: A
      real(DP)                           :: B
      call ieee_set_halting_mode(ieee_underflow, .false.)
      B = norm2(A(:))
      return
   end function operator_mag_dp

   pure module function operator_mag_el_sp(A) result(B)
      implicit none
      real(SP), dimension(:,:), intent(in) :: A
      real(SP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      call ieee_set_halting_mode(ieee_underflow, .false.)
      do i = 1,n 
         B(i) = norm2(A(:, i)) 
      end do
      return
   end function operator_mag_el_sp

   pure module function operator_mag_el_dp(A) result(B)
      implicit none
      real(DP), dimension(:,:), intent(in) :: A
      real(DP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      call ieee_set_halting_mode(ieee_underflow, .false.)
      do i = 1,n 
         B(i) = norm2(A(:, i)) 
      end do
      return 
   end function operator_mag_el_dp

#ifdef QUADPREC
   pure module function operator_mag_el_qp(A) result(B)
      implicit none
      real(QP), dimension(:,:), intent(in) :: A
      real(QP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      call ieee_set_halting_mode(ieee_underflow, .false.)
      do i = 1,n 
         B(i) = norm2(A(:, i)) 
      end do
      return 
   end function operator_mag_el_qp
#endif

end submodule s_operator_mag

