!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(operators) s_operator_cross
   use, intrinsic :: ieee_exceptions
   use swiftest
   !! author: David A. Minton
   !!
   !! Contains implementations for the .cross. operator for all defined integer and real types
   !! Computes the cross product of two (3) vectors or (3,:) arrays 
   !! Single vector implementations: C(1:3)   = A(1:3)   .cross. B(1:3) 
   !! Vector list implementations:   C(1:3, :) = A(1:3, :) .cross. B(1:3, :)
contains

   pure module function operator_cross_sp(A, B) result(C)
      implicit none
      real(SP), dimension(:), intent(in) :: A, B
      real(SP), dimension(3) :: C

      call ieee_set_halting_mode(ieee_underflow, .false.)
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_sp

   pure module function operator_cross_dp(A, B) result(C)
      implicit none
      real(DP), dimension(:), intent(in) :: A, B
      real(DP), dimension(3) :: C

      call ieee_set_halting_mode(ieee_underflow, .false.)
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_dp

#ifdef QUADPREC
   pure module function operator_cross_qp(A, B) result(C)
      implicit none
      real(QP), dimension(:), intent(in) :: A, B
      real(QP), dimension(3) :: C

      call ieee_set_halting_mode(ieee_underflow, .false.)
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_qp
#endif

   pure module function operator_cross_i1b(A, B) result(C)
      implicit none
      integer(I1B), dimension(:), intent(in) :: A, B
      integer(I1B), dimension(3) :: C

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i1b

   pure module function operator_cross_i2b(A, B) result(C)
      implicit none
      integer(I2B), dimension(:), intent(in) :: A, B
      integer(I2B), dimension(3) :: C

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i2b

   pure module function operator_cross_i4b(A, B) result(C)
      implicit none
      integer(I4B), dimension(:), intent(in) :: A, B
      integer(I4B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i4b

   pure module function operator_cross_i8b(A, B) result(C)     
      implicit none
      integer(I8B), dimension(:), intent(in) :: A, B
      integer(I8B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i8b

   pure module function operator_cross_el_sp(A, B) result(C)
      implicit none
      real(SP), dimension(:,:), intent(in)  :: A, B
      real(SP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_sp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_sp

   pure module function operator_cross_el_dp(A, B) result(C)
      implicit none
      real(DP), dimension(:,:), intent(in)  :: A, B
      real(DP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_dp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_dp

#ifdef QUADPREC
   pure module function operator_cross_el_qp(A, B) result(C)
      implicit none
      real(QP), dimension(:,:), intent(in)  :: A, B
      real(QP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_qp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_qp
#endif 

   pure module function operator_cross_el_i1b(A, B) result(C)
      implicit none
      integer(I1B), dimension(:,:), intent(in)  :: A, B
      integer(I1B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_i1b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i1b

   pure module function operator_cross_el_i2b(A, B) result(C)
      implicit none
      integer(I2B), dimension(:,:), intent(in)  :: A, B
      integer(I2B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_i2b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i2b

   pure module function operator_cross_el_i4b(A, B) result(C)
      implicit none
      integer(I4B), dimension(:,:), intent(in)  :: A, B
      integer(I4B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_i4b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i4b

   pure module function operator_cross_el_i8b(A, B) result(C)
      implicit none
      integer(I8B), dimension(:,:), intent(in)  :: A, B
      integer(I8B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do i = 1,n 
         C(:,i) = operator_cross_i8b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i8b

end submodule s_operator_cross