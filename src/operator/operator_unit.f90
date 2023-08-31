!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version NDIM of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(operators) s_operator_unit
   use, intrinsic :: ieee_exceptions
   !! author: David A. Minton
   !!
   !! Contains implementations for the .unit. operator for all defined real types
   !! Returns a unit vector or array of unit vectors from an input vector or array of vectors
   !! Single vector implementations:  B   = .unit. A(1:NDIM)
   !! Vector list implementations:   B(:) = .unit. A(1:NDIM, :)
contains

   pure module function operator_unit_sp(A) result(B)
      implicit none
      ! Arguments
      real(SP), dimension(:), intent(in) :: A
      real(SP), dimension(NDIM)          :: B
      ! Internals
      real(SP)  :: Amag

      call ieee_set_halting_mode(ieee_underflow, .false.)
      Amag = norm2(A(:))
      if (Amag > tiny(1._SP)) then
         B(:) = A(:) / Amag
      else
         B(:) = 0.0_SP
      end if

      return
   end function operator_unit_sp


   pure module function operator_unit_dp(A) result(B)
      implicit none
      ! Arguments
      real(DP), dimension(:), intent(in) :: A
      real(DP), dimension(NDIM)             :: B
      ! Internals
      real(DP)  :: Amag

      call ieee_set_halting_mode(ieee_underflow, .false.)
      Amag = norm2(A(:))
      if (Amag > tiny(1._DP)) then
         B(:) = A(:) / Amag
      else
         B(:) = 0.0_DP
      end if

      return
   end function operator_unit_dp

#ifdef QUADPREC
   pure module function operator_unit_qp(A) result(B)
      implicit none
      ! Arguments
      real(QP), dimension(:), intent(in) :: A
      real(QP), dimension(NDIM)             :: B
      ! Internals
      real(QP)  :: Amag

      call ieee_set_halting_mode(ieee_underflow, .false.)
      Amag = norm2(A(:))
      if (Amag > tiny(1._QP)) then
         B(:) = A(:) / Amag
      else
         B(:) = 0.0_QP
      end if

      return
   end function operator_unit_qp
#endif

   pure module function operator_unit_el_sp(A) result(B)
      implicit none
      ! Arguments
      real(SP), dimension(:,:), intent(in)   :: A
      real(SP), dimension(:,:), allocatable  :: B
      ! Internals
      integer(I4B)  :: i,n

      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(NDIM,n))

      do i=1,n
         B(:,i) = operator_unit_sp(A(:,i))
      end do

      return
   end function operator_unit_el_sp


   pure module function operator_unit_el_dp(A) result(B)
      implicit none
      ! Arguments
      real(DP), dimension(:,:), intent(in)   :: A
      real(DP), dimension(:,:), allocatable  :: B
      ! Internals
      integer(I4B)  :: i,n

      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(NDIM,n))

      do i=1,n
         B(:,i) = operator_unit_dp(A(:,i))
      end do

      return
   end function operator_unit_el_dp

#ifdef QUADPREC
   pure module function operator_unit_el_qp(A) result(B)
      implicit none
      ! Arguments
      real(QP), dimension(:,:), intent(in)   :: A
      real(QP), dimension(:,:), allocatable  :: B
      ! Internals
      integer(I4B)  :: i,n

      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(NDIM,n))

      do i=1,n
         B(:,i) = operator_unit_qp(A(:,i))
      end do

      return 
   end function operator_unit_el_qp
#endif

end submodule s_operator_unit

