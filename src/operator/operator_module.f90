! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module operators
   !! author: David A. Minton
   !!
   !! Custom operators, including
   !!   A .cross. B = Cross product of A(1:NDIM) and B(1:NDIM) 
   !!
   !! Each operator can also do element-wise computation on arrays of the form .mag. A(1:NDIM, 1:n)
   use globals
   implicit none
   public

   !********************************************************************************************************************************
   ! Interfaces for .cross. operator: Computes the cross product of two (NDIM) vectors or (NDIM,:) arrays 
   !********************************************************************************************************************************

   interface operator(.cross.)
      pure module function operator_cross_sp(A, B) result(C)
         !$omp declare simd(operator_cross_sp)
         implicit none
         real(SP), dimension(:), intent(in) :: A, B
         real(SP), dimension(NDIM) :: C
      end function operator_cross_sp

      pure module function operator_cross_dp(A, B) result(C)
         !$omp declare simd(operator_cross_dp)
         implicit none
         real(DP), dimension(:), intent(in) :: A, B
         real(DP), dimension(NDIM) :: C
      end function operator_cross_dp

#ifdef QUADPREC
      pure module function operator_cross_qp(A, B) result(C)
         !$omp declare simd(operator_cross_qp)
         implicit none
         real(QP), dimension(:), intent(in) :: A, B
         real(QP), dimension(NDIM) :: C
      end function operator_cross_qp
#endif

      pure module function operator_cross_i1b(A, B) result(C)
         !$omp declare simd(operator_cross_i1b)
         implicit none
         integer(I1B), dimension(:), intent(in) :: A, B
         integer(I1B), dimension(NDIM) :: C
      end function operator_cross_i1b

      pure module function operator_cross_i2b(A, B) result(C)
         !$omp declare simd(operator_cross_i2b)
         implicit none
         integer(I2B), dimension(:), intent(in) :: A, B
         integer(I2B), dimension(NDIM) :: C
      end function operator_cross_i2b

      pure module function operator_cross_i4b(A, B) result(C)
         !$omp declare simd(operator_cross_i4b)
         implicit none
         integer(I4B), dimension(:), intent(in) :: A, B
         integer(I4B), dimension(NDIM) :: C
      end function operator_cross_i4b

      pure module function operator_cross_i8b(A, B) result(C)
         !$omp declare simd(operator_cross_i8b)
         implicit none
         integer(I8B), dimension(:), intent(in) :: A, B
         integer(I8B), dimension(NDIM) :: C
      end function operator_cross_i8b

      pure module function operator_cross_el_sp(A, B) result(C)
         implicit none
         real(SP), dimension(:,:), intent(in) :: A, B
         real(SP), dimension(:,:), allocatable :: C
      end function operator_cross_el_sp

      pure module function operator_cross_el_dp(A, B) result(C)
         implicit none
         real(DP), dimension(:,:), intent(in) :: A, B
         real(DP), dimension(:,:), allocatable :: C
      end function operator_cross_el_dp

#ifdef QUADPREC
      pure module function operator_cross_el_qp(A, B) result(C)
         implicit none
         real(QP), dimension(:,:), intent(in) :: A, B
         real(QP), dimension(:,:), allocatable :: C
      end function operator_cross_el_qp
#endif

      pure module function operator_cross_el_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(:,:), intent(in) :: A, B
         integer(I1B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i1b

      pure module function operator_cross_el_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(:,:), intent(in) :: A, B
         integer(I2B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i2b

      pure module function operator_cross_el_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(:,:), intent(in) :: A, B
         integer(I4B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i4b

      pure module function operator_cross_el_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(:,:), intent(in) :: A, B
         integer(I8B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i8b
   end interface

   !********************************************************************************************************************************
   ! Interfaces for .mag. operator: Computes the magnitude of a vector or array of vectors using norm2
   !********************************************************************************************************************************

   interface operator(.mag.)
      pure module function operator_mag_sp(A) result(B)
         !$omp declare simd(operator_mag_sp)
         implicit none
         real(SP), dimension(:), intent(in) :: A
         real(SP)                           :: B
      end function operator_mag_sp

      pure module function operator_mag_dp(A) result(B)
         !$omp declare simd(operator_mag_dp)
         implicit none
         real(DP), dimension(:), intent(in) :: A
         real(DP)                           :: B
      end function operator_mag_dp

#ifdef QUADPREC
      pure module function operator_mag_qp(A) result(B)
         !$omp declare simd(operator_mag_qp)
         implicit none
         real(QP), dimension(:), intent(in) :: A
         real(QP)                           :: B
      end function operator_mag_qp
#endif
      pure module function operator_mag_el_sp(A) result(B)
         implicit none
         real(SP), dimension(:,:), intent(in) :: A
         real(SP), dimension(:), allocatable  :: B
      end function operator_mag_el_sp

      pure module function operator_mag_el_dp(A) result(B)
         implicit none
         real(DP), dimension(:,:), intent(in) :: A
         real(DP), dimension(:), allocatable  :: B
      end function operator_mag_el_dp

#ifdef QUADPREC
      pure module function operator_mag_el_qp(A) result(B)
         implicit none
         real(QP), dimension(:,:), intent(in) :: A
         real(QP), dimension(:), allocatable  :: B
      end function operator_mag_el_qp
#endif

   end interface


   !********************************************************************************************************************************
   ! Interfaces for .unit. operator: Returns a unit vector or array of unit vectors from an input vector or array of vectors
   !********************************************************************************************************************************

   interface operator(.unit.)
      pure module function operator_unit_sp(A) result(B)
         !$omp declare simd(operator_unit_sp)
         implicit none
         real(SP), dimension(:), intent(in)  :: A
         real(SP), dimension(NDIM)              :: B
      end function operator_unit_sp

      pure module function operator_unit_dp(A) result(B)
         !$omp declare simd(operator_unit_dp)
         implicit none
         real(DP), dimension(:), intent(in)  :: A
         real(DP), dimension(NDIM)              :: B
      end function operator_unit_dp

#ifdef QUADPREC
      pure module function operator_unit_qp(A) result(B)
         !$omp declare simd(operator_unit_qp)
         implicit none
         real(QP), dimension(:), intent(in)  :: A
         real(QP), dimension(NDIM)              :: B
      end function operator_unit_qp
#endif

      pure module function operator_unit_el_sp(A) result(B)
         implicit none
         real(SP), dimension(:,:), intent(in)  :: A
         real(SP), dimension(:,:), allocatable :: B
      end function operator_unit_el_sp

      pure module function operator_unit_el_dp(A) result(B)
         implicit none
         real(DP), dimension(:,:), intent(in) :: A
         real(DP), dimension(:,:), allocatable  :: B
      end function operator_unit_el_dp

#ifdef QUADPREC
      pure module function operator_unit_el_qp(A) result(B)
         implicit none
         real(QP), dimension(:,:), intent(in)  :: A
         real(QP), dimension(:,:), allocatable :: B
      end function operator_unit_el_qp
#endif
   end interface


end module operators
