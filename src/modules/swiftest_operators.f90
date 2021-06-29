module swiftest_operators
   !! author: David A. Minton
   !!
   !! Custom operators, including
   !!   A .cross. B = Cross product of A(1:3) and B(1:3) 
   !!
   !! Each operator can also do element-wise computation on arrays of the form .mag. A(1:3, 1:n)
   use swiftest_globals
   implicit none
   public

   !********************************************************************************************************************************
   ! Interfaces for .cross. operator
   !********************************************************************************************************************************

   interface operator(.cross.)
      module pure function operator_cross_sp(A, B) result(C)
         implicit none
         real(SP), dimension(:), intent(in) :: A, B
         real(SP), dimension(3) :: C
      end function operator_cross_sp

      module pure function operator_cross_dp(A, B) result(C)
         implicit none
         real(DP), dimension(:), intent(in) :: A, B
         real(DP), dimension(3) :: C
      end function operator_cross_dp

      module pure function operator_cross_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(:), intent(in) :: A, B
         integer(I1B), dimension(3) :: C
      end function operator_cross_i1b

      module pure function operator_cross_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(:), intent(in) :: A, B
         integer(I2B), dimension(3) :: C
      end function operator_cross_i2b

      module pure function operator_cross_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(:), intent(in) :: A, B
         integer(I4B), dimension(3) :: C
      end function operator_cross_i4b

      module pure function operator_cross_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(:), intent(in) :: A, B
         integer(I8B), dimension(3) :: C
      end function operator_cross_i8b

      module pure function operator_cross_el_sp(A, B) result(C)
         implicit none
         real(SP), dimension(:,:), intent(in) :: A, B
         real(SP), dimension(:,:), allocatable :: C
      end function operator_cross_el_sp

      module pure function operator_cross_el_dp(A, B) result(C)
         implicit none
         real(DP), dimension(:,:), intent(in) :: A, B
         real(DP), dimension(:,:), allocatable :: C
      end function operator_cross_el_dp

      module pure function operator_cross_el_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(:,:), intent(in) :: A, B
         integer(I1B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i1b

      module pure function operator_cross_el_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(:,:), intent(in) :: A, B
         integer(I2B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i2b

      module pure function operator_cross_el_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(:,:), intent(in) :: A, B
         integer(I4B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i4b

      module pure function operator_cross_el_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(:,:), intent(in) :: A, B
         integer(I8B), dimension(:,:), allocatable :: C
      end function operator_cross_el_i8b
   end interface

end module swiftest_operators
