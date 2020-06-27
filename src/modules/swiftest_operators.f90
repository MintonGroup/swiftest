module swiftest_operators
   !! author: David A. Minton
   !!
   !! Custom operators, including
   !!   A .cross. B = Cross product of A(1:3) and B(1:3) 
   !!       .mag. A = Vector magnitude (performance tests indicate some compilers do a better job with norm2(A) vs sqrt(dot_product(A)))
   !!
   !! Each operator can also do element-wise computation on arrays of the form .mag. A(1:1:n, 3)
   use swiftest_globals
   implicit none
   public

   !********************************************************************************************************************************
   ! Interfaces for .cross. operator
   !********************************************************************************************************************************

   interface operator(.cross.)
      module pure function operator_cross_sp(A, B) result(C)
         implicit none
         real(SP), dimension(3), intent(in) :: A, B
         real(SP), dimension(3) :: C
      end function operator_cross_sp

      module pure function operator_cross_dp(A, B) result(C)
         implicit none
         real(DP), dimension(3), intent(in) :: A, B
         real(DP), dimension(3) :: C
      end function operator_cross_dp

      module pure function operator_cross_qp(A, B) result(C)
         implicit none
         real(QP), dimension(3), intent(in) :: A, B
         real(QP), dimension(3) :: C
      end function operator_cross_qp

      module pure function operator_cross_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(3), intent(in) :: A, B
         integer(I1B), dimension(3) :: C
      end function operator_cross_i1b

      module pure function operator_cross_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(3), intent(in) :: A, B
         integer(I2B), dimension(3) :: C
      end function operator_cross_i2b

      module pure function operator_cross_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(3), intent(in) :: A, B
         integer(I4B), dimension(3) :: C
      end function operator_cross_i4b

      module pure function operator_cross_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(3), intent(in) :: A, B
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

      module pure function operator_cross_el_qp(A, B) result(C)
         implicit none
         real(QP), dimension(:,:), intent(in) :: A, B
         real(QP), dimension(:,:), allocatable :: C
      end function operator_cross_el_qp

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

   interface operator(.mag.)
      module pure function operator_mag_sp(A) result(B)
         implicit none
         real(SP), dimension(:), intent(in) :: A
         real(SP)                           :: B
      end function operator_mag_sp

      module pure function operator_mag_dp(A) result(B)
         implicit none
         real(DP), dimension(:), intent(in) :: A
         real(DP)                           :: B
      end function operator_mag_dp

      module pure function operator_mag_qp(A) result(B)
         implicit none
         real(QP), dimension(:), intent(in) :: A
         real(QP)                           :: B
      end function operator_mag_qp

      module pure function operator_mag_el_sp(A) result(B)
         implicit none
         real(SP), dimension(:,:), intent(in)   :: A
         real(SP), dimension(:),   allocatable  :: B
      end function operator_mag_el_sp

      module pure function operator_mag_el_dp(A) result(B)
         implicit none
         real(DP), dimension(:,:), intent(in)   :: A
         real(DP), dimension(:),   allocatable  :: B
      end function operator_mag_el_dp

      module pure function operator_mag_el_qp(A) result(B)
         implicit none
         real(QP), dimension(:,:), intent(in)   :: A
         real(QP), dimension(:),   allocatable  :: B
      end function operator_mag_el_qp
   end interface

   interface operator(.dot.)
      module pure function operator_dot_sp(A, B) result(C)
         implicit none
         real(SP), dimension(:), intent(in)  :: A, B
         real(SP)                            :: C
      end function operator_dot_sp

      module pure function operator_dot_dp(A, B) result(C)
         implicit none
         real(DP), dimension(:), intent(in)  :: A, B
         real(DP)                            :: C
      end function operator_dot_dp

      module pure function operator_dot_qp(A, B) result(C)
         implicit none
         real(QP), dimension(:), intent(in)  :: A, B
         real(QP)                            :: C
      end function operator_dot_qp

      module pure function operator_dot_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(:), intent(in)  :: A, B
         integer(I1B)                            :: C
      end function operator_dot_i1b

      module pure function operator_dot_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(:), intent(in)  :: A, B
         integer(I2B)                            :: C
      end function operator_dot_i2b

      module pure function operator_dot_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(:), intent(in)  :: A, B
         integer(I4B)                            :: C
      end function operator_dot_i4b

      module pure function operator_dot_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(:), intent(in)  :: A, B
         integer(I8B)                            :: C
      end function operator_dot_i8b

      module pure function operator_dot_el_sp(A, B) result(C)
         implicit none
         real(SP), dimension(:,:), intent(in)  :: A, B
         real(SP), dimension(:),   allocatable :: C
      end function operator_dot_el_sp

      module pure function operator_dot_el_dp(A, B) result(C)
         implicit none
         real(DP), dimension(:,:), intent(in)  :: A, B
         real(DP), dimension(:),   allocatable :: C
      end function operator_dot_el_dp

      module pure function operator_dot_el_qp(A, B) result(C)
         implicit none
         real(QP), dimension(:,:), intent(in)  :: A, B
         real(QP), dimension(:),   allocatable :: C
      end function operator_dot_el_qp

      module pure function operator_dot_el_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(:,:), intent(in)  :: A, B
         integer(I1B), dimension(:),   allocatable :: C
      end function operator_dot_el_i1b

      module pure function operator_dot_el_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(:,:), intent(in)  :: A, B
         integer(I2B), dimension(:),   allocatable :: C
      end function operator_dot_el_i2b

      module pure function operator_dot_el_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(:,:), intent(in)  :: A, B
         integer(I4B), dimension(:),   allocatable :: C   
      end function operator_dot_el_i4b

      module pure function operator_dot_el_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(:,:), intent(in)  :: A, B
         integer(I8B), dimension(:),   allocatable :: C   
      end function operator_dot_el_i8b
   end interface

end module swiftest_operators
