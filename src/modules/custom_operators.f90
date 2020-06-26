module custom_operators
   !! author: David A. Minton
   !!
   !! Custom operators, including
   !!    .cross.   = Cross product of two vectors 
   use swiftest_globals
   implicit none
   public

   interface operator(.cross.)
      module procedure cross_product_sp
      module procedure cross_product_dp
      module procedure cross_product_qp
      module procedure cross_product_i1b
      module procedure cross_product_i2b
      module procedure cross_product_i4b
      module procedure cross_product_i8b
   end interface

   contains
      pure function cross_product_sp(A, B) result(C)
         implicit none
         real(SP), dimension(3), intent(in) :: A, B
         real(SP), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_sp

      pure function cross_product_dp(A, B) result(C)
         implicit none
         real(DP), dimension(3), intent(in) :: A, B
         real(DP), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_dp

      pure function cross_product_qp(A, B) result(C)
         implicit none
         real(QP), dimension(3), intent(in) :: A, B
         real(QP), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_qp

      pure function cross_product_i1b(A, B) result(C)
         implicit none
         integer(I1B), dimension(3), intent(in) :: A, B
         integer(I1B), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_i1b

      pure function cross_product_i2b(A, B) result(C)
         implicit none
         integer(I2B), dimension(3), intent(in) :: A, B
         integer(I2B), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_i2b

      pure function cross_product_i4b(A, B) result(C)
         implicit none
         integer(I4B), dimension(3), intent(in) :: A, B
         integer(I4B), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_i4b

      pure function cross_product_i8b(A, B) result(C)
         implicit none
         integer(I8B), dimension(3), intent(in) :: A, B
         integer(I8B), dimension(3) :: C
         C(1) = A(2) * B(3) - A(3) * B(2)
         C(2) = A(3) * B(1) - A(1) * B(3)
         C(3) = A(1) * B(2) - A(2) * B(1)
         return
      end function cross_product_i8b


end module custom_operators
