module custom_operators
   !! author: David A. Minton
   !!
   !! Custom operators, including
   !!   A .cross. B = Cross product of A(1:3) and B(1:3) 
   !!    .mag. A    = Vector magnitude (performance tests indicate some compilers do a better job with norm2(A) vs sqrt(dot_product(A)))
   !!                 Can also do element-wise computation on arrays of the form .mag. A(1:3,1:n)
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

   interface operator(.mag.)
      module procedure vec_mag_sp
      module procedure vec_mag_dp
      module procedure vec_mag_qp
      module procedure vec_mag_el_sp
      module procedure vec_mag_el_dp
      module procedure vec_mag_el_qp
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

      pure function vec_mag_sp(A) result(B)
         implicit none
         real(SP), dimension(:), intent(in) :: A
         real(SP)                           :: B

         ! Pick your poison
         B = norm2(A(:))
         !y = sqrt(dot_product(x(:), x(:))  
         return
      end function vec_mag_sp

      pure function vec_mag_dp(x) result(y)
         implicit none
         real(DP), dimension(:), intent(in) :: x
         real(DP)                           :: y
   
         ! Pick your poison
         y = norm2(x(:))
         !y = sqrt(dot_product(x(:), x(:))  
         return
      end function vec_mag_dp

      pure function vec_mag_qp(x) result(y)
         implicit none
         real(QP), dimension(:), intent(in) :: x
         real(QP)                           :: y
      
         ! Pick your poison
         y = norm2(x(:))
         !y = sqrt(dot_product(x(:), x(:))  
         return
      end function vec_mag_qp


      ! The following allow the magnitude operator to work element-wise
      pure function vec_mag_el_sp(A) result(B)
         implicit none
         real(SP), dimension(:,:), intent(in) :: A
         real(SP), dimension(:), allocatable :: B
         B(:) = vecelem_sp(A(1,:), A(2, :), A(3, :))
         return
      end function vec_mag_el_sp

      pure elemental function vecelem_sp(Ax, Ay, Az) result(B)
         implicit none
         real(SP), intent(in) :: Ax, Ay, Az
         real(SP) :: B
         B = norm2([Ax, Ay, Az])
         return
      end function vecelem_sp

      pure function vec_mag_el_dp(A) result(B)
         implicit none
         real(DP), dimension(:,:), intent(in) :: A
         real(DP), dimension(:), allocatable :: B
         B(:) = vecelem_dp(A(1,:), A(2, :), A(3, :))
         return
      end function vec_mag_el_dp

      pure elemental function vecelem_dp(Ax, Ay, Az) result(B)
         implicit none
         real(DP), intent(in) :: Ax, Ay, Az
         real(DP) :: B
         B = norm2([Ax, Ay, Az])
         return
      end function vecelem_dp

      pure function vec_mag_el_qp(A) result(B)
         implicit none
         real(QP), dimension(:,:), intent(in) :: A
         real(QP), dimension(:), allocatable :: B
         B(:) = vecelem_qp(A(1,:), A(2, :), A(3, :))
         return
      end function vec_mag_el_qp

      pure elemental function vecelem_qp(Ax, Ay, Az) result(B)
         implicit none
         real(QP), intent(in) :: Ax, Ay, Az
         real(QP) :: B
         B = norm2([Ax, Ay, Az])
         return
      end function vecelem_qp



end module custom_operators
