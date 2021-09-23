submodule(swiftest_operators) s_operator_cross
   use swiftest
   !! author: David A. Minton
   !!
   !! Contains implementations for the .cross. operator for all defined integer and real types
   !! Single vector implementations: C(1:3)   = A(1:3)   .cross. B(1:3) 
   !! Vector list implementations:   C(1:3, :) = A(1:3, :) .cross. B(1:3, :)
contains

   module pure function operator_cross_sp(A, B) result(C)
      !$omp declare simd(operator_cross_sp)
      implicit none
      real(SP), dimension(:), intent(in) :: A, B
      real(SP), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_sp

   module pure function operator_cross_dp(A, B) result(C)
      !$omp declare simd(operator_cross_dp)
      implicit none
      real(DP), dimension(:), intent(in) :: A, B
      real(DP), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_dp

   module pure function operator_cross_qp(A, B) result(C)
      !$omp declare simd(operator_cross_qp)
      implicit none
      real(QP), dimension(:), intent(in) :: A, B
      real(QP), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_qp

   module pure function operator_cross_i1b(A, B) result(C)
      !$omp declare simd(operator_cross_i1b)
      implicit none
      integer(I1B), dimension(:), intent(in) :: A, B
      integer(I1B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i1b

   module pure function operator_cross_i2b(A, B) result(C)
      !$omp declare simd(operator_cross_i2b)
      implicit none
      integer(I2B), dimension(:), intent(in) :: A, B
      integer(I2B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i2b

   module pure function operator_cross_i4b(A, B) result(C)
      !$omp declare simd(operator_cross_i4b)
      implicit none
      integer(I4B), dimension(:), intent(in) :: A, B
      integer(I4B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i4b

   module pure function operator_cross_i8b(A, B) result(C)     
      !$omp declare simd(operator_cross_i8b)
      implicit none
      integer(I8B), dimension(:), intent(in) :: A, B
      integer(I8B), dimension(3) :: C
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end function operator_cross_i8b

   module pure function operator_cross_el_sp(A, B) result(C)
      implicit none
      real(SP), dimension(:,:), intent(in)  :: A, B
      real(SP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_sp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_sp

   module pure function operator_cross_el_dp(A, B) result(C)
      implicit none
      real(DP), dimension(:,:), intent(in)  :: A, B
      real(DP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_dp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_dp

   module pure function operator_cross_el_qp(A, B) result(C)
      implicit none
      real(QP), dimension(:,:), intent(in)  :: A, B
      real(QP), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_qp(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_qp

   module pure function operator_cross_el_i1b(A, B) result(C)
      implicit none
      integer(I1B), dimension(:,:), intent(in)  :: A, B
      integer(I1B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_i1b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i1b

   module pure function operator_cross_el_i2b(A, B) result(C)
      implicit none
      integer(I2B), dimension(:,:), intent(in)  :: A, B
      integer(I2B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_i2b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i2b

   module pure function operator_cross_el_i4b(A, B) result(C)
      implicit none
      integer(I4B), dimension(:,:), intent(in)  :: A, B
      integer(I4B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_i4b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i4b

   module pure function operator_cross_el_i8b(A, B) result(C)
      implicit none
      integer(I8B), dimension(:,:), intent(in)  :: A, B
      integer(I8B), dimension(:,:), allocatable :: C
      integer(I4B) :: i, n
      n = size(A, 2)
      if (allocated(C)) deallocate(C)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(:,i) = operator_cross_i8b(A(:,i), B(:,i))
      end do
      return
   end function operator_cross_el_i8b

end submodule s_operator_cross