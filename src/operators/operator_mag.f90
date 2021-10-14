submodule(swiftest_operators) s_operator_mag
   !! author: David A. Minton
   !!
   !! Contains implementations for the .mag. operator for all defined real types
   !! Single vector implementations:  B   = .mag. A(1:3)
   !! Vector list implementations:   B(:) = .mag. A(1:3, :)
   contains

   module pure function operator_mag_sp(A) result(B)
      implicit none
      real(SP), dimension(:), intent(in) :: A
      real(SP)                           :: B
      B = norm2(A(:))
      return
   end function operator_mag_sp

   module pure function operator_mag_dp(A) result(B)
      implicit none
      real(DP), dimension(:), intent(in) :: A
      real(DP)                           :: B
      B = norm2(A(:))
      return
   end function operator_mag_dp

   module pure function operator_mag_el_sp(A) result(B)
      implicit none
      real(SP), dimension(:,:), intent(in) :: A
      real(SP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      do concurrent (i=1:n)
         B(i) = norm2(A(:, i)) 
      end do
      return
   end function operator_mag_el_sp

   module pure function operator_mag_el_dp(A) result(B)
      implicit none
      real(DP), dimension(:,:), intent(in) :: A
      real(DP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      do concurrent (i=1:n)
         B(i) = norm2(A(:, i)) 
      end do
      return 
   end function operator_mag_el_dp

   module pure function operator_mag_el_qp(A) result(B)
      implicit none
      real(QP), dimension(:,:), intent(in) :: A
      real(QP), dimension(:), allocatable  :: B
      integer(I4B)  :: i,n
      n = size(A, 2)
      if (allocated(B)) deallocate(B)
      allocate(B(n))
      do concurrent (i=1:n)
         B(i) = norm2(A(:, i)) 
      end do
      return 
   end function operator_mag_el_qp

end submodule s_operator_mag

