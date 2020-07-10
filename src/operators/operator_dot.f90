submodule(swiftest_operators) operator_dot_implementation
   !! author: David A. Minton
   !!
   !! Contains implementations for the .dot. operator for all defined real types
   !! Single vector implementations:  C   = A(1:3) .dot. B(1:3)
   !! Vector list implementations:   C(:) = A(1:3, :) .dot. B(1:3, :)
   contains

   module procedure operator_dot_sp
      implicit none
      C = dot_product(A(:), B(:))
      return
   end procedure operator_dot_sp

   module procedure operator_dot_dp
      implicit none
      C = dot_product(A(:), B(:))
      return
   end procedure operator_dot_dp

   module procedure operator_dot_el_sp
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_sp

   module procedure operator_dot_el_dp
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_dp

   module procedure operator_dot_el_i1b
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_i1b

   module procedure operator_dot_el_i2b
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_i2b

   module procedure operator_dot_el_i4b
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_i4b

   module procedure operator_dot_el_i8b
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C(n))
      do concurrent (i = 1:n)
         C(i) = dot_product(A(:, i), B(:, i))
      end do
      return
   end procedure operator_dot_el_i8b


end submodule operator_dot_implementation

