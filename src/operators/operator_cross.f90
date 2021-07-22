submodule(swiftest_operators) operator_cross_implementation
   use swiftest
   !! author: David A. Minton
   !!
   !! Contains implementations for the .cross. operator for all defined integer and real types
   !! Single vector implementations: C(1:3)   = A(1:3)   .cross. B(1:3) 
   !! Vector list implementations:   C(1:3, :) = A(1:3, :) .cross. B(1:3, :)
contains

   module procedure operator_cross_sp
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_sp

   module procedure operator_cross_dp      
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_dp

   module procedure operator_cross_i1b      
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_i1b

   module procedure operator_cross_i2b      
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_i2b

   module procedure operator_cross_i4b      
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_i4b

   module procedure operator_cross_i8b      
      implicit none
      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)
      return
   end procedure operator_cross_i8b

   module procedure operator_cross_el_sp
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end procedure operator_cross_el_sp

   module procedure operator_cross_el_dp      
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end procedure operator_cross_el_dp

   module procedure operator_cross_el_i1b  
      implicit none    
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end procedure operator_cross_el_i1b

   module procedure operator_cross_el_i2b      
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end procedure operator_cross_el_i2b

   module procedure operator_cross_el_i4b      
      implicit none
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end procedure operator_cross_el_i4b

   module pure function operator_cross_el_i8b(A, B) result(C)
      implicit none
      ! Arguments
      integer(I8B), dimension(:,:), intent(in)  :: A, B
      integer(I8B), dimension(:,:), allocatable :: C
      ! Internals
      integer(I4B) :: i, n
      n = size(A, 2)
      allocate(C, mold = A)
      do concurrent (i = 1:n) 
         C(1, i) = A(2, i) * B(3, i) - A(3, i) * B(2, i)
         C(2, i) = A(3, i) * B(1, i) - A(1, i) * B(3, i)
         C(3, i) = A(1, i) * B(2, i) - A(2, i) * B(1, i)
      end do
      return
   end function operator_cross_el_i8b

end submodule operator_cross_implementation