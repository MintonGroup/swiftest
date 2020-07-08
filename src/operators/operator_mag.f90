submodule(swiftest_operators) operator_mag_implementation
   !! author: David A. Minton
   !!
   !! Contains implementations for the .mag. operator for all defined real types
   !! Single vector implementations:  B   = .mag. A(1:3)
   !! Vector list implementations:   B(:) = .mag. A(:, 1:3)
   contains

   module procedure operator_mag_sp
      implicit none
      B = norm2(A(:))
      return
   end procedure operator_mag_sp

   module procedure operator_mag_dp
      implicit none
      B = norm2(A(:))
      return
   end procedure operator_mag_dp

   module procedure operator_mag_el_sp
      implicit none
      integer(I4B)  :: i,n
      n = size(A,1)
      allocate(B(n))
      do concurrent (i=1:n)
         B(i) = sqrt(A(i,1)**2 + A(i,2)**2 + A(i,3)**2) 
      end do
      return
   end procedure operator_mag_el_sp

   module procedure operator_mag_el_dp
      implicit none
      integer(I4B)  :: i,n
      n = size(A,1)
      allocate(B(n))
      do concurrent (i=1:n)
         B(i) = sqrt(A(i,1)**2 + A(i,2)**2 + A(i,3)**2) 
      end do
      return 
   end procedure operator_mag_el_dp

end submodule operator_mag_implementation

