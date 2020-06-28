submodule(swiftest_operators) operator_dot_implementation
   !! author: David A. Minton
   !!
   !! Contains implementations for the .dot. operator for all defined real types
   !! Single vector implementations:  C   = A(1:3) .dot. B(1:3)
   !! Vector list implementations:   C(:) = A(:, 1:3) .dot. B(::, 1:3)
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

   !> Internal elemental function that receives the component-wise input vector (SP)
   pure elemental function operator_dotvecelem_sp(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_sp)
      implicit none
      real(SP), intent(in)  :: Ax, Ay, Az
      real(SP), intent(in)  :: Bx, By, Bz
      real(SP)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_sp

   !> Internal elemental function that receives the component-wise input vector (DP)
   pure elemental function operator_dotvecelem_dp(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_dp)
      implicit none
      real(DP), intent(in)  :: Ax, Ay, Az
      real(DP), intent(in)  :: Bx, By, Bz
      real(DP)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_dp

   !> Internal elemental function that receives the component-wise input vector (I1B)
   pure elemental function operator_dotvecelem_i1b(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_i1b)
      implicit none
      integer(I1B), intent(in)  :: Ax, Ay, Az
      integer(I1B), intent(in)  :: Bx, By, Bz
      integer(I1B)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_i1b

   !> Internal elemental function that receives the component-wise input vector (I2B)
   pure elemental function operator_dotvecelem_i2b(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_i2b)
      implicit none
      integer(I2B), intent(in)  :: Ax, Ay, Az
      integer(I2B), intent(in)  :: Bx, By, Bz
      integer(I2B)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_i2b

   !> Internal elemental function that receives the component-wise input vector (I4B)
   pure elemental function operator_dotvecelem_i4b(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_i4b)
      implicit none
      integer(I4B), intent(in)  :: Ax, Ay, Az
      integer(I4B), intent(in)  :: Bx, By, Bz
      integer(I4B)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_i4b

   !> Internal elemental function that receives the component-wise input vector (I8B)
   pure elemental function operator_dotvecelem_i8b(Ax, Ay, Az, Bx, By, Bz) result(C)
      !$omp declare simd(operator_dotvecelem_i8b)
      implicit none
      integer(I8B), intent(in)  :: Ax, Ay, Az
      integer(I8B), intent(in)  :: Bx, By, Bz
      integer(I8B)           :: C
      C = Ax * Bx + Ay * By + Az * Bz 
      return
   end function operator_dotvecelem_i8b

   module procedure operator_dot_el_sp
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_sp(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_sp

   module procedure operator_dot_el_dp
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_dp(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_dp

   module procedure operator_dot_el_i1b
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_i1b(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_i1b

   module procedure operator_dot_el_i2b
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_i2b(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_i2b

   module procedure operator_dot_el_i4b
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_i4b(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_i4b

   module procedure operator_dot_el_i8b
      implicit none
      allocate(C(size(A,1)))
      C(:) = operator_dotvecelem_i8b(A(:, 1), A(:, 2), A(:, 3), B(:, 1), B(:, 2), B(:, 3))
      return
   end procedure operator_dot_el_i8b


end submodule operator_dot_implementation

