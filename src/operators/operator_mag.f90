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

   !> Internal elemental function that receives the component-wise input vector (SP)
   pure elemental function operator_vecelem_sp(Ax, Ay, Az) result(B)
      !$omp declare simd(operator_vecelem_sp)
      implicit none
      real(SP), intent(in) :: Ax, Ay, Az
      real(SP) :: B
      B = sqrt(Ax**2 + Ay**2 + Az**2) 
      return
   end function operator_vecelem_sp

   !> Internal elemental function that receives the component-wise input vector (DP)
   pure elemental function operator_vecelem_dp(Ax, Ay, Az) result(B)
      !$omp declare simd(operator_vecelem_dp)
      implicit none
      real(DP), intent(in) :: Ax, Ay, Az
      real(DP) :: B
      B = sqrt(Ax**2 + Ay**2 + Az**2) 
      return
   end function operator_vecelem_dp

   module procedure operator_mag_el_sp
      implicit none
      allocate(B(size(A,1)))
      B(:) = operator_vecelem_sp(A(:, 1), A(:, 2), A(:, 3))
      return
   end procedure operator_mag_el_sp

   module procedure operator_mag_el_dp
      implicit none
      allocate(B(size(A,1)))
      B(:) = operator_vecelem_dp(A(:, 1), A(:, 2), A(:, 3))
      return
   end procedure operator_mag_el_dp

end submodule operator_mag_implementation

