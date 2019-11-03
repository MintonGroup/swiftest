!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_sigma_torque
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Changes the surface mass density in response to satellite torques. Uses method of Salmon & Canup (2019)
!                An improved version would calculate the full stress tensor and derive an effective viscosity
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_sigma_torque(swifter_pl1P,ring,seeds)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_sigma_torque(swifter_pl1P,ring,dt)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_torque
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(inout)    :: ring
   real(DP), intent(in)                   :: dt

! Internals
   integer(I4B)                           :: i
   real(DP),dimension(0:ring%N+1)         :: dGsigma
   real(DP)                               :: fac,dGsig

! Executable code

   dGsigma(:) = 0.0_DP
   fac = -2 * dt / (sqrt(swifter_pl1P%mass) * ring%deltaX)
   do concurrent (i = 1:ring%N)
      if (abs(ring%Torque(i)) > 1e-25_DP) dGsigma(i) = fac * ring%Torque(i) / ring%deltaA(i) 
   end do
   do i = 1, ring%N
      if (dGsig > 0.0_DP) then
         dGsig = min(dGsigma(i),ring%Gsigma(i - 1))
      else
         dGsig = -min(-dGsigma(i),ring%Gsigma(i))
      end if
      ring%Gsigma(i - 1) = ring%Gsigma(i - 1) - dGsig
      ring%Gsigma(i)     = ring%Gsigma(i)     + dGsig
   end do
   do concurrent (i = 1:ring%N)
      ring%Gm(i) = ring%Gsigma(i) * ring%deltaA(i)
   end do
   return
   
end subroutine ringmoons_sigma_torque
