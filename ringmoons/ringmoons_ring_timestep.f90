!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_timestep
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the maximum stable timestep for the surface mass density evolution
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
!  Invocation  : CALL ringmoons_ring_timestep(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_ring_timestep(swifter_pl1P,ring,dtin) result(dtout)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_timestep
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(in)       :: ring
      real(DP), intent(in)                   :: dtin
      real(DP)                               :: dtout

! Internals
      integer(I4B)                           :: i
      real(DP)                               :: sig_max
      real(DP)                               :: torque_term
      

! Executable code

      ! Start with viscous stability
      dtout = dtin

      sig_max = 16 * (12._DP / (ring%deltaX)**2) * maxval(ring%nu(:) / ring%X2(:),ring%Gm(:) >  N_DISK_FACTOR * ring%Gm_pdisk )
      
      if (sig_max > 0.0_DP) then
         dtout = min(dtin,(sig_max)**(-1))
      end if
      

      return
end function ringmoons_ring_timestep
