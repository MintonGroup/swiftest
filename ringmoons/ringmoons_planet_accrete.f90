!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_planet_accrete
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : solves
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
!  Invocation  : call ringmoons_planet_accrete(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
    subroutine ringmoons_planet_accrete(swifter_pl1P,ring)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_planet_accrete
      implicit NONE

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring),intent(inout) :: ring

! Internals
      integer(I4B) :: i,iin
      real(DP) :: rlo,rhi,GM_Planet
      real(DP) :: Lplanet, Lring, Ltot
      real(DP),dimension(1:ring%N) :: rhill
     
      do 
         iin = ring%inside
         do i = 1,iin 

            !Add ring mass to planet
            swifter_pl1P%radius = swifter_pl1P%radius * ((swifter_pl1P%mass + ring%Gm(i)) / swifter_pl1P%mass)**(1.0_DP / 3.0_DP)
            swifter_pl1P%mass = swifter_pl1P%mass + ring%Gm(i)
            ring%Gm(i) = 0.0_DP
            ring%Gsigma(i) = 0.0_DP
           
            !Conserve angular momentum 
            Ltot = ring%Iz(i) * ring%w(i) + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3)
            swifter_pl1P%rot(3) = Ltot / swifter_pl1P%Ip(3)
            ring%Iz(i) = 0.0_DP
         end do
         if (ring%rinner(iin) > swifter_pl1P%radius) exit !Find out if we need to update the inside bin
         ring%inside = ring%inside + 1
      end do
      GM_Planet = swifter_pl1P%mass
      ! update body-dependent parameters as needed
      ring%w = sqrt(GM_Planet / ring%r**3)
      rhill = ring%r * (2 * ring%Gm_pdisk /(3._DP * GM_Planet))**(1._DP/3._DP) ! See Salmon et al. 2010 for this
      ring%r_hstar = rhill / (2 * ring%r_pdisk)   

! Executable code

      

end subroutine ringmoons_planet_accrete
