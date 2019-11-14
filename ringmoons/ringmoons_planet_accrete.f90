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
    subroutine ringmoons_planet_accrete(swifter_pl1P,ring,seeds,dt)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_planet_accrete
      implicit NONE

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring),intent(inout) :: ring
      type(ringmoons_seeds),intent(inout) :: seeds
      real(DP),intent(in)                 :: dt

! Internals
      integer(I4B) :: i,j,iin
      real(DP) :: rlo,rhi,GMP, RP,rhoP, Mratio, Rratio, Mratiosqrt,MratioHill,deltaMp
      real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew, Lorig,Mring,dMtot
      real(DP),dimension(1:ring%N) :: rhill
     
! Executable code
      ! Save original mass and radius to use later
      Mratio = swifter_pl1P%mass
      Rratio = swifter_pl1P%radius
      dMtot = 0.0_DP
      !Lplanet = swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * (swifter_pl1P%radius)**2
      do 
         iin = ring%inside
         do i = iin,1,-1
            if (ring%Gm(i) / swifter_pl1P%mass < epsilon(1._DP)) exit
            !Conserve angular momentum 
            Lring = ring%Gm(i) * ring%Iz(i) * ring%w(i)
            ring%dLP = ring%dLP + Lring
            
            !Add ring mass to planet
            ring%dGMP = ring%dGMP + ring%Gm(i)
            dMtot = dMtot + ring%Gm(i)
            swifter_pl1P%mass = ring%GMPi + ring%dGMP 
            swifter_pl1P%radius = ring%RPi * (swifter_pl1P%mass / ring%GMPi)**(1.0_DP / 3.0_DP)
             
            swifter_pl1P%rot(3) = (ring%LPi + ring%dLP) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2)

            ring%Gm(i) = 0.0_DP
            ring%Gsigma(i) = 0.0_DP

         end do
         if (ring%rinner(iin) > swifter_pl1P%radius) exit !Find out if we need to update the inside bin
         ring%inside = ring%inside + 1
      end do
      Mratio = (ring%GMPi + ring%dGMP) / (ring%GMPi + ring%dGMP - dMtot)
      Rratio = Mratio**(1.0_DP / 3.0_DP)
      Mratiosqrt = 0.5_DP * dMtot / (ring%GMPi + ring%dGMP - dMtot)!sqrt(Mratio)
      MratioHill = Mratio**(-1._DP / 3._DP)
      ! update body-dependent parameters as needed
      Mring = sum(ring%Gm(:))
      if ((Mratiosqrt - 1._DP) > tiny(1._DP)) then
         ring%Torque(:) = ring%Torque(:) - ring%Gm(:) * ring%Iz(:) * ring%w(:) * (Mratiosqrt) / dt
         ring%w(:) = ring%w(:) * (Mratiosqrt + 1.0_DP)
         ring%r_hstar(:) = ring%r_hstar(:) * MRatioHill
      end if
      ring%FRL = ring%FRL * Rratio
      ring%RRL = ring%RRL * Rratio 


      ! Adjust bin locations of RLs as necessary
      ring%iFRL = ringmoons_ring_bin_finder(ring,ring%FRL)
      ring%iRRL = ringmoons_ring_bin_finder(ring,ring%RRL)

      
      !do (i = 1,seeds%N)
         !if (seeds%active(i)) then
      seeds%Rhill(:) = seeds%Rhill(:) *  MratioHill
      seeds%a(:) = seeds%a(:) * (ring%GMPi + ring%dGMP - ring%dGMP + seeds%Gm(:)) / (ring%GMPi + ring%dGMP + seeds%Gm(:))
         !end if
      !end do

      
      return
         


      

end subroutine ringmoons_planet_accrete
