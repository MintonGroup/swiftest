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
      real(DP) :: rlo,rhi,GMP, RP,rhoP, Mratio, Rratio, Mratiosqrt,MratioHill,rfac
      real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew, Lorig,Mring,dMtot,Lnow
      real(DP),dimension(size(seeds%a)) :: afac
      real(DP),dimension(ring%N)        :: Gmtmp
      real(DP) :: Lp0,Ls0,Lp1,Ls1,Lr0,Lr1
      
! Executable code
      ! Save original mass and radius to use later
      Lr0 = sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
      Ls0 = sum(pack(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active(:)))
      Lp0 = swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
      Lorig = Lr0 + Ls0 + Lp0
      
      ring%inside = ringmoons_ring_bin_finder(ring,swifter_pl1P%radius)
      GMP = ring%GMPi + ring%dGMp
      dMtot = sum(ring%Gm(0:ring%inside))
      if (dMtot / GMP < epsilon(1._DP)) return
            
      !Add ring mass to planet
      ring%dGMP = ring%dGMP + dMtot
      swifter_pl1P%mass = ring%GMPi + ring%dGMP 
      swifter_pl1P%radius = ring%RPi * (swifter_pl1P%mass / ring%GMPi)**(1.0_DP / 3.0_DP)

      ring%Gm(0:ring%inside) = 0.0_DP
      ring%Gsigma(0:ring%inside) = 0.0_DP

      ! update body-dependent parameters as needed
      rfac = 1._DP - dMtot / GMP
      afac(:) = 1._DP - dMtot / (GMP + seeds%Gm(:))
      seeds%a(:) = seeds%a(:) * afac(:)
      ring%r_F = ring%r_F * rfac
      ring%r_I = ring%r_I * rfac

      ! Save the mass so that we can correct for the change in geometry
      Gmtmp(:) = ring%Gm
      call ringmoons_ring_construct(swifter_pl1P,ring,seeds)
      ring%Gm(:) = Gmtmp(:)
      ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)

      !write(*,*) 'after planet_accrete'
      Lr1 = sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
      Ls1 = sum(pack(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active(:)))
      Lp1 = swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
      Lnow = Lr1 + Ls1 + Lp1 

      swifter_pl1P%rot(3) = (Lp1 - Lnow + Lorig) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2)

      Lp1 = swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
      ring%dLP = Lp1 - ring%LPi

      return
         


      

end subroutine ringmoons_planet_accrete
