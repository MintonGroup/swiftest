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
      real(DP),dimension(seeds%N) :: afac
      real(DP),dimension(0:ring%N+1)        :: Gmtmp,Lring_orig,Lring_now,dL
      real(DP) :: Lp0,Ls0,Lp1,Ls1,Lr0,Lr1
      
! Executable code
      ring%inside = ringmoons_ring_bin_finder(ring,swifter_pl1P%radius)
      GMP = ring%GMPi + ring%dGMp
      dMtot = sum(ring%Gm(0:ring%inside))
            
      !Add ring mass and angular momentum to planet
      Lring_orig(:) = ring%Gm(:) * ring%Iz(:) * ring%w(:) 
      Lring = sum(Lring_orig(0:ring%inside))
      ring%dGMP = ring%dGMP + dMtot
      ring%dLP = ring%dLP + Lring
      swifter_pl1P%mass = ring%GMPi + ring%dGMP 
      swifter_pl1P%radius = ring%RPi * (1._DP + ring%dGMP / ring%GMPi)**(1.0_DP / 3.0_DP)
      swifter_pl1P%rot(3) = (ring%LPI + ring%dLP) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2)

      ring%Gm(0:ring%inside) = 0.0_DP
      ring%Gsigma(0:ring%inside) = 0.0_DP

      ! update body-dependent parameters as needed
      rfac = 1._DP - dMtot / GMP
      afac(1:seeds%N) = 1._DP - dMtot / (GMP + seeds%Gm(1:seeds%N))
      seeds%a(1:seeds%N) = seeds%a(1:seeds%N) * afac(1:seeds%N)
      ring%r_F = ring%r_F * rfac
      ring%r_I = ring%r_I * rfac

      ! Save the mass so that we can correct for the change in geometry
      Gmtmp(:) = ring%Gm
      call ringmoons_ring_construct(swifter_pl1P,ring,seeds)
      ring%Gm(:) = Gmtmp(:)
      ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
      
      ! Any difference in angular momentum in each ring bin will result in a torque in that bin
      Lring_now(:) = ring%Gm(:) * ring%Iz(:) * ring%w(:) 
      dL(0:ring%inside) = 0.0_DP
      dL(ring%N+1) = 0.0_DP
      dL(ring%inside+1:ring%N) = (Lring_now(ring%inside+1:ring%N) - Lring_orig(ring%inside+1:ring%N)) / dt 
      where(abs(dL(:)) > epsilon(1._DP))
         ring%Torque(:) = ring%Torque(:) - dL(:)
      end where

      return

end subroutine ringmoons_planet_accrete
