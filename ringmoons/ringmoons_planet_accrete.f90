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
    subroutine ringmoons_planet_accrete(swifter_pl1P,ring,seeds)

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

! Internals
      integer(I4B) :: i,j,iin
      real(DP) :: rlo,rhi,GMP, RP,rhoP, Mratio, Rratio, Mratiosqrt,MratioHill
      real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew
      real(DP),dimension(1:ring%N) :: rhill
     
! Executable code
      ! Save original mass and radius to use later
      Mratio = swifter_pl1P%mass
      Rratio = swifter_pl1P%radius
      do 
         iin = ring%inside
         do i = 1,iin 

           
            !Conserve angular momentum 
            Ltot = ring%Gm(i) * ring%Iz(i) * ring%w(i) + &
                   swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%radius**2
            
            !Add ring mass to planet
            Mnew = swifter_pl1P%mass + ring%Gm(i)
            Rnew = swifter_pl1P%radius * (Mnew / swifter_pl1P%mass)**(1.0_DP / 3.0_DP)
            swifter_pl1P%mass = Mnew
            swifter_pl1P%radius = Rnew

            swifter_pl1P%rot(3) = Ltot / (swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%radius**2)
            
            ring%Gm(i) = 0.0_DP
            ring%Gsigma(i) = 0.0_DP
         end do
         if (ring%rinner(iin) > swifter_pl1P%radius) exit !Find out if we need to update the inside bin
         ring%inside = ring%inside + 1
      end do
      Mratio = swifter_pl1P%mass / Mratio
      Rratio = swifter_pl1P%radius / Rratio
      Mratiosqrt = sqrt(Mratio)
      MratioHill = Mratio**(-1._DP / 3._DP)
      ! update body-dependent parameters as needed
      !$OMP PARALLEL WORKSHARE SHARED(Mratiosqrt,MratioHill)
      ring%w(:) = ring%w(:) * Mratiosqrt
      ring%r_hstar(:) = ring%r_hstar(:) * MRatioHill
      !$OMP END PARALLEL WORKSHARE
      ring%FRL = ring%FRL * Rratio
      ring%RRL = ring%RRL * Rratio 


      ! Adjust bin locations of RLs as necessary
      iin = ring%iRRL 
      do i = iin,ring%N
         if (ring%RRL < ring%router(i)) exit
         ring%iRRL = i
      end do

      iin = ring%iFRL 
      do i = iin,ring%N
         if (ring%FRL < ring%router(i)) exit
         ring%iFRL = i
      end do

      !$OMP PARALLEL WORKSHARE SHARED(MratioHill)
      seeds%Rhill(:) = seeds%Rhill(:) *  MratioHill
      !$OMP END PARALLEL WORKSHARE

      
      return
         


      

end subroutine ringmoons_planet_accrete
