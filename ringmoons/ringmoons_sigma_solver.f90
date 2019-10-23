!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_sigma_solver
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
!  Invocation  : CALL ringmoons_sigma_solver(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
subroutine ringmoons_sigma_solver(swifter_pl1P,ring,dtin)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_solver
      implicit none

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring),intent(inout) :: ring
      real(DP),intent(in) :: dtin

! Internals
      real(DP) :: dtstab,dtleft,dt,fac, GM_Planet,Snew
      real(DP),dimension(0:ring%N+1) :: S
      integer(I4B) :: i,loop
      real(DP),parameter :: GMURN = 5.74811598e+21

! Executable code
      dtleft = dtin
      !TESTING
         call ringmoons_viscosity(swifter_pl1P%mass,ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         write(*,*) dtstab,ceiling(dtin/dtstab),(sum(ring%Gm) + (swifter_pl1P%mass - GMURN)) / GU
      !^^^^^^^^  
      do loop = 1, LOOPMAX
         call ringmoons_viscosity(swifter_pl1P%mass,ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         dt = min(dtleft,dtstab)
         S(0:ring%inside - 1) = 0.0_DP
         S(ring%N+1) = 0.0_DP
         S(1:ring%N) = ring%Gsigma(:) / GU * ring%X(:)

         fac = 12 * dt / ring%deltaX**2 
         !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
         !$OMP SHARED(ring,S,fac,GU)
         do i = 1,ring%N
            Snew = S(i) + fac / (ring%X2(i)) * (ring%nu(i + 1) * S(i + 1) - 2 * ring%nu(i) * S(i) + ring%nu(i - 1) * S(i - 1))
            ring%Gsigma(i) = GU * Snew / ring%X(i)
         end do
         !$OMP END PARALLEL DO
         ring%Gm = ring%Gsigma * ring%deltaA
         ring%Iz = 0.5_DP * ring%Gm / GU * (ring%rinner**2 + ring%router**2)
         call ringmoons_planet_accrete(swifter_pl1P,ring)
         dtleft = dtleft - dt
         if (dtleft <= 0.0_DP) exit
      end do 

      return

end subroutine ringmoons_sigma_solver
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton  
!
!  Revision Control System (RCS) Inforingation
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
