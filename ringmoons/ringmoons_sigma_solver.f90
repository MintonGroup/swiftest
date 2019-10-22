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
SUBROUTINE ringmoons_sigma_solver(GM_Planet,R_Planet,dtin,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_solver
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: GM_Planet,R_Planet,dtin
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: dtstab,dtleft,dt,fac
      real(DP),dimension(0:ring%N+1) :: S,Snew
      integer(I4B) :: i,loop

! Executable code
      S(1:ring%N) = ring%Gsigma(:) / GU * ring%X(:)
      dtleft = dtin
      !TESTING
         call ringmoons_viscosity(GM_Planet,R_Planet,ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         write(*,*) dtstab,ceiling(dtin/dtstab),sum(ring%Gm) / GU
      !^^^^^^^^  
      do loop = 1, LOOPMAX
         call ringmoons_viscosity(GM_Planet,R_Planet,ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         dt = min(dtleft,dtstab)
         S(0) = 0.0_DP
         S(ring%N+1) = 0.0_DP
         Snew = 0.0_DP
         fac = 12 * dt / ring%deltaX**2 
         !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
         !$OMP SHARED(ring,Snew,S,fac,GU)
         do i = 1,ring%N
            Snew(i) = S(i) + fac / (ring%X2(i)) * (ring%nu(i + 1) * S(i + 1) - 2 * ring%nu(i) * S(i) + ring%nu(i - 1) * S(i - 1))
            ring%Gsigma(i) = GU * Snew(i) / ring%X(i)
         end do
         !$OMP END PARALLEL DO
         S(:) = Snew(:)
         ring%Gm = ring%Gsigma * ring%deltaA
         dtleft = dtleft - dt
         if (dtleft <= 0.0_DP) exit
      end do 

      return

END SUBROUTINE ringmoons_sigma_solver
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
