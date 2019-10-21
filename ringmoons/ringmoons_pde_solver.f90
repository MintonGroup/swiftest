!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_pde_solver
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
!  Invocation  : CALL ringmoons_pde_solver(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_pde_solver(GM_Planet,R_Planet,dtin,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_pde_solver
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: GM_Planet,R_Planet,dtin
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: dtstab, dt,fac
      real(DP),dimension(ring%N) :: S,Snew
      integer(I4B) :: i,nloops,loop

! Executable code
      call ringmoons_viscosity(GM_Planet,R_Planet,ring)
      S(:) = ring%Gsigma(:) / GU * ring%X(:)
      dtstab = 0.5_DP * minval(ring%X) * ring%deltaX**2 / (12 * maxval(ring%nu))
      nloops = ceiling(dtin / dtstab)
      dt = dtin / nloops
      write(*,*) dtstab,nloops
      
      fac = 12 * dt / ring%deltaX**2 
      do loop = 1,nloops  
         Snew(1) = S(1) + fac / (ring%X(1)**2) * (ring%nu(1) * (S(2) - 2 * S(1)) &
                                                 + 0.5_DP * (S(2)) * (ring%nu(2)) &
                                                 + S(1) * (ring%nu(2) - 2 * ring%nu(1)))
         ring%Gsigma(1) = GU * Snew(1) / ring%X(1)
         !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
         !$OMP SHARED(ring,Snew,S,fac,GU)
         do i = 2,ring%N - 1
            Snew(i) = S(i) + fac / (ring%X(i)**2) * (ring%nu(i) * (S(i + 1) - 2 * S(i) + S(i - 1)) &
                                                    + 0.5_DP * (S(i + 1) - S(i - 1)) * (ring%nu(i + 1) - ring%nu(i - 1)) &
                                                    + S(i) * (ring%nu(i + 1) - 2 * ring%nu(i) + ring%nu(i - 1)))
            ring%Gsigma(i) = GU * Snew(i) / ring%X(i)
         end do
         !$OMP END PARALLEL DO
         i = ring%N
         Snew(i) = S(i) + fac / (ring%X(i)**2) * (ring%nu(i) * (-2 * S(i) + S(i - 1)) &
                                                    + 0.5_DP * ( -S(i - 1)) * ( -ring%nu(i - 1)) &
                                                    + S(i) * (- 2 * ring%nu(i) + ring%nu(i - 1)))
         S(:) = Snew(:)
         ring%m = ring%Gsigma * ring%deltaA
         call ringmoons_viscosity(GM_Planet,R_Planet,ring)
      end do 

      RETURN

END SUBROUTINE ringmoons_pde_solver
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
