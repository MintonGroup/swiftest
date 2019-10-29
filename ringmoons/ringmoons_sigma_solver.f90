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
subroutine ringmoons_sigma_solver(ring,dt)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_solver
      implicit none

! Arguments
      type(ringmoons_ring),intent(inout)  :: ring
      real(DP),intent(in)                 :: dt

! Internals

      real(DP),dimension(0:ring%N+1)      :: S
      integer(I4B)                        :: i
      real(DP)                            :: fac, Snew

! Executable code

      S(0:ring%inside - 1) = 0.0_DP
      S(ring%N+1) = 0.0_DP
      S(1:ring%N) = ring%Gsigma(:) * ring%X(:)

      fac = 12 * dt / ring%deltaX**2 

      !!$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
      !!$OMP SHARED(ring,S,fac)
      do concurrent(i = 1:ring%N)
         Snew = S(i) + fac / (ring%X2(i)) * (ring%nu(i + 1) * S(i + 1) - 2 * ring%nu(i) * S(i) + ring%nu(i - 1) * S(i - 1))
         ring%Gsigma(i) = Snew / ring%X(i)
         ring%Gm(i) = ring%Gsigma(i) * ring%deltaA(i)
      end do
      !!$OMP END PARALLEL DO

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
