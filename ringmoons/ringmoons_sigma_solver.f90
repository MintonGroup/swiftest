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
subroutine ringmoons_sigma_solver(ring,GMP,dt)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_solver
      implicit none

! Arguments
      type(ringmoons_ring),intent(inout)  :: ring
      real(DP),intent(in)                 :: GMP,dt

! Internals

      real(DP),dimension(0:ring%N+1)      :: S,Snew,fac,fac2,TdSX
      integer(I4B)                        :: i,N
      real(DP)                            :: Gin,Gout

! Executable code

      N = ring%N
      S(0:ring%inside - 1) = 0.0_DP
      S(1:N) = ring%Gsigma(1:N) * ring%X(1:N)
      S(N+1) = 0.0_DP

      fac(:)  = 12 * dt / (ring%deltaX)**2  / ring%X2(:)
      fac2(:) = fac * 2._DP / (3 * PI * sqrt(GMP))
      do concurrent (i = 1:N) 
         Snew(i) = S(i)     + fac(i) * (ring%nu(i + 1) * S(i + 1) - 2 * ring%nu(i) * S(i) + ring%nu(i - 1) * S(i - 1))
         Snew(i) = Snew(i) - fac2(i) * (ring%Torque(i + 1) - 2 * ring%Torque(i) + ring%Torque(i - 1))
         Snew(i) = max(Snew(i),0.0_DP) 
      end do

      ring%Gsigma(1:N) = Snew(1:N) / ring%X(1:N)
      ring%Gm(1:N) = ring%Gsigma(1:N) * ring%deltaA(1:N)
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
