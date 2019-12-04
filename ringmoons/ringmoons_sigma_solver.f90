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
subroutine ringmoons_sigma_solver(ring,GMP,dt,stepfail)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_sigma_solver
      implicit none

! Arguments
      type(ringmoons_ring),intent(inout)  :: ring
      real(DP),intent(in)                 :: GMP,dt
      logical(lgt),intent(out)            :: stepfail

! Internals

      real(DP),dimension(0:ring%N+1)      :: S,Snew,Sn1,Sn2,fac,artnu,L,dM1,dM2
      integer(I4B)                        :: i,N,j,loop

! Executable code
      stepfail = .false.

      N = ring%N

      S(0) = 0.0_DP
      S(1:N) = ring%Gsigma(1:N) * ring%X(1:N)
      S(N+1) = 0.0_DP
      ring%Torque(0) = 0.0_DP
      ring%Torque(N+1) = 0.0_DP
      ring%Gm(0) = 0.0_DP
      ring%Gm(N+1) = 0.0_DP
      ring%Gsigma(0) = 0.0_DP
      ring%Gsigma(N+1) = 0.0_DP


      if (any(S(:) < 0.0_DP)) then
         write(*,*) 'Negative mass already?'
         read(*,*)
      end if


      fac(:)  = 12 * dt / (ring%deltaX)**2  / ring%X2(:)

      Sn1(1:N) = ring%nu(2:N+1) * S(2:N+1) - 2 * ring%nu(1:N) * S(1:N) + ring%nu(0:N-1) * S(0:N-1)

      Sn2(1:N) = ring%Torque(2:N+1) - ring%Torque(0:N-1)

      Sn2(1:N) = Sn2(1:N) * (1._DP / (3 * PI * sqrt(GMP)))

      Snew(1:N) = S(1:N) + fac(1:N) * (Sn1(1:N) - Sn2(1:N))

      ! Prevent any bins from having negative mass by diffusing mass with an artificial viscosity
      loop = 1
      do while (any(Snew(1:N) < -epsilon(1._DP) * maxval(S(1:N))))
         Snew(0) = 0.0_DP
         Snew(N+1) = 0.0_DP
         artnu(:) = 0.0_DP
         where (Snew(1:N) < 0._DP)
            artnu(1:N) = 1._DP / (16 * fac(1:N))
            artnu(0:N-1) = 1._DP / (16 * fac(0:N-1))
            artnu(2:N+1) = 1._DP / (16 * fac(2:N+1))
         end where
         S(1:N) = Snew(1:N) 
         Sn1(1:N) = artnu(2:N+1) * S(2:N+1) - 2 * artnu(1:N) * S(1:N) + artnu(0:N-1) * S(0:N-1)
         Snew(1:N) = S(1:N) + fac(1:N) * Sn1(1:N) 
         loop = loop + 1
         if (loop > 100) then
            stepfail = .true.
            exit
         end if
      end do

      ring%Gsigma(1:N) = max(0.0_DP,Snew(1:N) / ring%X(1:N))
      ring%Gm(1:N) = ring%Gsigma(1:N) * ring%deltaA(1:N)
      ring%Torque(:) = 0.0_DP

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
