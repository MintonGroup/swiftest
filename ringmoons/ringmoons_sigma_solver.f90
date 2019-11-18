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

      real(DP),dimension(0:ring%N+1)      :: S,Snew,Sn1,Sn2,fac,fac2,L,dM1,dM2
      integer(I4B)                        :: i,N,j

! Executable code

      N = ring%N
      where(ring%Gsigma(:) * ring%X(:) > TINY(1._DP) )
         S(:) = ring%Gsigma(:) * ring%X(:)
      elsewhere
         S(:) = 0._DP
      end where

      S(0) = 0.0_DP
      S(N+1) = 0.0_DP
      ring%Torque(0) = 0.0_DP
      ring%Torque(N+1) = 0.0_DP

      fac(:)  = 12 * dt / (ring%deltaX)**2  / ring%X2(:)

      Sn1(1:N) = ring%nu(2:N+1) * S(2:N+1) - 2 * ring%nu(1:N) * S(1:N) + ring%nu(0:N-1) * S(0:N-1)

      Sn2(1:N) = ring%Torque(2:N+1) - ring%Torque(0:N+1)

      Sn2(1:N) = Sn2(1:N) * (1._DP / (3 * PI * sqrt(GMP)))

      Snew(1:N) = S(1:N) + fac(1:N) * (Sn1(1:N) - Sn2(1:N))

      ring%Gsigma(1:N) = Snew(1:N) / ring%X(1:N)
      ring%Gm(1:N) = ring%Gsigma(1:N) * ring%deltaA(1:N)

      ! Prevent any bins from having negative mass by shifting mass upward from interior bins  
      i = 1
      do while (any(ring%Gm(1:N) < 0.0_DP))
         !write(*,*) i,'Negative mass!'
         i = i + 1
         where(ring%Gm(:) < 0.0_DP)
            dM1(:) = ring%Gm(:)
         elsewhere
            dM1(:) = 0.0_DP
         end where
         L(:) = ring%Iz(:) * ring%w(:)
         dM2(:) = dM1(:) * (L(:) - cshift(L(:),1)) / (cshift(L(:),1) - cshift(L(:),2)) 
        ! Make sure we conserve both mass and angular momentum
         ring%Gm(1:N) = ring%Gm(1:N) - dM1(1:N) + cshift(dM1(1:N),1) + &
            cshift(dM2(1:N),1)  - cshift(dM2(1:N),2)
         ring%Gsigma(1:N) = ring%Gm(1:N) / ring%deltaA(1:N)
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
