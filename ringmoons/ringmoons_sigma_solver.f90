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
      integer(I4B)                        :: i,N,j

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

      fac(:)  = 12 * dt / (ring%deltaX)**2  / ring%X2(:)

      Sn1(1:N) = ring%nu(2:N+1) * S(2:N+1) - 2 * ring%nu(1:N) * S(1:N) + ring%nu(0:N-1) * S(0:N-1)

      Sn2(1:N) = ring%Torque(2:N+1) - ring%Torque(0:N-1)

      Sn2(1:N) = Sn2(1:N) * (1._DP / (3 * PI * sqrt(GMP)))

      Snew(1:N) = S(1:N) + fac(1:N) * (Sn1(1:N) - Sn2(1:N))

      
      !if (any(Snew(1:N) < 0.0_DP)) then
      !   do i = 1,N
      !      if (Snew(i) < 0.0_DP) write(*,*) i-1,Snew(i-1)
      !      if (Snew(i) < 0.0_DP) write(*,*) i,Snew(i)
      !      if (Snew(i) < 0.0_DP) write(*,*) i+1,Snew(i+1)
      !   end do
      !   write(*,*)
      !end if

      ! Prevent any bins from having negative mass by diffusing mass with an artificial viscosity
      do while (any(Snew(1:N) < -tiny(1._DP)))
         where (Snew(1:N) < 0.0_DP)
            artnu(1:N) = 1._DP / (16 * fac(1:N))
         elsewhere
            artnu(1:N) = 0.0_DP
         end where
         S(1:N) = Snew(1:N) 
         Sn1(1:N) = artnu(2:N+1) * S(2:N+1) - 2 * artnu(1:N) * S(1:N) + artnu(0:N-1) * S(0:N-1)
         Snew(1:N) = S(1:N) + fac(1:N) * Sn1(1:N) 

         !do i = 1,N
         !   if (Snew(i) < 0.0_DP) write(*,*) i-1,Snew(i-1),artnu(i-1)
         !   if (Snew(i) < 0.0_DP) write(*,*) i,Snew(i),artnu(i)
         !   if (Snew(i) < 0.0_DP) write(*,*) i+1,Snew(i+1),artnu(i+1)
         !end do
         !read(*,*)
      end do

      ring%Gsigma(1:N) = max(0.0_DP,Snew(1:N) / ring%X(1:N))
      ring%Gm(1:N) = ring%Gsigma(1:N) * ring%deltaA(1:N)

      !do while (any(ring%Gm(1:ring%N) < 0.0_DP))
      !   where(ring%Gm(:) < 0.0_DP)
      !      dM1(:) = ring%Gm(:)
      !   elsewhere
      !      dM1(:) = 0.0_DP
      !   end where
      !   L(:) = ring%Iz(:) * ring%w(:)
      !   dM2(:) = dM1(:) * (L(:) - cshift(L(:),1)) / (cshift(L(:),1) - cshift(L(:),2)) 
      !  ! Make sure we conserve both mass and angular momentum
      !   ring%Gm(1:ring%N) = ring%Gm(1:ring%N) - dM1(1:ring%N) + cshift(dM1(1:ring%N),1) + &
      !      cshift(dM2(1:ring%N),1)  - cshift(dM2(1:ring%N),2)
      !   ring%Gsigma(1:ring%N) = ring%Gm(1:ring%N) / ring%deltaA(1:ring%N)
      !end do 

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
