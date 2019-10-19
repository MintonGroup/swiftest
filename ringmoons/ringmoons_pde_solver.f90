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
!    Terminal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Terminal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_pde_solver(dt,rm,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_pde_solver(dtin,rm,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_pde_solver
      IMPLICIT NONE

! Arguments
      real(DP),intent(in) :: dtin
      TYPE(ringmoons_parameter),INTENT(IN) :: rm
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      real(DP) :: dtstab, dt,fac1
      real(DP),dimension(rm%N) :: S,Snew
      integer(I4B) :: i,nloops,loop

! Executable code
      S(:) = ring%sigma(:) * ring%X(:) 
      dtstab = 0.5_DP * maxval(ring%X) * rm%deltaX**2 / (12 * minval(ring%nu))
      nloops = ceiling(dtin / dtstab)
      dt = dtin / nloops
      
      fac = 12 * dt / rm%deltaX**2
      do loop = 1,nloops  
          
         !$OMP PARALLEL DO DEFAULT(PRIVATE) &
         !$OMP SHARED(rm,ring,Snew,S,fac)
         do i = 2,rm%N - 1
            Snew(i) = S(i) + fac / (ring%X(i)**2) * (ring%nu(i) * (S(i + 1) - 2 * S(i) + S(i - 1)) &
                                                    + 0.5 * (S(i + 1) - S(i - 1)) * (ring%nu(i + 1) - ring%nu(i - 1)) &
                                                    + S(i) * (ring%nu(i + 1) - 2 * ring%nu(i) - ring%nu(i - 1)))
         end do
         !$OMP END PARALLEL DO
         S(:) = Snew(:)
         
      end do 

      RETURN

END SUBROUTINE ringmoons_pde_solver
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton  
!
!  Revision Control System (RCS) Information
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
