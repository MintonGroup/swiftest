!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_io_init_ring
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in ring initial condition data
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
!  Invocation  : CALL ringmoons_io_init_ring()
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_io_init_ring(GM_Planet,R_Planet,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_io_init_ring
      IMPLICIT NONE

! Arguments
      real(DP),intent(in)                :: GM_Planet,R_Planet
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals
      character(STRMAX)                  :: ringfile
      integer(I4B),parameter             :: LUN = 22
      integer(I4B)                       :: i,ioerr
      real(DP)                           :: Xlo,Xhi

      ringfile='ring.in'
      open(unit=LUN,file=ringfile,status='old',iostat=ioerr)
      read(LUN,*) ring%N
      read(LUN,*) ring%r_I, ring%r_F
      read(LUN,*) ring%r_pdisk,ring%m_pdisk
      call ringmoons_allocate(ring)
      ring%deltaX = (2 * sqrt(ring%r_F) - 2 * sqrt(ring%r_I)) / ring%N
      do i = 1,ring%N
         read(LUN,*,iostat=ioerr) ring%sigma(i)
         if (ioerr /= 0) then 
            write(*,*) 'File read error in ',trim(adjustl(ringfile))
         end if
         Xlo = 2 * sqrt(ring%r_I) + ring%deltaX * (1._DP * i - 1._DP)
         Xhi = Xlo + ring%deltaX
         ring%deltar(i) = 0.5_DP * (Xhi**2 - Xlo**2)

         ring%X(i) = Xlo + 0.5_DP * ring%deltaX
         ring%r(i) = (0.5_DP * ring%X(i))**2
         ring%R_P(i) = ring%r(i) / R_Planet
         ring%deltaA(i) = 2 * PI * ring%deltar(i) * ring%r(i)
         ring%m(i) = ring%sigma(i) * ring%deltaA(i)
         ring%RR(i) = ring%r(i)**2 + 0.25_DP * ring%deltar(i)**2 
         ring%I(i) = ring%m(i) * ring%RR(i)
         ring%w(i) = sqrt(GM_Planet / ring%r(i)**3)
         ring%Torque_to_disk(i) = 0.0_DP
      end do
      call ringmoons_viscosity(GM_Planet,R_Planet,ring)

      
   

! Executable code


      RETURN

END SUBROUTINE ringmoons_io_init_ring
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
