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
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_io_init_ring(swifter_pl1P,ring,seeds)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_io_init_ring
   implicit none

! Arguments
   type(swifter_pl),pointer            :: swifter_pl1P
   type(ringmoons_ring),intent(inout) :: ring
   type(ringmoons_seeds),intent(inout) :: seeds

! Internals
   character(STRMAX)                   :: ringfile
   integer(I4B),parameter              :: LUN = 22
   integer(I4B)                        :: i,m,inner_outer_sign,ioerr
   real(DP)                            :: beta


! Executable code
   ringfile='ring.in'
   open(unit=LUN,file=ringfile,status='old',iostat=ioerr)
   read(LUN,*) ring%N, seeds%N
   read(LUN,*) ring%r_I, ring%r_F
   read(LUN,*) ring%r_pdisk,ring%Gm_pdisk
   call ringmoons_allocate(ring,seeds)
   do i = 1,ring%N
      read(LUN,*,iostat=ioerr) ring%Gsigma(i)
      if (ioerr /= 0) then 
         write(*,*) 'File read error in ',trim(adjustl(ringfile))
         write(*,*) 'Reading ring element ', i,' of ',ring%N
         call util_exit(FAILURE)
      end if
   end do

   do i = 1,seeds%N
      read(LUN,*,iostat=ioerr) seeds%a(i), seeds%Gm(i)
      if (ioerr /= 0) then 
         write(*,*) 'File read error in ',trim(adjustl(ringfile))
         write(*,*) 'Reading seed ', i,' of ',seeds%N
         call util_exit(FAILURE)
      end if
      seeds%active(i) = .true.

   end do

! For performance reasons, we compute a table of Laplace coefficient terms the first time through and then interpolate 
   do m = 2, m_max
      do inner_outer_sign = -1,1,2
         beta =  (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(-inner_outer_sign * 2._DP / 3._DP)
         lapm(inner_outer_sign,m)  = m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) 
         dlapm(inner_outer_sign,m) = 0.5_DP * beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1) 
      end do
   end do


   return

end subroutine ringmoons_io_init_ring
