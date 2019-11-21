!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_deallocate
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
!  Invocation  : CALL ringmoons_deallocate(ring,seeds)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
sUBROUTINE ringmoons_deallocate(ring,seeds)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_deallocate
      implicit none

! Arguments
      type(ringmoons_ring),intent(inout) :: ring
      type(ringmoons_seeds),intent(inout) :: seeds

! Internals

! Executable code
      deallocate(ring%r)
      deallocate(ring%X)
      deallocate(ring%X2)
      deallocate(ring%deltaA)
      deallocate(ring%Gm)
      deallocate(ring%Gsigma)
      deallocate(ring%nu)
      deallocate(ring%Iz)
      deallocate(ring%w)
      deallocate(ring%Torque)
      deallocate(ring%r_hstar)

      deallocate(seeds%active)
      deallocate(seeds%a)
      deallocate(seeds%Gm)
      deallocate(seeds%Rhill)
      deallocate(seeds%rbin)
      deallocate(seeds%fz_bin_inner)
      deallocate(seeds%fz_bin_outer)
      deallocate(seeds%Torque)
      deallocate(seeds%Ttide)
      

      return

end subroutine ringmoons_deallocate
