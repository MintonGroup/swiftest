!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swifter
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in Swifter package
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
module module_ringmoons_interfaces

      implicit none


      interface
         subroutine ringmoons_io_init_ring(GM_Planet,R_Planet,ring)
            use module_parameters
            use module_ringmoons
            implicit none
            real(DP),intent(in)     :: GM_Planet,R_Planet 
            type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_io_init_ring
      end interface

      interface
         subroutine ringmoons_step(lfirst, t, rmin, npl, nplmax, symba_pl1p, j2rp2, j4rp4, eoffset, dt,ring)
            use module_parameters
            use module_symba
            use module_ringmoons
            implicit none
            logical(lgt), intent(inout)                      :: lfirst
            integer(I4B), intent(in)                         :: npl, nplmax
            real(DP), intent(in)                             :: t, rmin, j2rp2, j4rp4, dt
            real(DP), intent(inout)                          :: eoffset
            type(symba_pl), pointer                          :: symba_pl1p
            type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_step 
      end interface

      interface
         subroutine ringmoons_pde_solver(GM_Planet,R_Planet,dtin,ring)
         use module_parameters
         use module_ringmoons
         implicit none
         real(DP),intent(in) :: GM_Planet,R_Planet,dtin
         TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         end subroutine ringmoons_pde_solver
      end interface

      interface
         subroutine ringmoons_allocate(ring)
         use module_parameters
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_allocate
      end interface

      interface
         subroutine ringmoons_deallocate(ring)
         use module_parameters
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_deallocate
      end interface

      interface
         subroutine ringmoons_viscosity(GM_Planet,R_Planet,ring)
         use module_parameters
         use module_ringmoons
         implicit none
         real(DP),intent(in) :: GM_Planet,R_Planet
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_viscosity
      end interface

      interface
         subroutine ringmoons_viscosity_s2010(GM_Planet,R_Planet,ring)
         use module_parameters
         use module_ringmoons
         implicit none
         real(DP),intent(in) :: GM_Planet,R_Planet
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_viscosity_s2010
      end interface

      interface
         function ringmoons_transition_function(y) result(kappa)
         use module_parameters
         implicit none
         real(DP),intent(in) :: y
         real(DP) :: kappa
         end function ringmoons_transition_function
      end interface



end module module_ringmoons_interfaces
