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
         subroutine ringmoons_io_init_ring(swifter_pl1P,ring,seeds)
            use module_parameters
            use module_ringmoons
            implicit none
            type(swifter_pl),pointer :: swifter_pl1P
            type(ringmoons_ring),intent(inout) :: ring
            type(ringmoons_seeds),intent(inout)  :: seeds
         end subroutine ringmoons_io_init_ring
      end interface


      interface
         subroutine ringmoons_io_write_frame(t, ring, seeds, ring_outfile, out_stat)
         use module_parameters
         use module_ringmoons
         implicit none
         real(DP), intent(in)            :: t
         type(ringmoons_ring),intent(in) :: ring
         type(ringmoons_seeds),intent(in) :: seeds
         character(*), intent(in)  :: ring_outfile, out_stat      
         end subroutine ringmoons_io_write_frame
      end interface

      interface
         subroutine ringmoons_step(swifter_pl1P,ring,seeds,dtin,lfirst)
            use module_parameters
            use module_swifter
            use module_ringmoons
            implicit none
            type(swifter_pl), pointer                        :: swifter_pl1P
            real(DP), intent(in)                             :: dtin
            type(ringmoons_ring),intent(inout)               :: ring
            type(ringmoons_seeds),intent(inout)              :: seeds
            logical(LGT), intent(inout)                      :: lfirst
         end subroutine ringmoons_step 
      end interface

      interface
         subroutine ringmoons_sigma_solver(ring,dt)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         real(DP),intent(in) :: dt
         end subroutine ringmoons_sigma_solver
      end interface

      interface
         subroutine ringmoons_allocate(ring,seeds)
         use module_parameters
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         end subroutine ringmoons_allocate
      end interface

      interface
         subroutine ringmoons_deallocate(ring,seeds)
         use module_parameters
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         end subroutine ringmoons_deallocate
      end interface

      interface
         subroutine ringmoons_viscosity(ring)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_viscosity
      end interface

      interface
         function ringmoons_transition_function(y) result(kappa)
         use module_parameters
         implicit none
         real(DP),intent(in) :: y
         real(DP) :: kappa
         end function ringmoons_transition_function
      end interface

      interface
         subroutine ringmoons_planet_accrete(swifter_pl1P,ring,seeds)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         end subroutine ringmoons_planet_accrete
      end interface

      interface
         subroutine ringmoons_ring_construct(swifter_pl1P,ring)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         end subroutine ringmoons_ring_construct
      end interface

      interface
         subroutine ringmoons_seed_construct(swifter_pl1P,ring,seeds)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         end subroutine ringmoons_seed_construct
      end interface


      interface
         subroutine ringmoons_seed_grow(swifter_pl1P,ring,seeds,dt)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         real(DP),intent(in)                 :: dt
         end subroutine ringmoons_seed_grow
      end interface

      interface
         elemental function ringmoons_seed_dMdt(ring,GMP,Gsigma,Gmseed,a) result(Gmdot)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(ringmoons_ring), intent(in)       :: ring
         real(DP), intent(in)                   :: GMP,Gsigma,Gmseed,a
         real(DP)                               :: Gmdot
         end function ringmoons_seed_dMdt
      end interface




end module module_ringmoons_interfaces
