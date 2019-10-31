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
            use module_swifter
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
         subroutine ringmoons_step(t,swifter_pl1P,ring,seeds,dtin,lfirst,Merror,Lerror)
            use module_parameters
            use module_swifter
            use module_ringmoons
            implicit none
            type(swifter_pl), pointer                        :: swifter_pl1P
            real(DP), intent(in)                             :: t,dtin
            type(ringmoons_ring),intent(inout)               :: ring
            type(ringmoons_seeds),intent(inout)              :: seeds
            logical(LGT), intent(inout)                      :: lfirst
            real(DP), intent(out)                            :: Merror,Lerror
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
         pure function ringmoons_transition_function(y) result(kappa)
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
         subroutine ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dt)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         real(DP),intent(in)                 :: dt
         end subroutine ringmoons_seed_evolve
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


      interface
         elemental function ringmoons_seed_dadt(GMP,Gmseed,a,Torque) result(adot)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         real(DP), intent(in)                   :: GMP,Gmseed,a,Torque
         real(DP)                               :: adot
         end function ringmoons_seed_dadt
      end interface

      interface
         elemental function ringmoons_ring_bin_finder(ring,r) result(bin)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(ringmoons_ring), intent(in)       :: ring
         real(DP), intent(in)                   :: r
         integer(I4B)                           :: bin
         end function ringmoons_ring_bin_finder
      end interface

      interface
         function ringmoons_timestep(swifter_pl1P,ring,seeds,dtin) result(dtout)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer               :: swifter_pl1P
         type(ringmoons_ring),intent(in)        :: ring
         type(ringmoons_seeds),intent(in)       :: seeds
         real(DP),intent(in)                    :: dtin
         real(DP)                               :: dtout
         end function ringmoons_timestep
      end interface


      interface
         recursive function ringmoons_laplace_coefficient(alpha,j,s,n) result(ans)
         use module_parameters
         implicit none
         real(DP),intent(in)  :: alpha,s
         integer, intent(in)  :: j,n
         real(DP)             :: ans         
         end function ringmoons_laplace_coefficient
      end interface

      interface
         function ringmoons_lindblad_torque(swifter_pl1P,ring,Gm,a,e,inc) result(Torque)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer                  :: swifter_pl1P
         type(ringmoons_ring), intent(inout)       :: ring
         real(DP),intent(in)                       :: Gm,a,e,inc
         real(DP)                                  :: Torque         
         end function ringmoons_lindblad_torque
      end interface

      interface
         function ringmoons_tidal_torque(swifter_pl1P,Gm,n,a,e,inc) result(Torque)
         use module_parameters
         use module_swifter
         implicit none
         type(swifter_pl),pointer                  :: swifter_pl1P
         real(DP),intent(in)                       :: Gm,n,a,e,inc
         real(DP)                                  :: Torque         
         end function ringmoons_tidal_torque
      end interface


      interface
         subroutine ringmoons_calc_torques(swifter_pl1P,ring,seeds)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer :: swifter_pl1P
         type(ringmoons_ring),intent(inout) :: ring
         type(ringmoons_seeds),intent(inout) :: seeds
         end subroutine ringmoons_calc_torques
      end interface

      interface
         subroutine ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,Gm)
         use module_parameters
         use module_swifter
         use module_ringmoons
         implicit none
         type(swifter_pl),pointer               :: swifter_pl1P
         type(ringmoons_ring), intent(inout)    :: ring
         type(ringmoons_seeds), intent(inout)   :: seeds
         real(DP), intent(in)                   :: a, Gm
         end subroutine
      end interface

end module module_ringmoons_interfaces
