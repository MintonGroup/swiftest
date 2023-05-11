!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module whm
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest
   implicit none
   public

   !> Swiftest central body particle class
   type, extends(swiftest_cb) :: whm_cb
   contains
   end type whm_cb


   !> WHM massive body particle class
   type, extends(swiftest_pl) :: whm_pl
      real(DP), dimension(:),   allocatable :: eta    !! Jacobi mass
      real(DP), dimension(:,:), allocatable :: xj     !! Jacobi position
      real(DP), dimension(:,:), allocatable :: vj     !! Jacobi velocity
      real(DP), dimension(:),   allocatable :: muj    !! Jacobi mu: GMcb * eta(i) / eta(i - 1) 
      real(DP), dimension(:),   allocatable :: ir3j    !! Third term of heliocentric acceleration
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as whm_util_setup_pl and whm_util_spill_pl
   contains
      procedure :: h2j         => whm_coord_h2j_pl           !! Convert position and velcoity vectors from heliocentric to Jacobi coordinates 
      procedure :: j2h         => whm_coord_j2h_pl           !! Convert position and velcoity vectors from Jacobi to helliocentric coordinates 
      procedure :: vh2vj       => whm_coord_vh2vj_pl         !! Convert velocity vectors from heliocentric to Jacobi coordinates 
      procedure :: drift       => whm_drift_pl               !! Loop through massive bodies and call Danby drift routine to jacobi coordinates
      procedure :: accel_gr    => whm_gr_kick_getacch_pl     !! Acceleration term arising from the post-Newtonian correction
      procedure :: gr_pos_kick => whm_gr_p4_pl               !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: accel       => whm_kick_getacch_pl        !! Compute heliocentric accelerations of massive bodies
      procedure :: kick        => whm_kick_vh_pl             !! Kick heliocentric velocities of massive bodies
      procedure :: append      => whm_util_append_pl         !! Appends elements from one structure to another
      procedure :: dealloc     => whm_util_dealloc_pl        !! Deallocates all allocatable arrays
      procedure :: fill        => whm_util_fill_pl           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize      => whm_util_resize_pl         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: set_ir3     => whm_util_set_ir3j          !! Sets both the heliocentric and jacobi inverse radius terms (1/rj**3 and 1/rh**3)
      procedure :: set_mu      => whm_util_set_mu_eta_pl     !! Sets the Jacobi mass value for all massive bodies.
      procedure :: sort        => whm_util_sort_pl           !! Sort a WHM massive body object in-place. 
      procedure :: rearrange   => whm_util_sort_rearrange_pl !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill       => whm_util_spill_pl          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      procedure :: setup       => whm_util_setup_pl          !! Constructor method - Allocates space for the input number of bodiess
      procedure :: step        => whm_step_pl                !! Steps the body forward one stepsize
      final     ::                whm_final_pl               !! Finalizes the WHM massive body object - deallocates all allocatables
#ifdef COARRAY
      procedure :: coclone      => whm_coarray_coclone_pl       !! Clones the image 1 body object to all other images in the coarray structure.
#endif 
   end type whm_pl


   !! WHM test particle class
   type, extends(swiftest_tp) :: whm_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as whm_util_spill_tp
   contains
      procedure :: accel_gr    => whm_gr_kick_getacch_tp !! Acceleration term arising from the post-Newtonian correction
      procedure :: gr_pos_kick => whm_gr_p4_tp           !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: accel       => whm_kick_getacch_tp    !! Compute heliocentric accelerations of test particles
      procedure :: kick        => whm_kick_vh_tp         !! Kick heliocentric velocities of test particles
      procedure :: step        => whm_step_tp            !! Steps the particle forward one stepsize
      final     ::                whm_final_tp      !! Finalizes the WHM test particle object - deallocates all allocatables 
   end type whm_tp

   !> An abstract class for the WHM integrator nbody system 
   type, extends(swiftest_nbody_system) :: whm_nbody_system
   contains
      !> Replace the abstract procedures with concrete ones
      procedure :: initialize   => whm_util_setup_initialize_system !! Performs WHM-specific initilization steps, like calculating the Jacobi masses
      procedure :: step         => whm_step_system             !! Advance the WHM nbody system forward in time by one step
      final     ::                 whm_final_system       !! Finalizes the WHM nbody_system object - deallocates all allocatables 
   end type whm_nbody_system

   interface
      module subroutine whm_coord_h2j_pl(self, cb)
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_h2j_pl

      module subroutine whm_coord_j2h_pl(self, cb)
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_j2h_pl

      module subroutine whm_coord_vh2vj_pl(self, cb)
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_vh2vj_pl

      module subroutine whm_drift_pl(self, nbody_system, param, dt)
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine whm_drift_pl

      !> Get heliocentric accelration of massive bodies
      module subroutine whm_kick_getacch_pl(self, nbody_system, param, t, lbeg)
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! WHM nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine whm_kick_getacch_pl

      !> Get heliocentric accelration of the test particle
      module subroutine whm_kick_getacch_tp(self, nbody_system, param, t, lbeg)
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! WHM nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine whm_kick_getacch_tp

      module subroutine whm_kick_vh_pl(self, nbody_system, param, t, dt, lbeg)
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine whm_kick_vh_pl

      module subroutine whm_kick_vh_tp(self, nbody_system, param, t, dt, lbeg)
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine whm_kick_vh_tp

      pure module subroutine whm_gr_kick_getacch_pl(self, param)
         implicit none
         class(whm_pl),              intent(inout) :: self  !! WHM massive body particle data structure
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine whm_gr_kick_getacch_pl

      pure module subroutine whm_gr_kick_getacch_tp(self, param)
         implicit none
         class(whm_tp),              intent(inout) :: self  !! WHM test particle data structure
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine whm_gr_kick_getacch_tp

      pure module subroutine whm_gr_p4_pl(self, nbody_system, param, dt)
         implicit none
         class(whm_pl),                intent(inout) :: self  !! WHM massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt    !! Step size
      end subroutine whm_gr_p4_pl

      pure module subroutine whm_gr_p4_tp(self, nbody_system, param, dt)
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Step size
      end subroutine whm_gr_p4_tp

      !> Reads WHM massive body object in from file
      module subroutine whm_util_setup_pl(self, n, param)
         implicit none
         class(whm_pl),             intent(inout) :: self  !! WHM massive body objectobject
         integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine whm_util_setup_pl

      module subroutine whm_util_setup_initialize_system(self, system_history, param)
         implicit none
         class(whm_nbody_system),                 intent(inout) :: self            !! WHM nbody system object
         class(swiftest_storage),    allocatable, intent(inout) :: system_history  !! Stores the system history between output dumps
         class(swiftest_parameters),              intent(inout) :: param           !! Current run configuration parameters 
      end subroutine whm_util_setup_initialize_system

      module subroutine whm_step_pl(self, nbody_system, param, t, dt)
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t       !! Simulation time
         real(DP),                     intent(in)    :: dt     !! Current stepsize
      end subroutine whm_step_pl

      module subroutine whm_step_system(self, param, t, dt)
         implicit none
         class(whm_nbody_system),    intent(inout) :: self    !! WHM nbody_system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine whm_step_system

      module subroutine whm_step_tp(self, nbody_system, param, t, dt)
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine whm_step_tp

      module subroutine whm_util_append_pl(self, source, lsource_mask)
         implicit none
         class(whm_pl),                   intent(inout) :: self         !! WHM massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine whm_util_append_pl

      module subroutine whm_util_dealloc_pl(self)
         implicit none
         class(whm_pl),  intent(inout) :: self !! WHM massive body object
      end subroutine whm_util_dealloc_pl

      module subroutine whm_util_spill_pl(self, discards, lspill_list, ldestructive)
         implicit none
         class(whm_pl),         intent(inout) :: self        !! WHM massive body object
         class(swiftest_body),      intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine whm_util_spill_pl

      module subroutine whm_util_fill_pl(self, inserts, lfill_list)
         implicit none
         class(whm_pl),         intent(inout) :: self       !! WHM massive body object
         class(swiftest_body),  intent(in)    :: inserts    !! inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine whm_util_fill_pl

      module subroutine whm_util_resize_pl(self, nnew)
         implicit none
         class(whm_pl), intent(inout) :: self  !! WHM massive body object
         integer(I4B),  intent(in)    :: nnew  !! New size neded
      end subroutine whm_util_resize_pl

      module subroutine whm_util_set_ir3j(self)
         implicit none
         class(whm_pl),                intent(inout) :: self    !! WHM massive body object
      end subroutine whm_util_set_ir3j

      module subroutine whm_util_set_mu_eta_pl(self, cb)
         implicit none
         class(whm_pl),                intent(inout) :: self    !! WHM massive body object
         class(swiftest_cb),           intent(inout) :: cb     !! Swiftest central body object
      end subroutine whm_util_set_mu_eta_pl

      module subroutine whm_util_sort_pl(self, sortby, ascending)
         implicit none
         class(whm_pl), intent(inout) :: self        !! WHM massive body object
         character(*),  intent(in)    :: sortby    !! Sorting attribute
         logical,       intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine whm_util_sort_pl

      module subroutine whm_util_sort_rearrange_pl(self, ind)
         implicit none
         class(whm_pl),               intent(inout) :: self !! WHM massive body object
         integer(I4B),  dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine whm_util_sort_rearrange_pl
   end interface

#ifdef COARRAY
   interface
      module subroutine whm_coarray_coclone_pl(self)
         implicit none
         class(whm_pl),intent(inout),codimension[*]  :: self  !! WHM pl object
      end subroutine whm_coarray_coclone_pl
   end interface
#endif

   contains

      subroutine whm_final_pl(self)
         !! author: David A. Minton
         !!
         !! Finalize the WHM massive body object - deallocates all allocatables
         implicit none
         ! Argument
         type(whm_pl),  intent(inout) :: self !! WHM massive body object

         call self%dealloc()

         return
      end subroutine whm_final_pl


      subroutine whm_final_system(self)
         !! author: David A. Minton
         !!
         !! Finalize the WHM nbody system object - deallocates all allocatables
         implicit none
         ! Arguments
         type(whm_nbody_system),  intent(inout) :: self !! WHM nbody system object

         call self%dealloc()

         return
      end subroutine whm_final_system


      subroutine whm_final_tp(self)
         !! author: David A. Minton
         !!
         !! Finalize the WHM test particle object - deallocates all allocatables
         implicit none
         ! Arguments
         type(whm_tp),  intent(inout) :: self !! WHM test particle object

         call self%dealloc()

         return
      end subroutine whm_final_tp

end module whm
