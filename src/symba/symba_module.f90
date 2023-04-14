!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module symba
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the SyMBA integrator
   !! Adapted from David E. Kaufmann's Swifter routine: module_symba.f90
   use swiftest
   use helio
   implicit none
   public

   integer(I4B), private, parameter :: NENMAX  = 32767
   integer(I4B), private, parameter :: NTENC   = 3
   real(DP),     private, parameter :: RHSCALE = 6.5_DP
   real(DP),     private, parameter :: RSHELL  = 0.48075_DP

   !> SyMBA central body particle class
   type, extends(helio_cb) :: symba_cb
   end type symba_cb


   !> SyMBA massive body class
   type, extends(helio_pl) :: symba_pl
      integer(I4B),              dimension(:), allocatable :: levelg     !! level at which this body should be moved
      integer(I4B),              dimension(:), allocatable :: levelm     !! deepest encounter level achieved this time step
   contains
      procedure :: flatten         => symba_util_flatten_eucl_plpl      !! Sets up the (i, j) -> k indexing used for the single-loop blocking Euclidean distance matrix
      procedure :: discard         => symba_discard_pl                  !! Process massive body discards
      procedure :: drift           => symba_drift_pl                    !! Method for Danby drift in Democratic Heliocentric coordinates. Sets the mask to the current recursion level
      procedure :: encounter_check => symba_encounter_check_pl          !! Checks if massive bodies are going through close encounters with each other
      procedure :: gr_pos_kick     => symba_gr_p4_pl                    !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: accel_int       => symba_kick_getacch_int_pl         !! Compute direct cross (third) term heliocentric accelerations of massive bodiess, with no mutual interactions between bodies below GMTINY
      procedure :: accel           => symba_kick_getacch_pl             !! Compute heliocentric accelerations of massive bodies
      procedure :: setup           => symba_util_setup_pl                    !! Constructor method - Allocates space for the input number of bodies
      procedure :: append          => symba_util_append_pl              !! Appends elements from one structure to another
      procedure :: dealloc         => symba_util_dealloc_pl             !! Deallocates all allocatable arrays
      procedure :: fill            => symba_util_fill_pl                !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize          => symba_util_resize_pl              !! Checks the current size of a SyMBA massive body against the requested size and resizes it if it is too small.
      procedure :: set_renc_I4B    => symba_util_set_renc               !! Sets the critical radius for encounter given an input recursion depth
      procedure :: sort            => symba_util_sort_pl                !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => symba_util_sort_rearrange_pl      !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => symba_util_spill_pl               !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type symba_pl


   !> SyMBA test particle class
   type, extends(helio_tp) :: symba_tp
      integer(I4B), dimension(:), allocatable :: levelg  !! level at which this particle should be moved
      integer(I4B), dimension(:), allocatable :: levelm  !! deepest encounter level achieved this time step
   contains
      procedure :: drift           => symba_drift_tp               !! Method for Danby drift in Democratic Heliocentric coordinates. Sets the mask to the current recursion level
      procedure :: encounter_check => symba_encounter_check_tp     !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure :: gr_pos_kick     => symba_gr_p4_tp               !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: accel           => symba_kick_getacch_tp        !! Compute heliocentric accelerations of test particles
      procedure :: setup           => symba_util_setup_tp               !! Constructor method - Allocates space for the input number of bodies
      procedure :: append          => symba_util_append_tp         !! Appends elements from one structure to another
      procedure :: dealloc         => symba_util_dealloc_tp        !! Deallocates all allocatable arrays
      procedure :: fill            => symba_util_fill_tp           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize          => symba_util_resize_tp         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: sort            => symba_util_sort_tp           !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => symba_util_sort_rearrange_tp !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => symba_util_spill_tp          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type symba_tp


   !> SyMBA class for tracking close encounters in a step
   type, extends(collision_list_plpl) :: symba_list_plpl
   contains
      procedure :: encounter_check => symba_encounter_check_list_plpl !! Checks if massive bodies are going through close encounters with each other
      procedure :: kick            => symba_kick_list_plpl            !! Kick barycentric velocities of active massive bodies within SyMBA recursion
   end type symba_list_plpl


   !> SyMBA class for tracking close encounters in a step
   type, extends(collision_list_pltp) :: symba_list_pltp
   contains
      procedure :: encounter_check => symba_encounter_check_list_pltp !! Checks if massive bodies are going through close encounters with test particles
      procedure :: kick            => symba_kick_list_pltp            !! Kick barycentric velocities of active test particles within SyMBA recursion
   end type symba_list_pltp


   type, extends(helio_nbody_system) :: symba_nbody_system
      integer(I4B)  :: irec = -1 !! nbody_system recursion level
   contains
      procedure :: dealloic         => symba_util_dealloc_system          !! Deallocates all allocatables
      procedure :: initialize       => symba_util_setup_initialize_system !! Performs SyMBA-specific initilization steps
      procedure :: step             => symba_step_system                  !! Advance the SyMBA nbody system forward in time by one step
      procedure :: interp           => symba_step_interp_system           !! Perform an interpolation step on the SymBA nbody system 
      procedure :: set_recur_levels => symba_step_set_recur_levels_system !! Sets recursion levels of bodies and encounter lists to the current nbody_system level
      procedure :: recursive_step   => symba_step_recur_system            !! Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current recursion level, if applicable, and descend to the next deeper level if necessary
      procedure :: reset            => symba_step_reset_system            !! Resets pl, tp,and encounter structures at the start of a new step 
   end type symba_nbody_system

   interface
      module subroutine symba_discard_pl(self, nbody_system, param)
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      end subroutine symba_discard_pl

      module subroutine symba_drift_pl(self, nbody_system, param, dt)
         implicit none
         class(symba_pl),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine symba_drift_pl
   
      module subroutine symba_drift_tp(self, nbody_system, param, dt)
         implicit none
         class(symba_tp),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine symba_drift_tp

      module function symba_encounter_check_pl(self, param, nbody_system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_pl),            intent(inout) :: self           !! SyMBA test particle object  
         class(swiftest_parameters), intent(inout) :: param          !! Current Swiftest run configuration parameters
         class(symba_nbody_system),  intent(inout) :: nbody_system         !! SyMBA nbody system object
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level
         logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_pl

      module function symba_encounter_check_list_plpl(self, param, nbody_system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_list_plpl),     intent(inout) :: self           !! SyMBA pl-pl encounter list object
         class(swiftest_parameters), intent(inout) :: param          !! Current Swiftest run configuration parameters
         class(symba_nbody_system),  intent(inout) :: nbody_system         !! SyMBA nbody system object
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level 
         logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_list_plpl

      module function symba_encounter_check_list_pltp(self, param, nbody_system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_list_pltp),     intent(inout) :: self           !! SyMBA pl-tp encounter list object
         class(swiftest_parameters), intent(inout) :: param          !! Current Swiftest run configuration parameters
         class(symba_nbody_system),  intent(inout) :: nbody_system         !! SyMBA nbody system object
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level 
         logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_list_pltp

      module function symba_encounter_check_tp(self, param, nbody_system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_tp),            intent(inout) :: self           !! SyMBA test particle object  
         class(swiftest_parameters), intent(inout) :: param          !! Current Swiftest run configuration parameters
         class(symba_nbody_system),  intent(inout) :: nbody_system         !! SyMBA nbody system object
         real(DP),                   intent(in)    :: dt             !! step size
         integer(I4B),               intent(in)    :: irec           !! Current recursion level 
         logical                                   :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_tp

      pure module subroutine symba_gr_p4_pl(self, nbody_system, param, dt)
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Step size
      end subroutine symba_gr_p4_pl
   
      pure module subroutine symba_gr_p4_tp(self, nbody_system, param, dt)
         implicit none
         class(symba_tp),              intent(inout) :: self   !! SyMBA test particle object
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Step size
      end subroutine symba_gr_p4_tp

      module subroutine symba_util_set_renc(self, scale)
         implicit none
         class(symba_pl), intent(inout) :: self !! SyMBA massive body object
         integer(I4B),    intent(in)    :: scale !! Current recursion depth
      end subroutine symba_util_set_renc
   
      module subroutine symba_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(swiftest_parameters),intent(in)    :: self      !! Current run configuration parameters with SyMBA additions
         integer,                intent(in)    :: unit      !! File unit number
         character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                            !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer,                intent(in)    :: v_list(:) !! Not used in this procedure
         integer,                intent(out)   :: iostat    !! IO status code
         character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine symba_io_param_writer

      module subroutine symba_kick_getacch_int_pl(self, param)
         implicit none
         class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current Swiftest run configuration parameters
      end subroutine symba_kick_getacch_int_pl

      module subroutine symba_kick_getacch_pl(self, nbody_system, param, t, lbeg)
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine symba_kick_getacch_pl

      module subroutine symba_kick_getacch_tp(self, nbody_system, param, t, lbeg)
         implicit none
         class(symba_tp),              intent(inout) :: self   !! SyMBA test particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine symba_kick_getacch_tp

      module subroutine symba_kick_list_plpl(self, nbody_system, dt, irec, sgn)
         implicit none
         class(symba_list_plpl),    intent(in)    :: self   !! SyMBA pl-tp encounter list object
         class(symba_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt     !! step size
         integer(I4B),              intent(in)    :: irec   !! Current recursion level
         integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      end subroutine symba_kick_list_plpl

      module subroutine symba_kick_list_pltp(self, nbody_system, dt, irec, sgn)
         implicit none
         class(symba_list_pltp),    intent(in)    :: self   !! SyMBA pl-tp encounter list object
         class(symba_nbody_system), intent(inout) :: nbody_system !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt     !! step size
         integer(I4B),              intent(in)    :: irec   !! Current recursion level
         integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      end subroutine symba_kick_list_pltp

      module subroutine symba_util_dealloc_system(self)
         implicit none
         class(symba_nbody_system), intent(inout) :: self
      end subroutine symba_util_dealloc_system


      module subroutine symba_util_setup_initialize_system(self, system_history, param)
         implicit none
         class(symba_nbody_system),               intent(inout) :: self           !! SyMBA nbody_system object
         class(swiftest_storage),    allocatable, intent(inout) :: system_history !! Stores the system history between output dumps
         class(swiftest_parameters),              intent(inout) :: param          !! Current run configuration parameters 
      end subroutine symba_util_setup_initialize_system

      module subroutine symba_util_setup_pl(self, n, param)
         implicit none
         class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine symba_util_setup_pl

      module subroutine symba_util_setup_tp(self, n, param)
         implicit none
         class(symba_tp),            intent(inout) :: self  !! SyMBA test particle object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      end subroutine symba_util_setup_tp

      module subroutine symba_step_system(self, param, t, dt)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
      end subroutine symba_step_system

      module subroutine symba_step_interp_system(self, param, t, dt)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
      end subroutine symba_step_interp_system

      module subroutine symba_step_set_recur_levels_system(self, ireci)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system objec
         integer(I4B),               intent(in)    :: ireci !! Input recursion level
      end subroutine symba_step_set_recur_levels_system

      recursive module subroutine symba_step_recur_system(self, param, t, ireci)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Current simulation time
         integer(I4B),               intent(in)    :: ireci !! input recursion level
      end subroutine symba_step_recur_system

      module subroutine symba_step_reset_system(self, param)
         implicit none
         class(symba_nbody_system), intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters with SyMBA additions
      end subroutine symba_step_reset_system

      module subroutine symba_util_append_pl(self, source, lsource_mask)
         implicit none
         class(symba_pl),                 intent(inout) :: self         !! SyMBA massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_pl

      module subroutine symba_util_append_tp(self, source, lsource_mask)
         implicit none
         class(symba_tp),                 intent(inout) :: self        !! SyMBA test particle object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_tp

      module subroutine symba_util_dealloc_pl(self)
         implicit none
         class(symba_pl),  intent(inout) :: self !! SyMBA massive body object
      end subroutine symba_util_dealloc_pl

      module subroutine symba_util_dealloc_tp(self)
         implicit none
         class(symba_tp),  intent(inout) :: self !! SyMBA test particle object
      end subroutine symba_util_dealloc_tp
   end interface 


   interface
      module subroutine symba_util_fill_pl(self, inserts, lfill_list)
         implicit none
         class(symba_pl),       intent(inout) :: self       !! SyMBA massive body object
         class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_pl

      module subroutine symba_util_fill_tp(self, inserts, lfill_list)
         implicit none
         class(symba_tp),       intent(inout) :: self       !! SyMBA test particle object
         class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_tp

      module subroutine symba_util_flatten_eucl_plpl(self, param)
         implicit none
         class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      end subroutine symba_util_flatten_eucl_plpl

   end interface

   interface
      module subroutine symba_util_resize_pl(self, nnew)
         implicit none
         class(symba_pl), intent(inout) :: self  !! SyMBA massive body object
         integer(I4B),    intent(in)    :: nnew  !! New size neded
      end subroutine symba_util_resize_pl

      module subroutine symba_util_resize_tp(self, nnew)
         implicit none
         class(symba_tp), intent(inout) :: self  !! SyMBA massive body object
         integer(I4B),    intent(in)    :: nnew  !! New size neded
      end subroutine symba_util_resize_tp

      module subroutine symba_util_sort_pl(self, sortby, ascending)
         implicit none
         class(symba_pl), intent(inout) :: self      !! SyMBA massive body object
         character(*),    intent(in)    :: sortby    !! Sorting attribute
         logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine symba_util_sort_pl 

      module subroutine symba_util_sort_tp(self, sortby, ascending)
         implicit none
         class(symba_tp), intent(inout) :: self      !! SyMBA test particle object
         character(*),    intent(in)    :: sortby    !! Sorting attribute
         logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine symba_util_sort_tp
   end interface

   interface
      module subroutine symba_util_sort_rearrange_pl(self, ind)
         implicit none
         class(symba_pl),               intent(inout) :: self !! SyMBA massive body object
         integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine symba_util_sort_rearrange_pl

      module subroutine symba_util_sort_rearrange_tp(self, ind)
         implicit none
         class(symba_tp),               intent(inout) :: self !! SyMBA massive body object
         integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine symba_util_sort_rearrange_tp
   end interface

   interface
      module subroutine symba_util_spill_pl(self, discards, lspill_list, ldestructive)
         implicit none
         class(symba_pl),       intent(inout) :: self         !! SyMBA massive body object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_pl

      module subroutine symba_util_spill_tp(self, discards, lspill_list, ldestructive)
         implicit none
         class(symba_tp),       intent(inout) :: self         !! SyMBA test particle object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_tp
   end interface

end module symba