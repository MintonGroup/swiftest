module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Regularized Mixed Variable Symplectic (RMVS) integrator
   !! Partially adapted from David E. Kaufmann's Swifter module: module_rmvs.f90
   use swiftest_globals
   use whm_classes, only : whm_cb, whm_pl, whm_tp, whm_nbody_system
   implicit none
   public

   integer(I4B), private, parameter :: NTENC = 10
   integer(I4B), private, parameter :: NTPHENC = 3
   integer(I4B), private, parameter :: NTPENC = NTENC * NTPHENC
   real(DP),     private, parameter :: RHSCALE = 3.5_DP
   real(DP),     private, parameter :: RHPSCALE = 1.0_DP
   real(DP),     private, parameter :: FACQDT = 2.0_DP

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, extends(whm_nbody_system) :: rmvs_nbody_system
      !> In the RMVS integrator, only test particles are discarded
      logical                               :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
      real(DP)                              :: rts                        !! fraction of Hill's sphere radius to use as radius of encounter region
      real(DP), dimension(:,:), allocatable :: vbeg                       !! Planet velocities at beginning ot step
   contains
      !> Replace the abstract procedures with concrete ones
      procedure :: initialize => rmvs_setup_initialize_system  !! Performs RMVS-specific initilization steps, including generating the close encounter planetocentric structures
      procedure :: step       => rmvs_step_system              !! Advance the RMVS nbody system forward in time by one step
   end type rmvs_nbody_system

   type, private :: rmvs_interp
      real(DP), dimension(:, :), allocatable :: x     !! interpolated heliocentric planet position for outer encounter
      real(DP), dimension(:, :), allocatable :: v     !! interpolated heliocentric planet velocity for outer encounter
      real(DP), dimension(:, :), allocatable :: aobl  !! Encountering planet's oblateness acceleration value
      real(DP), dimension(:, :), allocatable :: atide !! Encountering planet's tidal acceleration value
   end type rmvs_interp

   !********************************************************************************************************************************
   ! rmvs_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> RMVS central body particle class 
   type, extends(whm_cb) :: rmvs_cb
      type(rmvs_interp), dimension(:), allocatable :: outer !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp), dimension(:), allocatable :: inner !! interpolated heliocentric central body position for inner encounters
      logical                                      :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   end type rmvs_cb

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_util_spill_tp
      ! encounter steps)
      logical,      dimension(:),   allocatable :: lperi  !! planetocentric pericenter passage flag (persistent for a full rmvs time step) over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plperP !! index of planet associated with pericenter distance peri (persistent over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plencP !! index of planet that test particle is encountering (not persistent for a full RMVS time step)

      ! The following are used to correctly set the oblateness values of the acceleration during an inner encounter with a planet
      type(rmvs_cb)                             :: cb_heliocentric !! Copy of original central body object passed to close encounter (used for oblateness acceleration during planetocentric encoountters)
      real(DP),     dimension(:,:), allocatable :: xheliocentric   !! original heliocentric position (used for oblateness calculation during close encounters)
      integer(I4B)                              :: index           !!  inner substep number within current set
      integer(I4B)                              :: ipleP           !!  index value of encountering planet
      logical                                   :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   contains
      procedure :: discard         => rmvs_discard_tp             !! Check to see if test particles should be discarded based on pericenter passage distances with respect to planets encountered
      procedure :: encounter_check => rmvs_encounter_check_tp     !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure :: accel           => rmvs_kick_getacch_tp        !! Calculates either the standard or modified version of the acceleration depending if the
                                                                  !!    if the test particle is undergoing a close encounter or not
      procedure :: setup           => rmvs_setup_tp               !! Constructor method - Allocates space for number of particles
      procedure :: append          => rmvs_util_append_tp         !! Appends elements from one structure to another
      procedure :: fill            => rmvs_util_fill_tp           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize          => rmvs_util_resize_tp         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: sort            => rmvs_util_sort_tp           !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => rmvs_util_sort_rearrange_tp !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => rmvs_util_spill_tp          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_tp

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************
 
   !> RMVS massive body particle class
   type, extends(whm_pl) :: rmvs_pl
      integer(I4B),             dimension(:), allocatable :: nenc                      !! number of test particles encountering planet this full rmvs time step
      integer(I4B),             dimension(:), allocatable :: tpenc1P                   !! index of first test particle encountering planet
      integer(I4B),             dimension(:), allocatable :: plind                     !! Connects the planetocentric indices back to the heliocentric planet list
      type(rmvs_interp),        dimension(:), allocatable :: outer                     !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp),        dimension(:), allocatable :: inner                     !! interpolated heliocentric central body position for inner encounters
      class(rmvs_nbody_system), dimension(:), allocatable :: planetocentric            !! Planetocentric version of the massive body objects (one for each massive body)
      logical                                             :: lplanetocentric = .false. !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   contains
      procedure :: setup     => rmvs_setup_pl               !! Constructor method - Allocates space for number of particles
      procedure :: append    => rmvs_util_append_pl         !! Appends elements from one structure to another
      procedure :: fill      => rmvs_util_fill_pl           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize    => rmvs_util_resize_pl         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: sort      => rmvs_util_sort_pl           !! Sorts body arrays by a sortable componen
      procedure :: rearrange => rmvs_util_sort_rearrange_pl !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill     => rmvs_util_spill_pl          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_pl

   interface
      module elemental function rmvs_chk_ind(r2, v2, vdotr, dt, r2crit) result(lflag)
         implicit none
         real(DP), intent(in)       :: r2, v2, vdotr, dt, r2crit
         logical                    :: lflag
      end function rmvs_chk_ind

      module subroutine rmvs_discard_tp(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine rmvs_discard_tp

      module function rmvs_encounter_check_tp(self, system, dt) result(lencounter)
         implicit none
         class(rmvs_tp),           intent(inout) :: self       !! RMVS test particle object  
         class(rmvs_nbody_system), intent(inout) :: system     !! RMVS nbody system object
         real(DP),                 intent(in)    :: dt         !! step size
         logical                                 :: lencounter !! Returns true if there is at least one close encounter      
      end function rmvs_encounter_check_tp

      module subroutine rmvs_kick_getacch_tp(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structuree 
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine rmvs_kick_getacch_tp

      module subroutine rmvs_setup_pl(self, n, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_pl),             intent(inout) :: self  !! RMVS massive body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine rmvs_setup_pl

      module subroutine rmvs_setup_initialize_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),   intent(inout) :: self    !! RMVS system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      end subroutine rmvs_setup_initialize_system

      module subroutine rmvs_setup_tp(self, n, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_tp),            intent(inout) :: self !! RMVS test particle object
         integer(I4B),              intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parametere
      end subroutine rmvs_setup_tp

      module subroutine rmvs_util_append_pl(self, source, lsource_mask)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),                  intent(inout) :: self         !! RMVS massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine rmvs_util_append_pl

      module subroutine rmvs_util_append_tp(self, source, lsource_mask)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),                  intent(inout) :: self         !! RMVS test particle object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine rmvs_util_append_tp

      module subroutine rmvs_util_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object 
         class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_util_fill_pl

      module subroutine rmvs_util_fill_tp(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self        !! RMVS massive body object
         class(swiftest_body),  intent(in)    :: inserts     !!  Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_util_fill_tp

      module subroutine rmvs_util_resize_pl(self, nnew)
         implicit none
         class(rmvs_pl), intent(inout) :: self  !! RMVS massive body object
         integer(I4B),   intent(in)    :: nnew  !! New size neded
      end subroutine rmvs_util_resize_pl

      module subroutine rmvs_util_resize_tp(self, nnew)
         implicit none
         class(rmvs_tp), intent(inout) :: self  !! RMVS test particle object
         integer(I4B),   intent(in)    :: nnew  !! New size neded
      end subroutine rmvs_util_resize_tp

      module subroutine rmvs_util_sort_pl(self, sortby, ascending)
         implicit none
         class(rmvs_pl), intent(inout) :: self      !! RMVS massive body object
         character(*),   intent(in)    :: sortby    !! Sorting attribute
         logical,        intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine rmvs_util_sort_pl   

      module subroutine rmvs_util_sort_tp(self, sortby, ascending)
         implicit none
         class(rmvs_tp), intent(inout) :: self      !! RMVS test particle object
         character(*),   intent(in)    :: sortby    !! Sorting attribute
         logical,        intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      end subroutine rmvs_util_sort_tp

      module subroutine rmvs_util_sort_rearrange_pl(self, ind)
         implicit none
         class(rmvs_pl),               intent(inout) :: self !! RMVS massive body object
         integer(I4B),   dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine rmvs_util_sort_rearrange_pl

      module subroutine rmvs_util_sort_rearrange_tp(self, ind)
         implicit none
         class(rmvs_tp),                intent(inout) :: self !! RMVS test particle object
         integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      end subroutine rmvs_util_sort_rearrange_tp

      module subroutine rmvs_util_spill_pl(self, discards, lspill_list, ldestructive)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),        intent(inout) :: self         !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine rmvs_util_spill_pl

      module subroutine rmvs_util_spill_tp(self, discards, lspill_list, ldestructive)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self         !! RMVS test particle object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine rmvs_util_spill_tp

      module subroutine rmvs_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),   intent(inout) :: self    !! RMVS nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine rmvs_step_system

   end interface

end module rmvs_classes
