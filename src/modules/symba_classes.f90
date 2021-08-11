module symba_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic SyMBAcentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use swiftest_classes, only : swiftest_parameters, swiftest_base, swiftest_encounter
   use helio_classes,    only : helio_cb, helio_pl, helio_tp, helio_nbody_system
   use rmvs_classes,     only : rmvs_chk_ind
   implicit none
   public

   integer(I4B), private, parameter :: NENMAX           = 32767
   integer(I4B), private, parameter :: NTENC            = 3
   real(DP),     private, parameter :: RHSCALE          = 6.5_DP
   real(DP),     private, parameter :: RSHELL           = 0.48075_DP
   character(*),          parameter :: PARTICLE_OUTFILE = 'particle.dat'
   integer(I4B),          parameter :: PARTICLEUNIT     = 44 !! File unit number for the binary particle info output file

   type, extends(swiftest_parameters) :: symba_parameters
      character(STRMAX)                       :: particle_out  = PARTICLE_OUTFILE !! Name of output particle information file
      real(DP)                                :: GMTINY          = -1.0_DP          !! Smallest mass that is fully gravitating
      integer(I4B), dimension(:), allocatable :: seed                              !! Random seeds
      logical                                 :: lfragmentation = .false.          !! Do fragmentation modeling instead of simple merger.
   contains
      procedure :: reader => symba_io_param_reader
      procedure :: writer => symba_io_param_writer
   end type symba_parameters

   !********************************************************************************************************************************
   !                                    symba_particle_info class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the particle origin information object. This object is used to track time, location, and collisional regime
   !> of fragments produced in collisional events.
   type :: symba_particle_info
      sequence
      character(len=32)         :: origin_type !! String containing a description of the origin of the particle (e.g. Initial Conditions, Supercatastrophic, Disruption, etc.)
      real(DP)                  :: origin_time !! The time of the particle's formation
      real(DP), dimension(NDIM) :: origin_xh   !! The heliocentric distance vector at the time of the particle's formation
      real(DP), dimension(NDIM) :: origin_vh   !! The heliocentric velocity vector at the time of the particle's formation
   end type symba_particle_info

   !********************************************************************************************************************************
   !                                    symba_kinship class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the kinship relationships used in bookkeeping multiple collisions bodies in a single time step.
   type symba_kinship
      integer(I4B)                            :: parent !! Index of parent particle
      integer(I4B)                            :: nchild !! number of children in merger list
      integer(I4B), dimension(:), allocatable :: child  !! Index of children particles
   end type symba_kinship

   !********************************************************************************************************************************
   ! symba_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA central body particle class
   type, extends(helio_cb) :: symba_cb
      real(DP) :: M0  = 0.0_DP !! Initial mass of the central body
      real(DP) :: dM  = 0.0_DP !! Change in mass of the central body
      real(DP) :: R0  = 0.0_DP !! Initial radius of the central body
      real(DP) :: dR  = 0.0_DP !! Change in the radius of the central body
      type(symba_particle_info) :: info
   contains
   end type symba_cb

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA massive body class
   type, extends(helio_pl) :: symba_pl
      logical,                   dimension(:), allocatable :: lcollision !! flag indicating whether body has merged with another this time step
      logical,                   dimension(:), allocatable :: lencounter !! flag indicating whether body is part of an encounter this time step
      logical,                   dimension(:), allocatable :: lmtiny     !! flag indicating whether this body is below the GMTINY cutoff value
      integer(I4B)                                         :: nplm       !! number of bodies above the GMTINY limit
      integer(I8B)                                         :: nplplm     !! Number of body (all massive)-body (only those above GMTINY) comparisons in the flattened upper triangular matrix 
      integer(I4B),              dimension(:), allocatable :: nplenc     !! number of encounters with other planets this time step
      integer(I4B),              dimension(:), allocatable :: ntpenc     !! number of encounters with test particles this time step
      integer(I4B),              dimension(:), allocatable :: levelg     !! level at which this body should be moved
      integer(I4B),              dimension(:), allocatable :: levelm     !! deepest encounter level achieved this time step
      integer(I4B),              dimension(:), allocatable :: isperi     !! perihelion passage flag
      real(DP),                  dimension(:), allocatable :: peri       !! perihelion distance
      real(DP),                  dimension(:), allocatable :: atp        !! semimajor axis following perihelion passage
      type(symba_kinship),       dimension(:), allocatable :: kin        !! Array of merger relationship structures that can account for multiple pairwise mergers in a single step
      type(symba_particle_info), dimension(:), allocatable :: info
   contains
      procedure :: make_family     => symba_collision_make_family_pl !! When a single body is involved in more than one collision in a single step, it becomes part of a family
      procedure :: discard         => symba_discard_pl               !! Process massive body discards
      procedure :: drift           => symba_drift_pl                 !! Method for Danby drift in Democratic Heliocentric coordinates. Sets the mask to the current recursion level
      procedure :: encounter_check => symba_encounter_check_pl       !! Checks if massive bodies are going through close encounters with each other
      procedure :: accel           => symba_kick_getacch_pl          !! Compute heliocentric accelerations of massive bodies
      procedure :: setup           => symba_setup_pl                 !! Constructor method - Allocates space for the input number of bodies
      procedure :: append          => symba_util_append_pl           !! Appends elements from one structure to another
      procedure :: fill            => symba_util_fill_pl             !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: get_peri        => symba_util_peri_pl             !! Determine system pericenter passages for massive bodies
      procedure :: rearray         => symba_util_rearray_pl          !! Clean up the massive body structures to remove discarded bodies and add new bodies
      procedure :: resize          => symba_util_resize_pl           !! Checks the current size of a SyMBA massive body against the requested size and resizes it if it is too small.
      procedure :: sort            => symba_util_sort_pl             !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => symba_util_sort_rearrange_pl   !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => symba_util_spill_pl            !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type symba_pl

   type, extends(symba_pl) :: symba_merger
      integer(I4B), dimension(:), allocatable :: ncomp
   contains
      procedure :: append          => symba_util_append_merger       !! Appends elements from one structure to another
      procedure :: resize          => symba_util_resize_merger       !! Checks the current size of a SyMBA merger list against the requested size and resizes it if it is too small.
      procedure :: setup           => symba_setup_merger             !! Constructor method - Allocates space for the input number of bodies
   end type symba_merger

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA test particle class
   type, extends(helio_tp) :: symba_tp
      integer(I4B), dimension(:), allocatable :: nplenc  !! number of encounters with planets this time step
      integer(I4B), dimension(:), allocatable :: levelg  !! level at which this particle should be moved
      integer(I4B), dimension(:), allocatable :: levelm  !! deepest encounter level achieved this time step
      type(symba_particle_info), dimension(:), allocatable :: info
   contains
      procedure :: drift           => symba_drift_tp               !! Method for Danby drift in Democratic Heliocentric coordinates. Sets the mask to the current recursion level
      procedure :: encounter_check => symba_encounter_check_tp     !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure :: accel           => symba_kick_getacch_tp        !! Compute heliocentric accelerations of test particles
      procedure :: setup           => symba_setup_tp               !! Constructor method - Allocates space for the input number of bodies
      procedure :: append          => symba_util_append_tp         !! Appends elements from one structure to another
      procedure :: fill            => symba_util_fill_tp           !! "Fills" bodies from one object into another depending on the results of a mask (uses the UNPACK intrinsic)
      procedure :: resize          => symba_util_resize_tp         !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      procedure :: sort            => symba_util_sort_tp           !! Sorts body arrays by a sortable componen
      procedure :: rearrange       => symba_util_sort_rearrange_tp !! Rearranges the order of array elements of body based on an input index array. Used in sorting methods
      procedure :: spill           => symba_util_spill_tp          !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type symba_tp

   !********************************************************************************************************************************
   !                                    symba_pltpenc class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA class for tracking pl-tp close encounters in a step
   type, extends(swiftest_encounter) :: symba_pltpenc
      integer(I4B), dimension(:),   allocatable :: level  !! encounter recursion level
   contains
      procedure :: collision_check => symba_collision_check_pltpenc !! Checks if a test particle is going to collide with a massive body
      procedure :: encounter_check => symba_encounter_check_pltpenc !! Checks if massive bodies are going through close encounters with each other
      procedure :: kick            => symba_kick_pltpenc            !! Kick barycentric velocities of active test particles within SyMBA recursion
      procedure :: setup           => symba_setup_pltpenc           !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      procedure :: spill           => symba_util_spill_pltpenc      !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type symba_pltpenc

   !********************************************************************************************************************************
   !                                    symba_plplenc class definitions and method interfaces
   !*******************************************************************************************************************************
   !> SyMBA class for tracking pl-pl close encounters in a step
   type, extends(symba_pltpenc) :: symba_plplenc
   contains
      procedure :: extract_collisions     => symba_collision_encounter_extract_collisions        !! Processes the pl-pl encounter list remove only those encounters that led to a collision
      procedure :: resolve_fragmentations => symba_collision_resolve_fragmentations !! Process list of collisions, determine the collisional regime, and then create fragments
      procedure :: resolve_mergers        => symba_collision_resolve_mergers        !! Process list of collisions and merge colliding bodies together
   end type symba_plplenc

   !********************************************************************************************************************************
   !  symba_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, extends(helio_nbody_system) :: symba_nbody_system
      class(symba_merger),  allocatable :: pl_adds !! List of added bodies in mergers or collisions
      class(symba_pltpenc), allocatable :: pltpenc_list  !! List of massive body-test particle encounters in a single step 
      class(symba_plplenc), allocatable :: plplenc_list  !! List of massive body-massive body encounters in a single step
      class(symba_plplenc), allocatable :: plplcollision_list  !! List of massive body-massive body collisions in a single step
      integer(I4B)                      :: irec          !! System recursion level
   contains
      procedure :: write_discard    => symba_io_write_discard             !! Write out information about discarded and merged planets and test particles in SyMBA
      procedure :: initialize       => symba_setup_initialize_system      !! Performs SyMBA-specific initilization steps
      procedure :: step             => symba_step_system                  !! Advance the SyMBA nbody system forward in time by one step
      procedure :: interp           => symba_step_interp_system           !! Perform an interpolation step on the SymBA nbody system 
      procedure :: set_recur_levels => symba_step_set_recur_levels_system !! Sets recursion levels of bodies and encounter lists to the current system level
      procedure :: recursive_step   => symba_step_recur_system            !! Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current recursion level, if applicable, and descend to the next deeper level if necessary
      procedure :: reset            => symba_step_reset_system            !! Resets pl, tp,and encounter structures at the start of a new step 
   end type symba_nbody_system

   interface
      module subroutine symba_collision_check_pltpenc(self, system, param, t, dt, irec)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_pltpenc),       intent(inout) :: self   !! SyMBA pl-tp encounter list object
         class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
         class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! current time
         real(DP),                   intent(in)    :: dt     !! step size
         integer(I4B),               intent(in)    :: irec   !! Current recursion level
      end subroutine symba_collision_check_pltpenc

      module subroutine symba_collision_encounter_extract_collisions(self, system, param)
         implicit none
         class(symba_plplenc),       intent(inout) :: self   !! SyMBA pl-pl encounter list
         class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
         class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      end subroutine

      module subroutine symba_collision_make_family_pl(self,idx)
         implicit none
         class(symba_pl),            intent(inout) :: self !! SyMBA massive body object
         integer(I4B), dimension(2), intent(in)    :: idx !! Array holding the indices of the two bodies involved in the collision
      end subroutine symba_collision_make_family_pl

      module subroutine symba_collision_resolve_fragmentations(self, system, param)
         implicit none
         class(symba_plplenc),      intent(inout) :: self   !! SyMBA pl-pl encounter list
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
         class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters with SyMBA additions
      end subroutine symba_collision_resolve_fragmentations
   
      module subroutine symba_collision_resolve_mergers(self, system, param)
         implicit none
         class(symba_plplenc),      intent(inout) :: self   !! SyMBA pl-pl encounter list
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
         class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters with SyMBA additions
      end subroutine symba_collision_resolve_mergers

      module subroutine symba_discard_pl(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      end subroutine symba_discard_pl

      module subroutine symba_drift_pl(self, system, param, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine symba_drift_pl
   
      module subroutine symba_drift_tp(self, system, param, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_tp),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine symba_drift_tp

      module pure elemental subroutine symba_encounter_check_one(xr, yr, zr, vxr, vyr, vzr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
         implicit none
         real(DP),     intent(in)  :: xr, yr, zr, vxr, vyr, vzr
         real(DP),     intent(in)  :: rhill1, rhill2, dt
         integer(I4B), intent(in)  :: irec
         logical,      intent(out) :: lencounter, lvdotr
      end subroutine symba_encounter_check_one

      module function symba_encounter_check_pl(self, system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_pl),           intent(inout) :: self       !! SyMBA test particle object  
         class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt         !! step size
         integer(I4B),              intent(in)    :: irec       !! Current recursion level 
         logical                                  :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_pl

      module function symba_encounter_check_pltpenc(self, system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_pltpenc),      intent(inout) :: self           !! SyMBA pl-pl encounter list object
         class(symba_nbody_system), intent(inout) :: system         !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt             !! step size
         integer(I4B),              intent(in)    :: irec           !! Current recursion level 
         logical                                  :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_pltpenc

      module function symba_encounter_check_tp(self, system, dt, irec) result(lany_encounter)
         implicit none
         class(symba_tp),           intent(inout) :: self           !! SyMBA test particle object  
         class(symba_nbody_system), intent(inout) :: system         !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt             !! step size
         integer(I4B),              intent(in)    :: irec           !! Current recursion level 
         logical                                  :: lany_encounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_tp

      module function symba_fragmentation_casedisruption(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
         implicit none
         class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
         class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
         integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
         real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
         real(DP),                        intent(inout) :: Qloss            !! Energy lost during collisionn
         integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      end function symba_fragmentation_casedisruption
   
      module function symba_fragmentation_casehitandrun(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
         implicit none
         class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
         class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
         integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
         real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
         real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
         integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      end function symba_fragmentation_casehitandrun

      module function symba_fragmentation_casemerge(system, param, family, x, v, mass, radius, L_spin, Ip)  result(status)
         implicit none
         class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
         class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
         integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
         real(DP),        dimension(:,:), intent(in)    :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(in)    :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collisio
         integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      end function symba_fragmentation_casemerge

      module function symba_fragmentation_casesupercatastrophic(system, param, family, x, v, mass, radius, L_spin, Ip, mass_res, Qloss)  result(status)
         implicit none
         class(symba_nbody_system),       intent(inout) :: system           !! SyMBA nbody system object
         class(symba_parameters),         intent(in)    :: param            !! Current run configuration parameters with SyMBA additions
         integer(I4B),    dimension(:),   intent(in)    :: family           !! List of indices of all bodies inovlved in the collision
         real(DP),        dimension(:,:), intent(inout) :: x, v, L_spin, Ip !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass, radius     !! Input values that represent a 2-body equivalent of a possibly 2+ body collision
         real(DP),        dimension(:),   intent(inout) :: mass_res         !! The distribution of fragment mass obtained by the regime calculation 
         real(DP),                        intent(inout) :: Qloss            !! Energy lost during collision
         integer(I4B)                                   :: status           !! Status flag assigned to this outcome
      end function symba_fragmentation_casesupercatastrophic

      module subroutine symba_io_write_discard(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine symba_io_write_discard

      module subroutine symba_io_dump_particle_info(system, param, lincludecb, tpidx, plidx) 
         implicit none
         class(symba_nbody_system),             intent(inout) :: system !! SyMBA nbody system object
         class(symba_parameters),               intent(in)    :: param  !! Current run configuration parameters with SyMBA extensions
         logical,                     optional, intent(in)    :: lincludecb  !! Set to true to include the central body (default is false)
         integer(I4B), dimension(:),  optional, intent(in)    :: tpidx  !! Array of test particle indices to append to the particle file
         integer(I4B), dimension(:),  optional, intent(in)    :: plidx  !! Array of massive body indices to append to the particle file
      end subroutine symba_io_dump_particle_info
   
      module subroutine symba_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(symba_parameters), intent(inout) :: self       !! Current run configuration parameters with SyMBA additionss
         integer,                 intent(in)    :: unit       !! File unit number
         character(len=*),        intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                              !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer,                 intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
         integer,                 intent(out)   :: iostat     !! IO status code
         character(len=*),        intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      end subroutine symba_io_param_reader
   
      module subroutine symba_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         implicit none
         class(symba_parameters),intent(in)    :: self      !! Current run configuration parameters with SyMBA additions
         integer,                intent(in)    :: unit      !! File unit number
         character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                            !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer,                intent(in)    :: v_list(:) !! Not used in this procedure
         integer,                intent(out)   :: iostat    !! IO status code
         character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine symba_io_param_writer

      module subroutine symba_io_read_particle(system, param)
         implicit none
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system file
         class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters with SyMBA extensions
      end subroutine symba_io_read_particle

      module subroutine symba_kick_getacch_pl(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine symba_kick_getacch_pl

      module subroutine symba_kick_getacch_tp(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_tp),              intent(inout) :: self   !! SyMBA test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine symba_kick_getacch_tp

      module subroutine symba_kick_pltpenc(self, system, dt, irec, sgn)
         implicit none
         class(symba_pltpenc),      intent(in)    :: self   !! SyMBA pl-tp encounter list object
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
         real(DP),                  intent(in)    :: dt     !! step size
         integer(I4B),              intent(in)    :: irec   !! Current recursion level
         integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
      end subroutine symba_kick_pltpenc
   
      module subroutine symba_setup_initialize_particle_info(system, param) 
         implicit none
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
         class(symba_parameters),   intent(inout) :: param  !! Current run configuration parameters with SyMBA extensions
      end subroutine symba_setup_initialize_particle_info

      module subroutine symba_setup_initialize_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      end subroutine symba_setup_initialize_system

      module subroutine symba_setup_merger(self, n, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_merger),        intent(inout) :: self  !! SyMBA merger list object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine symba_setup_merger

      module subroutine symba_setup_pl(self, n, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine symba_setup_pl

      module subroutine symba_setup_pltpenc(self,n)
         implicit none
         class(symba_pltpenc), intent(inout) :: self !! SyMBA pl-tp encounter structure
         integer(I4B),         intent(in)    :: n    !! Number of encounters to allocate space for
      end subroutine symba_setup_pltpenc

      module subroutine symba_setup_tp(self, n, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_tp),            intent(inout) :: self  !! SyMBA test particle object
         integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      end subroutine symba_setup_tp

      module subroutine symba_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t     !! Simulation time
         real(DP),                   intent(in)    :: dt    !! Current stepsize
      end subroutine symba_step_system

      module subroutine symba_step_interp_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
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

      module recursive subroutine symba_step_recur_system(self, param, t, ireci)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
         class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
         real(DP),                   value         :: t
         integer(I4B),               value         :: ireci !! input recursion level
      end subroutine symba_step_recur_system

      module subroutine symba_step_reset_system(self)
         implicit none
         class(symba_nbody_system),  intent(inout) :: self !! SyMBA nbody system object
      end subroutine symba_step_reset_system
   end interface

   interface util_append
      module subroutine symba_util_append_arr_info(arr, source, nold, nsrc, lsource_mask)
         implicit none
         type(symba_particle_info), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         type(symba_particle_info), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                                         intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,                   dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_arr_info
     
      module subroutine symba_util_append_arr_kin(arr, source, nold, nsrc, lsource_mask)
         implicit none
         type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
         type(symba_kinship), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
         integer(I4B),                                   intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
         logical,             dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_arr_kin
   end interface

   interface
      module subroutine symba_util_append_merger(self, source, lsource_mask)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_merger),             intent(inout) :: self         !! SyMBA massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_merger

      module subroutine symba_util_append_pl(self, source, lsource_mask)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_pl),                 intent(inout) :: self         !! SyMBA massive body object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_pl

      module subroutine symba_util_append_tp(self, source, lsource_mask)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_tp),                 intent(inout) :: self        !! SyMBA test particle object
         class(swiftest_body),            intent(in)    :: source       !! Source object to append
         logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine symba_util_append_tp
   end interface 

   interface util_fill
      module subroutine symba_util_fill_arr_info(keeps, inserts, lfill_list)
         implicit none
         type(symba_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         type(symba_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,                   dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_arr_info

      module subroutine symba_util_fill_arr_kin(keeps, inserts, lfill_list)
         implicit none
         type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
         type(symba_kinship), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
         logical,             dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_arr_kin
   end interface

   interface
      module subroutine symba_util_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_pl),       intent(inout) :: self       !! SyMBA massive body object
         class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_pl

      module subroutine symba_util_fill_tp(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_tp),       intent(inout) :: self       !! SyMBA test particle object
         class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine symba_util_fill_tp

      module subroutine symba_util_peri_pl(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      end subroutine symba_util_peri_pl

      module subroutine symba_util_rearray_pl(self, system, param)
         implicit none
         class(symba_pl),           intent(inout) :: self   !! SyMBA massive body object
         class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
         class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters with SyMBA additions
      end subroutine symba_util_rearray_pl
   end interface

   interface util_resize
      module subroutine symba_util_resize_arr_info(arr, nnew)
         implicit none
         type(symba_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                         intent(in)    :: nnew !! New size
      end subroutine symba_util_resize_arr_info

      module subroutine symba_util_resize_arr_kin(arr, nnew)
         implicit none
         type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
         integer(I4B),                                   intent(in)    :: nnew !! New size
      end subroutine symba_util_resize_arr_kin
   end interface
   
   interface
      module subroutine symba_util_resize_merger(self, nnew)
         implicit none
         class(symba_merger), intent(inout) :: self  !! SyMBA merger list object
         integer(I4B),        intent(in)    :: nnew  !! New size neded
      end subroutine symba_util_resize_merger

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

   interface util_spill
      module subroutine symba_util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
         implicit none
         type(symba_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         type(symba_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,                    dimension(:),             intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                              intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_arr_info

      module subroutine symba_util_spill_arr_kin(keeps, discards, lspill_list, ldestructive)
         implicit none
         type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
         type(symba_kinship), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
         logical,               dimension(:),            intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
         logical,                                        intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_arr_kin
   end interface

   interface
      module subroutine symba_util_spill_pl(self, discards, lspill_list, ldestructive)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_pl),       intent(inout) :: self         !! SyMBA massive body object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_pl

      module subroutine symba_util_spill_pltpenc(self, discards, lspill_list, ldestructive)
         use swiftest_classes, only : swiftest_encounter
         implicit none
         class(symba_pltpenc),      intent(inout) :: self         !! SyMBA pl-tp encounter list
         class(swiftest_encounter), intent(inout) :: discards     !! Discarded object 
         logical, dimension(:),     intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,                   intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      end subroutine symba_util_spill_pltpenc

      module subroutine symba_util_spill_tp(self, discards, lspill_list, ldestructive)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(symba_tp),       intent(inout) :: self         !! SyMBA test particle object
         class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not
      end subroutine symba_util_spill_tp
   end interface

end module symba_classes