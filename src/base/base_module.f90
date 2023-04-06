!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

module base
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Base type definitions. This allows the collision and encounter modules to be defined before the swiftest module.
   !!
   use globals
#ifdef COARRAY
   use coarray
#endif
   implicit none
   public

   !> User defined parameters that are read in from the parameters input file. 
   !>    Each paramter is initialized to a default values. 
   type, abstract :: base_parameters
      character(STRMAX)                       :: integrator                             !! Name of the nbody integrator used
      character(STRMAX)                       :: param_file_name                        !! The name of the parameter file
      real(DP)                                :: t0                   =  0.0_DP         !! Integration reference time
      real(DP)                                :: tstart               = -1.0_DP         !! Integration start time
      real(DP)                                :: tstop                = -1.0_DP         !! Integration stop time
      real(DP)                                :: dt                   = -1.0_DP         !! Time step
      integer(I8B)                            :: iloop                = 0_I8B           !! Main loop counter
      integer(I8B)                            :: nloops               = 0_I8B           !! Total number of loops to execute
      character(STRMAX)                       :: incbfile             = CB_INFILE       !! Name of input file for the central body
      character(STRMAX)                       :: inplfile             = PL_INFILE       !! Name of input file for massive bodies
      character(STRMAX)                       :: intpfile             = TP_INFILE       !! Name of input file for test particles
      character(STRMAX)                       :: nc_in                = NC_INFILE       !! Name of system input file for NetCDF input
      character(STRMAX)                       :: in_type              = "NETCDF_DOUBLE" !! Data representation type of input data files
      character(STRMAX)                       :: in_form              = "XV"            !! Format of input data files ("EL" or ["XV"])
      integer(I4B)                            :: istep_out            = -1              !! Number of time steps between saved outputs
      integer(I4B)                            :: nstep_out            = -1              !! Total number of saved outputs
      real(DP)                                :: fstep_out            = 1.0_DP          !! The output step time stretching factor
      logical                                 :: ltstretch            = .false.         !! Whether to employ time stretching or not
      character(STRMAX)                       :: outfile              = BIN_OUTFILE     !! Name of output binary file
      character(STRMAX)                       :: out_type             = "NETCDF_DOUBLE" !! Binary format of output file
      character(STRMAX)                       :: out_form             = "XVEL"          !! Data to write to output file
      character(STRMAX)                       :: out_stat             = 'NEW'           !! Open status for output binary file
      integer(I4B)                            :: dump_cadence         =  10             !! Number of output steps between dumping simulation data to file
      real(DP)                                :: rmin                 = -1.0_DP         !! Minimum heliocentric radius for test particle
      real(DP)                                :: rmax                 = -1.0_DP         !! Maximum heliocentric radius for test particle
      real(DP)                                :: rmaxu                = -1.0_DP         !! Maximum unbound heliocentric radius for test particle
      real(DP)                                :: qmin                 = -1.0_DP         !! Minimum pericenter distance for test particle
      character(STRMAX)                       :: qmin_coord           = "HELIO"         !! Coordinate frame to use for qmin (["HELIO"] or "BARY")
      real(DP)                                :: qmin_alo             = -1.0_DP         !! Minimum semimajor axis for qmin
      real(DP)                                :: qmin_ahi             = -1.0_DP         !! Maximum semimajor axis for qmin
      real(QP)                                :: MU2KG                = -1.0_QP         !! Converts mass units to grams
      real(QP)                                :: TU2S                 = -1.0_QP         !! Converts time units to seconds
      real(QP)                                :: DU2M                 = -1.0_QP         !! Converts distance unit to centimeters
      real(DP)                                :: GU                   = -1.0_DP         !! Universal gravitational constant in the system units
      real(DP)                                :: inv_c2               = -1.0_DP         !! Inverse speed of light squared in the system units
      real(DP)                                :: GMTINY               = -1.0_DP         !! Smallest G*mass that is fully gravitating
      real(DP)                                :: min_GMfrag           = -1.0_DP         !! Smallest G*mass that can be produced in a fragmentation event
      real(DP)                                :: nfrag_reduction      =  30.0_DP        !! Reduction factor for limiting the number of fragments in a collision
      integer(I4B), dimension(:), allocatable :: seed                                   !! Random seeds for fragmentation modeling
      logical                                 :: lmtiny_pl            = .false.         !! Include semi-interacting massive bodies
      character(STRMAX)                       :: collision_model      = "MERGE"         !! The Coll
      character(STRMAX)                       :: encounter_save       = "NONE"          !! Indicate if and how encounter data should be saved
      logical                                 :: lenc_save_trajectory = .false.         !! Indicates that when encounters are saved, the full trajectory through recursion steps are saved
      logical                                 :: lenc_save_closest    = .false.         !! Indicates that when encounters are saved, the closest approach distance between pairs of bodies is saved
      character(NAMELEN)                      :: interaction_loops    = "ADAPTIVE"      !! Method used to compute interaction loops. Options are "TRIANGULAR", "FLAT", or "ADAPTIVE" 
      character(NAMELEN)                      :: encounter_check_plpl = "ADAPTIVE"      !! Method used to compute pl-pl encounter checks. Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE" 
      character(NAMELEN)                      :: encounter_check_pltp = "ADAPTIVE"      !! Method used to compute pl-tp encounter checks. Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE" 

      ! The following are used internally, and are not set by the user, but instead are determined by the input value of INTERACTION_LOOPS
      logical :: lflatten_interactions     = .false. !! Use the flattened upper triangular matrix for pl-pl interaction loops
      logical :: lencounter_sas_plpl       = .false. !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters
      logical :: lencounter_sas_pltp       = .false. !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters

      ! Logical flags to turn on or off various features of the code
      logical :: lrhill_present = .false. !! Hill radii are given as an input rather than calculated by the code (can be used to inflate close encounter regions manually)
      logical :: lextra_force   = .false. !! User defined force function turned on
      logical :: lbig_discard   = .false. !! Save big bodies on every discard
      logical :: lclose         = .false. !! Turn on close encounters
      logical :: lenergy        = .false. !! Track the total energy of the system
      logical :: loblatecb      = .false. !! Calculate acceleration from oblate central body (automatically turns true if nonzero J2 is input)
      logical :: lrotation      = .false. !! Include rotation states of big bodies
      logical :: ltides         = .false. !! Include tidal dissipation 

      ! Initial values to pass to the energy report subroutine (usually only used in the case of a restart, otherwise these will be updated with initial conditions values)
      real(DP)                  :: E_orbit_orig  = 0.0_DP  !! Initial orbital energy
      real(DP)                  :: GMtot_orig   = 0.0_DP  !! Initial system mass
      real(DP), dimension(NDIM) :: L_total_orig    = 0.0_DP  !! Initial total angular momentum vector
      real(DP), dimension(NDIM) :: L_orbit_orig  = 0.0_DP  !! Initial orbital angular momentum
      real(DP), dimension(NDIM) :: L_spin_orig   = 0.0_DP  !! Initial spin angular momentum vector
      real(DP), dimension(NDIM) :: L_escape      = 0.0_DP  !! Angular momentum of bodies that escaped the system (used for bookeeping)
      real(DP)                  :: GMescape     = 0.0_DP  !! Mass of bodies that escaped the system (used for bookeeping)
      real(DP)                  :: E_collisions  = 0.0_DP  !! Energy lost from system due to collisions
      real(DP)                  :: E_untracked   = 0.0_DP  !! Energy gained from system due to escaped bodies
      logical                   :: lfirstenergy = .true.  !! This is the first time computing energe
      logical                   :: lfirstkick   = .true.  !! Initiate the first kick in a symplectic step
      logical                   :: lrestart     = .false. !! Indicates whether or not this is a restarted run

      character(NAMELEN)        :: display_style         !! Style of the output display {"STANDARD", "COMPACT"}). Default is "STANDARD"
      integer(I4B)             :: display_unit          !! File unit number for display (either to stdout or to a log file)
      logical                  :: log_output  = .false. !! Logs the output to file instead of displaying it on the terminal

      ! Future features not implemented or in development
      logical :: lgr        = .false. !! Turn on GR
      logical :: lyarkovsky = .false. !! Turn on Yarkovsky effect
      logical :: lyorp      = .false. !! Turn on YORP effect
   contains
      procedure :: dealloc => base_util_dealloc_param
      procedure(abstract_io_dump_param),      deferred :: dump
      procedure(abstract_io_param_reader),    deferred :: reader
      procedure(abstract_io_param_writer),    deferred :: writer    
      procedure(abstract_io_read_in_param),   deferred :: read_in 
#ifdef COARRAY
      procedure :: coclone => base_coclone_param
#endif
   end type base_parameters

   abstract interface
      subroutine abstract_io_dump_param(self, param_file_name)
         import base_parameters
         implicit none
         class(base_parameters),intent(in)    :: self            !! Output collection of parameters
         character(len=*),      intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine abstract_io_dump_param

      subroutine abstract_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         import base_parameters, I4B
         implicit none
         class(base_parameters), intent(inout) :: self       !! Collection of parameters
         integer(I4B),           intent(in)    :: unit       !! File unit number
         character(len=*),       intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                             !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         character(len=*),       intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
         integer(I4B),           intent(out)   :: iostat     !! IO status code
         character(len=*),       intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      end subroutine abstract_io_param_reader

      subroutine abstract_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         import base_parameters, I4B
         implicit none
         class(base_parameters), intent(in)    :: self      !! Collection of parameters
         integer(I4B),           intent(in)    :: unit      !! File unit number
         character(len=*),       intent(in)    :: iotype    !! Dummy argument passed to the  input/output procedure contains the text from the char-literal-constant, prefixed with DT. 
                                                            !!    If you do not include a char-literal-constant, the iotype argument contains only DT.
         integer(I4B),           intent(in)    :: v_list(:) !! Not used in this procedure
         integer(I4B),           intent(out)   :: iostat    !! IO status code
         character(len=*),       intent(inout) :: iomsg     !! Message to pass if iostat /= 0
      end subroutine abstract_io_param_writer

      subroutine abstract_io_read_in_param(self, param_file_name) 
         import base_parameters
         implicit none
         class(base_parameters), intent(inout) :: self            !! Current run configuration parameters
         character(len=*),       intent(in)    :: param_file_name !! Parameter input file name (i.e. param.in)
      end subroutine abstract_io_read_in_param
   end interface


   type :: base_storage_frame
      class(*), allocatable :: item
   contains
      procedure :: store         => base_util_copy_store       !! Stores a snapshot of the nbody system so that later it can be retrieved for saving to file.
      generic   :: assignment(=) => store
      final     ::                  base_final_storage_frame
   end type


   type, abstract :: base_storage
      !! An class that establishes the pattern for various storage objects
      integer(I4B) :: nframes !! Total number of frames that can be stored

      !! An class that establishes the pattern for various storage objects
      type(base_storage_frame), dimension(:), allocatable :: frame          !! Array of stored frames
      integer(I4B)                                        :: iframe = 0     !! Index of the last frame stored in the system
      integer(I4B)                                        :: nid            !! Number of unique id values in all saved snapshots
      integer(I4B),             dimension(:), allocatable :: idvals         !! The set of unique id values contained in the snapshots
      integer(I4B),             dimension(:), allocatable :: idmap          !! The id value -> index map  
      integer(I4B)                                        :: nt             !! Number of unique time values in all saved snapshots
      real(DP),                 dimension(:), allocatable :: tvals          !! The set of unique time values contained in the snapshots
      integer(I4B),             dimension(:), allocatable :: tmap           !! The t value -> index map
   contains
      procedure :: dealloc => base_util_dealloc_storage !! Deallocates all allocatables
      procedure :: reset   => base_util_reset_storage   !! Resets the storage object back to its original state by removing all of the saved items from the storage frames
      procedure :: resize  => base_util_resize_storage  !! Resizes storage if it is too small 
      procedure :: setup   => base_util_setup_storage   !! Sets up a storage system with a set number of frames
      procedure :: save    => base_util_snapshot_save   !! Takes a snapshot of the current system
   end type base_storage


   !> Class definition for the particle origin information object. This object is used to track time, location, and collisional regime
   !> of fragments produced in collisional events.
   type, abstract :: base_particle_info
   end type base_particle_info


   !> An abstract class for a generic collection of Swiftest bodies
   type, abstract :: base_object
   contains
      procedure(abstract_util_dealloc_object), deferred :: dealloc
   end type base_object

   abstract interface
      subroutine abstract_util_dealloc_object(self)
         import base_object
         implicit none
         class(base_object), intent(inout) :: self !! Generic Swiftest object type
      end subroutine abstract_util_dealloc_object
   end interface


   !> Class definition for the kinship relationships used in bookkeeping multiple collisions bodies in a single time step.
   type, abstract :: base_kinship
   end type base_kinship
      

   !> An abstract class for a basic Swiftest nbody system 
   type, abstract :: base_nbody_system
   end type base_nbody_system

   contains

      subroutine base_util_copy_store(self, source)
         !! author: David A. Minton
         !!
         !! Stores a snapshot of the nbody system so that later it can be retrieved for saving to file.
         implicit none
         class(base_storage_frame),  intent(inout) :: self   !! Swiftest storage frame object
         class(*),                   intent(in)    :: source !! Swiftest n-body system object

         if (allocated(self%item)) deallocate(self%item)
         allocate(self%item, source=source)
         
         return
      end subroutine base_util_copy_store 


      subroutine base_util_dealloc_param(self)
         !! author: David A. Minton
         !!
         !! Deallocates all allocatables
         implicit none
         ! Arguments
         class(base_parameters),intent(inout)  :: self  !! Collection of parameters

         if (allocated(self%seed)) deallocate(self%seed)

         return
      end subroutine base_util_dealloc_param


      subroutine base_util_dealloc_storage(self)
         !! author: David A. Minton
         !!
         !! Resets a storage object by deallocating all items and resetting the frame counter to 0
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self !! Swiftest storage object

         call self%reset()
         if (allocated(self%frame)) deallocate(self%frame)
         self%nframes = 0
   
         return
      end subroutine base_util_dealloc_storage


      subroutine base_util_exit(code)
         !! author: David A. Minton
         !!
         !! Print termination message and exit program
         !!
         !! Adapted from David E. Kaufmann's Swifter routine: base_util_exit.f90
         !! Adapted from Hal Levison's Swift routine base_util_exit.f
         implicit none
         ! Arguments
         integer(I4B), intent(in) :: code
         ! Internals
         character(*), parameter :: BAR = '("------------------------------------------------")'
         character(*), parameter :: SUCCESS_MSG = '(/, "Normal termination of Swiftest (version ", f3.1, ")")'
         character(*), parameter :: FAIL_MSG = '(/, "Terminating Swiftest (version ", f3.1, ") due to error!!")'
         character(*), parameter :: USAGE_MSG = '("Usage: swiftest [bs|helio|ra15|rmvs|symba|tu4|whm] <paramfile> [standard|compact|progress|NONE]")'
         character(*), parameter :: HELP_MSG  = USAGE_MSG
   
         select case(code)
         case(SUCCESS)
            write(*, SUCCESS_MSG) VERSION_NUMBER
            write(*, BAR)
         case(USAGE) 
            write(*, USAGE_MSG)
         case(HELP)
            write(*, HELP_MSG)
         case default
            write(*, FAIL_MSG) VERSION_NUMBER
            write(*, BAR)
            error stop
         end select
   
         stop
   
      end subroutine base_util_exit


      subroutine base_util_reset_storage(self)
         !! author: David A. Minton
         !!
         !! Resets the storage object back to its original state by removing all of the saved items from the storage frames, but does not deallocate the frames
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self
         ! Internals
         integer(I4B) :: i
 
         if (allocated(self%frame)) then
            do i = 1, self%nframes
               if (allocated(self%frame(i)%item)) deallocate(self%frame(i)%item)
            end do
         end if

         if (allocated(self%idmap)) deallocate(self%idmap)
         if (allocated(self%idvals)) deallocate(self%idvals)
         if (allocated(self%tmap)) deallocate(self%tmap)
         if (allocated(self%tvals)) deallocate(self%tvals)
         self%nid = 0
         self%nt = 0
         self%iframe = 0

         return
      end subroutine base_util_reset_storage 

      
      subroutine base_util_resize_storage(self, nnew)
         !! author: David A. Minton
         !!
         !! Checks the current size of a Swiftest against the requested size and resizes it if it is too small.
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self !! Storage object
         integer(I4B),        intent(in)    :: nnew !! New size
         ! Internals
         class(base_storage_frame), dimension(:), allocatable :: tmp
         integer(I4B) :: i, nold, nbig
  
         nold = self%nframes
         if (nnew <= nold) return
   
         nbig = nold
         do while (nbig < nnew)
            nbig = nbig * 2
         end do
         call move_alloc(self%frame, tmp)
         allocate(self%frame(nbig))
         self%nframes = nbig
         do i = 1, nold
            if (allocated(tmp(i)%item)) call move_alloc(tmp(i)%item, self%frame(i)%item)
         end do
   
         return
      end subroutine base_util_resize_storage 


      subroutine base_util_setup_storage(self, n)
         !! author: David A. Minton
         !!
         !! Checks the current size of a Swiftest against the requested size and resizes it if it is too small.
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self !! Storage object
         integer(I4B),        intent(in)    :: n    !! New size

         if (allocated(self%frame)) deallocate(self%frame)
         allocate(self%frame(n))
         self%nframes = n
   
         return
      end subroutine base_util_setup_storage 
  
      
      subroutine base_util_snapshot_save(self, snapshot)
         !! author: David A. Minton
         !!
         !! Checks the current size of the storage object against the required size and extends it by a factor of 2 more than requested if it is too small.
         !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
         !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
         !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self     !! Storage ncounter storage object
         class(*),            intent(in)    :: snapshot !! Object to snapshot
         ! Internals
         integer(I4B) :: nnew, nold

         ! Advance the snapshot frame counter
         self%iframe = self%iframe + 1

         nnew = self%iframe
         nold = self%nframes

         ! Check to make sure the current storage  object is big enough. If not, grow it by a factor of 2
         if (nnew > nold) call self%resize(nnew)

         self%frame(nnew) = snapshot
      
         return
      end subroutine base_util_snapshot_save
   

      subroutine base_final_storage(self)
         !! author: David A. Minton
         !!
         !! Finalizer for the storage object
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self

         call self%dealloc()
         return
      end subroutine base_final_storage


      subroutine base_final_storage_frame(self)
         !! author: David A. Minton
         !!
         !! Finalizer for the storage frame data type
         implicit none
         type(base_storage_frame) :: self
   
         if (allocated(self%item)) deallocate(self%item)
   
         return
      end subroutine base_final_storage_frame

#ifdef COARRAY
      subroutine base_coclone_param(self)
         !! author: David A. Minton
         !!
         !! Broadcasts the image 1 parameter to all other images in a parameter coarray
         implicit none
         ! Arguments
         class(base_parameters),intent(inout),codimension[*]  :: self  !! Collection of parameters
         ! Internals
         integer(I4B) :: i

         call coclone(self%integrator)
         call coclone(self%param_file_name)
         call coclone(self%t0)
         call coclone(self%tstart)
         call coclone(self%tstop)
         call coclone(self%dt)
         call coclone(self%iloop)
         call coclone(self%nloops)
         call coclone(self%incbfile)
         call coclone(self%inplfile)
         call coclone(self%intpfile)
         call coclone(self%nc_in)
         call coclone(self%in_type)
         call coclone(self%in_form)
         call coclone(self%istep_out)
         call coclone(self%nstep_out)
         call coclone(self%fstep_out)
         call coclone(self%ltstretch)
         call coclone(self%outfile)
         call coclone(self%out_type)
         call coclone(self%out_form)
         call coclone(self%out_stat)
         call coclone(self%dump_cadence)
         call coclone(self%rmin)
         call coclone(self%rmax)
         call coclone(self%rmaxu)
         call coclone(self%qmin)
         call coclone(self%qmin_coord)
         call coclone(self%qmin_alo)
         call coclone(self%qmin_ahi)
         call coclone(self%MU2KG)
         call coclone(self%TU2S)
         call coclone(self%DU2M)
         call coclone(self%GU)
         call coclone(self%inv_c2)
         call coclone(self%GMTINY)
         call coclone(self%min_GMfrag)
         call coclone(self%nfrag_reduction)
         call coclone(self%lmtiny_pl)
         call coclone(self%collision_model)
         call coclone(self%encounter_save)
         call coclone(self%lenc_save_trajectory)
         call coclone(self%lenc_save_closest   )
         call coclone(self%interaction_loops   )
         call coclone(self%encounter_check_plpl)
         call coclone(self%encounter_check_pltp)
         call coclone(self%lflatten_interactions)
         call coclone(self%lencounter_sas_plpl)
         call coclone(self%lencounter_sas_pltp      )
         call coclone(self%lrhill_present)
         call coclone(self%lextra_force  )
         call coclone(self%lbig_discard  )
         call coclone(self%lclose        )
         call coclone(self%lenergy       )
         call coclone(self%loblatecb     )
         call coclone(self%lrotation     )
         call coclone(self%ltides        )
         call coclone(self%E_orbit_orig )
         call coclone(self%GMtot_orig  )
         do i = 1, NDIM
            call coclone(self%L_total_orig(i))
            call coclone(self%L_orbit_orig(i))
            call coclone(self%L_spin_orig(i))
            call coclone(self%L_escape(i))
         end do
         call coclone(self%GMescape    )
         call coclone(self%E_collisions )
         call coclone(self%E_untracked  )
         call coclone(self%lfirstenergy)
         call coclone(self%lfirstkick  )
         call coclone(self%lrestart    )
         call coclone(self%display_style)
         call coclone(self%display_unit )
         call coclone(self%log_output )
         call coclone(self%lgr       )
         call coclone(self%lyarkovsky)
         call coclone(self%lyorp     )
         call coclone(self%seed)

         return
      end subroutine base_coclone_param 

#endif
    

end module base
