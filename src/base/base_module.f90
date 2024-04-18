! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

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
      character(STRMAX) :: integrator 
         !! Name of the nbody integrator used  
      character(STRMAX) :: param_file_name                        
         !! The name of the parameter file  
      real(DP)  :: t0 =  0.0_DP 
         !! Integration reference time  
      real(DP) :: tstart = -1.0_DP
         !! Integration start time  
      real(DP) :: tstop = -1.0_DP
         !! Integration stop time  
      real(DP) :: dt = -1.0_DP
         !! Time step 
      integer(I8B) :: iloop = 0_I8B
         !! Main loop counter  
      integer(I8B) :: nloops = 0_I8B
         !! Total number of loops to execute  
      integer(I8B) :: istart = 0_I8B
         !! Starting index for loop counter  
      integer(I4B) :: iout = 0
         !! Output cadence counter  
      integer(I4B) :: idump = 0
         !! Dump cadence counter  
      integer(I4B)      :: nout                 = 0               
         !! Current output step  
      integer(I4B)      :: istep                = 0               
         !! Current value of istep (used for time stretching)  
      character(STRMAX) :: incbfile             = CB_INFILE       
         !! Name of input file for the central body  
      character(STRMAX) :: inplfile             = PL_INFILE       
         !! Name of input file for massive bodies  
      character(STRMAX) :: intpfile             = TP_INFILE       
         !! Name of input file for test particles  
      character(STRMAX) :: nc_in                = NC_INFILE       
         !! Name of system input file for NetCDF input  
      character(STRMAX) :: in_type              = "NETCDF_DOUBLE" 
         !! Data representation type of input data files  
      character(STRMAX) :: in_form              = "XV"            
         !! Format of input data files ("EL" or ["XV"])  
      integer(I4B)      :: istep_out            = -1              
         !! Number of time steps between saved outputs  
      integer(I4B)      :: nstep_out            = -1              
         !! Total number of saved outputs  
      real(DP)          :: fstep_out            = 1.0_DP          
         !! The output step time stretching factor  
      logical           :: ltstretch            = .false.         
         !! Whether to employ time stretching or not  
      character(STRMAX) :: outfile              = BIN_OUTFILE     
         !! Name of output binary file  
      character(STRMAX) :: out_type             = "NETCDF_DOUBLE" 
         !! Binary format of output file  
      character(STRMAX) :: out_form             = "XVEL"          
         !! Data to write to output file  
      character(STRMAX) :: out_stat             = 'NEW'           
         !! Open status for output binary file  
      integer(I4B)      :: dump_cadence         =  10             
         !! Number of output steps between dumping simulation data to file  
      real(DP)          :: rmin                 = -1.0_DP         
         !! Minimum heliocentric radius for test particle  
      real(DP)          :: rmax                 = -1.0_DP         
         !! Maximum heliocentric radius for test particle  
      real(DP)          :: rmaxu                = -1.0_DP         
         !! Maximum unbound heliocentric radius for test particle  
      real(DP)          :: qmin                 = -1.0_DP         
         !! Minimum pericenter distance for test particle 
      character(STRMAX) :: qmin_coord           = "HELIO"         
         !! Coordinate frame to use for qmin (["HELIO"] or "BARY") 
      real(DP)          :: qmin_alo             = -1.0_DP         
         !! Minimum semimajor axis for qmin 
      real(DP)          :: qmin_ahi             = -1.0_DP         
         !! Maximum semimajor axis for qmin 
      real(QP)          :: MU2KG                = -1.0_QP         
         !! Converts mass units to grams 
      real(QP)          :: TU2S                 = -1.0_QP         
         !! Converts time units to seconds 
      real(QP)          :: DU2M                 = -1.0_QP         
         !! Converts distance unit to centimeters 
      real(DP)          :: GU                   = -1.0_DP         
         !! Universal gravitational constant in the system units 
      real(DP)          :: inv_c2               = -1.0_DP         
         !! Inverse speed of light squared in the system units 
      real(DP)          :: GMTINY               = -1.0_DP         
         !! Smallest G*mass that is fully gravitating 
      real(DP)          :: min_GMfrag           = -1.0_DP         
         !! Smallest G*mass that can be produced in a fragmentation event 
      real(DP)          :: nfrag_reduction      =  30.0_DP        
         !! Reduction factor for limiting the number of collision fragments  
      integer(I4B), dimension(:), allocatable :: seed             
         !! Random seeds for fragmentation modeling 
      logical           :: lmtiny_pl            = .false.         
         !! Include semi-interacting massive bodies 
      character(STRMAX) :: collision_model      = "MERGE"         
         !! The Coll 
      character(STRMAX) :: encounter_save       = "NONE"          
         !! Indicate if and how encounter data should be saved 
      logical           :: lenc_save_trajectory = .false.         
         !! Indicates that when encounters are saved, the full trajectory through recursion steps are saved 
      logical           :: lenc_save_closest    = .false.         
         !! Indicates that when encounters are saved, the closest approach distance between pairs of bodies is saved 
      character(NAMELEN):: interaction_loops    = "ADAPTIVE"      
         !! Method used to compute interaction loops.  
         !!    Options are "TRIANGULAR", "FLAT", or "ADAPTIVE"  
      character(NAMELEN):: encounter_check_plpl = "ADAPTIVE"      
         !! Method used to compute pl-pl encounter checks.  
         !!    Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE"  
      character(NAMELEN):: encounter_check_pltp = "ADAPTIVE"      
         !! Method used to compute pl-tp encounter checks.  
         !!    Options are "TRIANGULAR", "SORTSWEEP", or "ADAPTIVE"  
      logical           :: lcoarray             = .false.         
         !! Use Coarrays for test particle parallelization. 

      ! The following are not set by the user, but instead are determined by the input value of INTERACTION_LOOPS
      logical :: lflatten_interactions     = .false. 
         !! Use the flattened upper triangular matrix for pl-pl interaction loops 
      logical :: lencounter_sas_plpl       = .false. 
         !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters 
      logical :: lencounter_sas_pltp       = .false. 
         !! Use the Sort and Sweep algorithm to prune the encounter list before checking for close encounters 

      ! Logical flags to turn on or off various features of the code
      logical :: lrhill_present = .false. 
         !! Hill radii are given as an input rather than calculated by the code (can be used to  inflate close encounter regions 
         !! manually) 
      logical :: lextra_force   = .false. 
         !! User defined force function turned on 
      logical :: lbig_discard   = .false. 
         !! Save big bodies on every discard 
      logical :: lclose         = .false. 
         !! Turn on close encounters 
      logical :: lenergy        = .false. 
         !! Track the total energy of the system 
      logical :: lnon_spherical_cb = .false. 
         !! Calculate acceleration from oblate central body (automatically turns true if nonzero J2, J4, or c_lm is input) 
      logical :: lrotation      = .false. 
         !! Include rotation states of big bodies 
      logical :: ltides         = .false. 
         !! Include tidal dissipation  

      ! Initial values to pass to the energy report subroutine (usually only used in the case of a restart, otherwise these will be 
      ! updated with initial conditions values)
      real(DP)                  :: E_orbit_orig = 0.0_DP  
         !! Initial orbital energy 
      real(DP)                  :: GMtot_orig   = 0.0_DP  
         !! Initial system mass 
      real(DP), dimension(NDIM) :: L_total_orig = 0.0_DP  
         !! Initial total angular momentum vector 
      real(DP), dimension(NDIM) :: L_orbit_orig = 0.0_DP  
         !! Initial orbital angular momentum 
      real(DP), dimension(NDIM) :: L_spin_orig  = 0.0_DP  
         !! Initial spin angular momentum vector 
      real(DP), dimension(NDIM) :: L_escape     = 0.0_DP  
         !! Angular momentum of escaped bodies (used for bookeeping) 
      real(DP)                  :: GMescape     = 0.0_DP  
         !! Mass of bodies that escaped the system (used for bookeeping) 
      real(DP)                  :: E_collisions = 0.0_DP  
         !! Energy lost from system due to collisions 
      real(DP)                  :: E_untracked  = 0.0_DP  
         !! Energy gained from system due to escaped bodies 
      logical                   :: lfirstenergy = .true.  
         !! This is the first time computing energe 
      logical                   :: lfirstkick   = .true.  
         !! Initiate the first kick in a symplectic step 
      logical                   :: lrestart     = .false. 
         !! Indicates whether or not this is a restarted run 

      character(NAMELEN)       :: display_style        
         !! Style of the output display {["STANDARD"], "COMPACT"}).  
      integer(I4B)             :: display_unit = OUTPUT_UNIT  
         !! File unit number for display (either to stdout or to a log file) 
      logical                  :: log_output  = .false. 
         !! Logs the output to file instead of displaying it on the terminal 

      ! Future features not implemented or in development
      logical :: lgr        = .false. 
         !! Turn on GR 
      logical :: lyarkovsky = .false. 
         !! Turn on Yarkovsky effect 
      logical :: lyorp      = .false. 
         !! Turn on YORP effect 
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
         class(base_parameters),intent(in)    :: self            
            !! Output collection of parameters 
         character(len=*),      intent(in)    :: param_file_name 
            !! Parameter input file name (i.e. param.in) 
      end subroutine abstract_io_dump_param

      subroutine abstract_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
         import base_parameters, I4B
         implicit none
         class(base_parameters), intent(inout) :: self       
            !! Collection of parameters 
         integer(I4B),           intent(in)    :: unit       
            !! File unit number 
         character(len=*),       intent(in)    :: iotype     
            !! Dummy argument passed to the  input/output procedure contains the  text from the char-literal-constant, prefixed with
            !! DT. If you do not include a char-literal-constant, the iotype argument contains  only DT. 
         character(len=*),       intent(in)    :: v_list(:)  
            !! The first element passes the integrator code to the reader 
         integer(I4B),           intent(out)   :: iostat     
               !! IO status code 
         character(len=*),       intent(inout) :: iomsg      
            !! Message to pass if iostat /= 0 
      end subroutine abstract_io_param_reader

      subroutine abstract_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
         import base_parameters, I4B
         implicit none
         class(base_parameters), intent(in)    :: self      
            !! Collection of parameters 
         integer(I4B),           intent(in)    :: unit      
            !! File unit number 
         character(len=*),       intent(in)    :: iotype    
            !! Dummy argument passed to the  input/output procedure contains the  
                                                            
            !!    text from the char-literal-constant, prefixed with DT. If you do  
                                                            
            !!    not include a char-literal-constant, the iotype argument contains  
                                                            
            !!    only DT. 
         integer(I4B),           intent(in)    :: v_list(:) 
            !! Not used in this procedure 
         integer(I4B),           intent(out)   :: iostat    
            !! IO status code 
         character(len=*),       intent(inout) :: iomsg     
            !! Message to pass if iostat /= 0 
      end subroutine abstract_io_param_writer

      subroutine abstract_io_read_in_param(self, param_file_name) 
         import base_parameters
         implicit none
         class(base_parameters), intent(inout) :: self            
            !! Current run configuration parameters 
         character(len=*),       intent(in)    :: param_file_name 
            !! Parameter input file name (i.e. param.in) 
      end subroutine abstract_io_read_in_param
   end interface


   type :: base_storage_frame
      class(*), allocatable :: item
   contains
      procedure :: store         => base_util_copy_store       
            !! Stores a snapshot of the nbody system so that later it can be  
                                                               
            !!    retrieved for saving to file. 
      generic   :: assignment(=) => store
      final     ::                  base_final_storage_frame
   end type


   type, abstract :: base_storage
      
         !! An class that establishes the pattern for various storage objects 
      integer(I4B) :: nframes 
         !! Total number of frames that can be stored An class that establishes the pattern for various storage objects 
      type(base_storage_frame), dimension(:), allocatable :: frame      
         !! Array of stored frames 
      integer(I4B)                                        :: iframe = 0 
         !! Index of the last frame stored in the system 
      integer(I4B)                                        :: nid        
         !! Number of unique id values in all saved snapshots 
      integer(I4B),             dimension(:), allocatable :: idvals     
         !! The set of unique id values contained in the snapshots 
      integer(I4B),             dimension(:), allocatable :: idmap      
         !! The id value -> index map   
      integer(I4B)                                        :: nt         
         !! Number of unique time values in all saved snapshots 
      real(DP),                 dimension(:), allocatable :: tvals      
         !! The set of unique time values contained in the snapshots 
      integer(I4B),             dimension(:), allocatable :: tmap       
         !! The t value -> index map 
   contains
      procedure :: dealloc => base_util_dealloc_storage 
         !! Deallocates all allocatables 
      procedure :: reset   => base_util_reset_storage   
         !! Resets the storage object back to its original state by removing all of the saved items from the storage frames 
      procedure :: resize  => base_util_resize_storage  
         !! Resizes storage if it is too small  
      procedure :: setup   => base_util_setup_storage   
         !! Sets up a storage system with a set number of frames 
      procedure :: save    => base_util_snapshot_save   
         !! Takes a snapshot of the current system 
   end type base_storage


   !> Class definition for the particle origin information object. This object is used to track time, location, and collisional 
   !> regime of fragments produced in collisional events.
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
         class(base_object), intent(inout) :: self 
            !! Generic Swiftest object type 
      end subroutine abstract_util_dealloc_object
   end interface


   !> Class definition for the kinship relationships used in bookkeeping multiple collisions bodies in a single time step.
   type, abstract :: base_kinship
   end type base_kinship
      

   !> An abstract class for a basic Swiftest nbody system 
   type, abstract :: base_nbody_system
   end type base_nbody_system


   interface util_append
      module procedure base_util_append_arr_char_string
      module procedure base_util_append_arr_DP
      module procedure base_util_append_arr_DPvec
      module procedure base_util_append_arr_I4B
      module procedure base_util_append_arr_logical
   end interface

   interface util_fill
      module procedure base_util_fill_arr_char_string
      module procedure base_util_fill_arr_DP
      module procedure base_util_fill_arr_DPvec
      module procedure base_util_fill_arr_I4B
      module procedure base_util_fill_arr_logical
   end interface

   interface util_resize
      module procedure base_util_resize_arr_char_string
      module procedure base_util_resize_arr_DP
      module procedure base_util_resize_arr_DPvec
      module procedure base_util_resize_arr_I4B
      module procedure base_util_resize_arr_logical
   end interface

   interface util_sort      
      module procedure base_util_sort_i4b
      module procedure base_util_sort_index_i4b
      module procedure base_util_sort_index_I4B_I8Bind
      module procedure base_util_sort_index_I8B_I8Bind
      module procedure base_util_sort_sp
      module procedure base_util_sort_index_sp
      module procedure base_util_sort_dp
      module procedure base_util_sort_index_dp
   end interface 

   interface util_sort_rearrange
      module procedure base_util_sort_rearrange_arr_char_string
      module procedure base_util_sort_rearrange_arr_DP
      module procedure base_util_sort_rearrange_arr_DPvec
      module procedure base_util_sort_rearrange_arr_I4B
      module procedure base_util_sort_rearrange_arr_I4B_I8Bind
      module procedure base_util_sort_rearrange_arr_logical
      module procedure base_util_sort_rearrange_arr_logical_I8Bind
   end interface

   interface util_spill
      module procedure base_util_spill_arr_char_string
      module procedure base_util_spill_arr_DP
      module procedure base_util_spill_arr_DPvec
      module procedure base_util_spill_arr_I4B
      module procedure base_util_spill_arr_I8B
      module procedure base_util_spill_arr_logical
   end interface

   interface util_unique
      module procedure base_util_unique_DP
      module procedure base_util_unique_I4B
   end interface

   contains

      subroutine base_util_append_arr_char_string(arr, source, nold, lsource_mask)
         
         !! author: David A. Minton 
         !! 
         !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big  
         !! enough, this will allocate space for it. 
         implicit none
         ! Arguments
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr    
            !! Destination array  
         character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source 
            !! Array to append  
         integer(I4B), intent(in), optional :: nold 
            !! Extent of original array. If passed, the source array will begin at  arr(nold+1). 
            !! Otherwise, the size of arr will be used. 
         logical, dimension(:), intent(in), optional :: lsource_mask 
            !! Logical mask indicating which elements to append to 
         ! Internals
         integer(I4B) :: nnew, nsrc, nend_orig

         if (.not.allocated(source)) return

         if (present(lsource_mask)) then
            nsrc = count(lsource_mask(:))
         else
            nsrc = size(source)
         end if
         if (nsrc == 0) return

         if (.not.allocated(arr)) then
            nend_orig = 0
            allocate(arr(nsrc))
         else
            if (present(nold)) then
               nend_orig = nold
            else
               nend_orig = size(arr)
            end if
            call util_resize(arr, nend_orig + nsrc)
         end if
         nnew = nend_orig + nsrc

         if (present(lsource_mask)) then
            arr(nend_orig + 1:nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))
         else
            arr(nend_orig + 1:nnew) = source(1:nsrc)
         end if

         return
      end subroutine base_util_append_arr_char_string


      subroutine base_util_append_arr_DP(arr, source, nold, lsource_mask)
         !! author: David A. Minton 
         !! 
         !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big  
         !! enough, this will allocate space for it. 
         implicit none
         ! Arguments
         real(DP), dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         real(DP), dimension(:), allocatable, intent(in) :: source 
            !! Array to append  
         integer(I4B), intent(in), optional :: nold 
            !! Extent of original array. If passed, the source array will begin at  
                                                    
            !!   arr(nold+1). Otherwise, the size of arr will be used. 
         logical,  dimension(:), intent(in), optional :: lsource_mask 
            !! Logical mask indicating which elements to append to 
         ! Internals
         integer(I4B) :: nnew, nsrc, nend_orig

         if (.not.allocated(source)) return

         if (present(lsource_mask)) then
            nsrc = count(lsource_mask(:))
         else
            nsrc = size(source)
         end if
         if (nsrc == 0) return

         if (.not.allocated(arr)) then
            nend_orig = 0
            allocate(arr(nsrc))
         else
            if (present(nold)) then
               nend_orig = nold
            else
               nend_orig = size(arr)
            end if
            call util_resize(arr, nend_orig + nsrc)
         end if
         nnew = nend_orig + nsrc

         if (present(lsource_mask)) then
            arr(nend_orig + 1:nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))
         else
            arr(nend_orig + 1:nnew) = source(1:nsrc)
         end if

         return
      end subroutine base_util_append_arr_DP


      subroutine base_util_append_arr_DPvec(arr, source, nold, lsource_mask)
         !! author: David A. Minton 
         !! 
         !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not  
         !! allocated, or is not big enough, this will allocate space for it. 
         implicit none
         ! Arguments
         real(DP), dimension(:,:), allocatable, intent(inout) :: arr       
            !! Destination array  
         real(DP), dimension(:,:), allocatable, intent(in) :: source       
            !! Array to append  
         integer(I4B), intent(in), optional :: nold 
            !! Extent of original array. If passed, the source array will begin at arr(nold+1). 
            !! Otherwise, the size of arr will be used. 
         logical,  dimension(:), intent(in), optional :: lsource_mask 
            !! Logical mask indicating which elements to append to 
         ! Internals
         integer(I4B) :: nnew, nsrc, nend_orig

         if (.not.allocated(source)) return

         if (present(lsource_mask)) then
            nsrc = count(lsource_mask(:))
         else
            nsrc = size(source,dim=2)
         end if
         if (nsrc == 0) return

         if (.not.allocated(arr)) then
            nend_orig = 0
            allocate(arr(NDIM,nsrc))
         else
            if (present(nold)) then
               nend_orig = nold
            else
               nend_orig = size(arr,dim=2)
            end if
            call util_resize(arr, nend_orig + nsrc)
         end if
         nnew = nend_orig + nsrc

         if (present(lsource_mask)) then
            arr(1, nend_orig + 1:nnew) = pack(source(1,1:nsrc), lsource_mask(1:nsrc))
            arr(2, nend_orig + 1:nnew) = pack(source(2,1:nsrc), lsource_mask(1:nsrc))
            arr(3, nend_orig + 1:nnew) = pack(source(3,1:nsrc), lsource_mask(1:nsrc))
         else
            arr(:,nend_orig + 1:nnew) = source(:,1:nsrc)
         end if

         return
      end subroutine base_util_append_arr_DPvec


      subroutine base_util_append_arr_I4B(arr, source, nold, lsource_mask)
         !! author: David A. Minton 
         !! 
         !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough,
         !! this will allocate space for it. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr          
            !! Destination array  
         integer(I4B), dimension(:), allocatable, intent(in)    :: source       
            !! Array to append  
         integer(I4B), intent(in), optional :: nold 
            !! Extent of original array. If passed, the source array will begin at arr(nold+1).
            !! Otherwise, the size of arr will be used. 
         logical, dimension(:), intent(in), optional :: lsource_mask 
            !! Logical mask indicating which elements to append to 
         ! Internals
         integer(I4B) :: nnew, nsrc, nend_orig

         if (.not.allocated(source)) return

         if (present(lsource_mask)) then
            nsrc = count(lsource_mask(:))
         else
            nsrc = size(source)
         end if
         if (nsrc == 0) return

         if (.not.allocated(arr)) then
            nend_orig = 0
            allocate(arr(nsrc))
         else
            if (present(nold)) then
               nend_orig = nold
            else
               nend_orig = size(arr)
            end if
            call util_resize(arr, nend_orig + nsrc)
         end if
         nnew = nend_orig + nsrc

         if (present(lsource_mask)) then
            arr(nend_orig + 1:nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))
         else
            arr(nend_orig + 1:nnew) = source(1:nsrc)
         end if

         return
      end subroutine base_util_append_arr_I4B


      subroutine base_util_append_arr_logical(arr, source, nold, lsource_mask)
         !! author: David A. Minton 
         !! 
         !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough,
         !! this will allocate space for it. 
         implicit none
         ! Arguments
         logical, dimension(:), allocatable, intent(inout) :: arr          
            !! Destination array  
         logical, dimension(:), allocatable, intent(in)    :: source       
            !! Array to append  
         integer(I4B), intent(in), optional :: nold   
            !! Extent of original array. If passed, the source array will begin at arr(nold+1). 
            !! Otherwise, the size of arr will be used. 
         logical, dimension(:), intent(in), optional :: lsource_mask 
            !! Logical mask indicating which elements to append to 
         ! Internals
         integer(I4B) :: nnew, nsrc, nend_orig

         if (.not.allocated(source)) return

         if (present(lsource_mask)) then
            nsrc = count(lsource_mask(:))
         else
            nsrc = size(source)
         end if
         if (nsrc == 0) return

         if (.not.allocated(arr)) then
            nend_orig = 0
            allocate(arr(nsrc))
         else
            if (present(nold)) then
               nend_orig = nold
            else
               nend_orig = size(arr)
            end if
            call util_resize(arr, nend_orig + nsrc)
         end if
         nnew = nend_orig + nsrc
         
         if (present(lsource_mask)) then
            arr(nend_orig + 1:nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))
         else
            arr(nend_orig + 1:nnew) = source(:)
         end if

         return
      end subroutine base_util_append_arr_logical


      subroutine base_util_copy_store(self, source)
         !! author: David A. Minton 
         !! 
         !! Stores a snapshot of the nbody system so that later it can be retrieved for saving to file. 
         implicit none
         class(base_storage_frame),  intent(inout) :: self   
            !! Swiftest storage frame object 
         class(*),                   intent(in)    :: source 
            !! Swiftest n-body system object 

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
         class(base_parameters),intent(inout)  :: self  
            !! Collection of parameters 

         if (allocated(self%seed)) deallocate(self%seed)

         return
      end subroutine base_util_dealloc_param


      subroutine base_util_dealloc_storage(self)
         !! author: David A. Minton 
         !! 
         !! Resets a storage object by deallocating all items and resetting the frame counter to 0 
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self 
            !! Swiftest storage object 

         call self%reset()
         if (allocated(self%frame)) deallocate(self%frame)
         self%nframes = 0
   
         return
      end subroutine base_util_dealloc_storage


      subroutine base_util_exit(code,unit)
         !! author: David A. Minton 
         !! 
         !! Print termination message and exit program 
         !! 
         !! Adapted from David E. Kaufmann's Swifter routine: base_util_exit.f90 
         !! Adapted from Hal Levison's Swift routine base_util_exit.f 
         implicit none
         ! Arguments
         integer(I4B), intent(in) :: code
         integer(I4B), intent(in), optional :: unit
         ! Internals
         character(*), parameter :: BAR = '("---------------------------------------------------")'
         character(*), parameter :: SUCCESS_MSG = '(/, "Normal termination of Swiftest (version ", A, ")")'
         character(*), parameter :: FAIL_MSG = '(/, "Terminating Swiftest (version ", A, ") due to error!")' 
         character(*), parameter :: USAGE_MSG = '("Usage: swiftest <whm|helio|rmvs|symba> <paramfile> ' // &
                                                '[{standard}|compact|progress]")'
         character(*), parameter :: HELP_MSG  = USAGE_MSG
         integer(I4B) :: iu

         if (present(unit)) then
            iu = unit
         else
            iu = OUTPUT_UNIT
         end if
   
         select case(code)
         case(SUCCESS)
            write(iu, SUCCESS_MSG) VERSION
            write(iu, BAR)
         case(USAGE) 
            write(iu, USAGE_MSG)
         case(HELP)
            write(iu, HELP_MSG)
         case default
            write(iu, FAIL_MSG) VERSION
            write(iu, BAR)
            stop 
         end select
   
         stop
   
      end subroutine base_util_exit


      subroutine base_util_fill_arr_char_string(keeps, inserts, lfill_list)
         !! author: David A. Minton 
         !! 
         !! Performs a fill operation on a single array of type character strings.
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps      
         !! Array of values to keep  
         character(len=STRMAX), dimension(:), allocatable, intent(in)    :: inserts    
         !! Array of values to insert into keep 
         logical,               dimension(:),              intent(in)    :: lfill_list 
         !! Logical array of bodies to merge into the  
                                                                                       
         !!    keeps 
   
         if (.not.allocated(keeps) .or. .not.allocated(inserts)) return
   
         keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
         keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
         return
      end subroutine base_util_fill_arr_char_string
   
   
      subroutine base_util_fill_arr_DP(keeps, inserts, lfill_list)
         !! author: David A. Minton 
         !! 
         !! Performs a fill operation on a single array of type DP.
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         real(DP), dimension(:), allocatable, intent(inout) :: keeps      
            !! Array of values to keep  
         real(DP), dimension(:), allocatable, intent(in)    :: inserts    
            !! Array of values to insert into keep 
         logical,  dimension(:),              intent(in)    :: lfill_list 
            !! Logical array of bodies to merge into the keeps 
   
         if (.not.allocated(keeps) .or. .not.allocated(inserts)) return
   
         keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
         keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
         return
      end subroutine base_util_fill_arr_DP
   
   
      subroutine base_util_fill_arr_DPvec(keeps, inserts, lfill_list)
         !! author: David A. Minton 
         !! 
         !! Performs a fill operation on a single array of DP vectors with shape (NDIM, n) 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         real(DP), dimension(:,:), allocatable, intent(inout) :: keeps      
            !! Array of values to keep  
         real(DP), dimension(:,:), allocatable, intent(in)    :: inserts    
            !! Array of values to insert into keep 
         logical,  dimension(:),                intent(in)    :: lfill_list 
            !! Logical array of bodies to merge into the keeps 
         ! Internals
         integer(I4B) :: i
   
         if (.not.allocated(keeps) .or. .not.allocated(inserts)) return
   
         do i = 1, NDIM
            keeps(i,:) = unpack(keeps(i,:),   .not.lfill_list(:), keeps(i,:))
            keeps(i,:) = unpack(inserts(i,:),      lfill_list(:), keeps(i,:))
         end do
   
         return
      end subroutine base_util_fill_arr_DPvec
   
   
      subroutine base_util_fill_arr_I4B(keeps, inserts, lfill_list)
         !! author: David A. Minton 
         !! 
         !! Performs a fill operation on a single array of type I4B 
                  !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: keeps      
            !! Array of values to keep  
         integer(I4B), dimension(:), allocatable, intent(in)    :: inserts    
            !! Array of values to insert into keep 
         logical,      dimension(:),              intent(in)    :: lfill_list 
            !! Logical array of bodies to merge into the keeps 
   
         if (.not.allocated(keeps) .or. .not.allocated(inserts)) return
   
         keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
         keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
         return
      end subroutine base_util_fill_arr_I4B
   
   
      subroutine base_util_fill_arr_logical(keeps, inserts, lfill_list)
         !! author: David A. Minton 
         !! 
         !! Performs a fill operation on a single array of logicals 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         logical, dimension(:), allocatable, intent(inout) :: keeps      
            !! Array of values to keep  
         logical, dimension(:), allocatable, intent(in)    :: inserts    
            !! Array of values to insert into keep 
         logical, dimension(:),              intent(in)    :: lfill_list 
            !! Logical array of bodies to merge into the keeps 
   
         if (.not.allocated(keeps) .or. .not.allocated(inserts)) return
   
         keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
         keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
         return
      end subroutine base_util_fill_arr_logical      

      subroutine base_util_reset_storage(self)
         !! author: David A. Minton 
         !! 
         !! Resets the storage object back to its original state by removing all of the saved items from the storage frames, but  
         !! does not deallocate the frames 
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


      subroutine base_util_resize_arr_char_string(arr, nnew)
         !! author: David A. Minton 
         !! 
         !! Resizes an array component of type character string. nnew = 0 will deallocate. 
         implicit none
         ! Arguments
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr  
            !! Array to resize 
         integer(I4B),                                     intent(in)    :: nnew 
            !! New size 
         ! Internals
         character(len=STRMAX), dimension(:), allocatable :: tmp 
            !! Temp. storage array in case the input array is already allocated 
         integer(I4B) :: nold 
            !! Old size 
   
         if (nnew < 0) return
   
         if (nnew == 0) then
            if (allocated(arr)) deallocate(arr)
            return
         end if
         
         if (allocated(arr)) then
            nold = size(arr)
         else
            nold = 0
         end if
   
         if (nnew == nold) return
         
         allocate(tmp(nnew))
         if (nold > 0) then
            if (nnew > nold) then
               tmp(1:nold) = arr(1:nold)
               tmp(nold+1:nnew) = ""
            else
               tmp(1:nnew) = arr(1:nnew)
            end if
         else
            tmp(1:nnew) = ""
         end if
         call move_alloc(tmp, arr)
   
         return
      end subroutine base_util_resize_arr_char_string
   
   
      subroutine base_util_resize_arr_DP(arr, nnew)
         !! author: David A. Minton 
         !! 
         !! Resizes an array component of double precision type. Passing nnew = 0 will deallocate. 
         implicit none
         ! Arguments
         real(DP), dimension(:), allocatable, intent(inout) :: arr  
            !! Array to resize 
         integer(I4B),                        intent(in)    :: nnew 
            !! New size 
         ! Internals
         real(DP), dimension(:), allocatable :: tmp 
            !! Temporary storage array in case the input array is already allocated 
         integer(I4B) :: nold 
            !! Old size 
         real(DP), parameter :: init_val = 0.0_DP
   
         if (nnew < 0) return
   
         if (nnew == 0) then
            if (allocated(arr)) deallocate(arr)
            return
         end if
         
         if (allocated(arr)) then
            nold = size(arr)
         else
            nold = 0
         end if
   
         if (nnew == nold) return
         
         allocate(tmp(nnew))
         if (nold > 0) then
            if (nnew > nold) then
               tmp(1:nold) = arr(1:nold)
               tmp(nold+1:nnew) = init_val
            else
               tmp(1:nnew) = arr(1:nnew)
            end if
         else
            tmp(1:nnew) = init_val
         end if
         call move_alloc(tmp, arr)
   
         return
      end subroutine base_util_resize_arr_DP
   
   
      subroutine base_util_resize_arr_DPvec(arr, nnew)
         !! author: David A. Minton 
         !! 
         !! Resizes an array component of double precision vectors of size (NDIM, n). Passing nnew = 0 will deallocate. 
         implicit none
         ! Arguments
         real(DP), dimension(:,:), allocatable, intent(inout) :: arr  
            !! Array to resize 
         integer(I4B),                          intent(in)    :: nnew 
            !! New size 
         ! Internals
         real(DP), dimension(:,:), allocatable :: tmp 
            !! Temporary storage array in case the input array is already allocated 
         integer(I4B) :: nold 
            !! Old size 
         real(DP), dimension(NDIM), parameter :: init_val = 0.0_DP
         integer(I4B) :: i
   
         if (nnew < 0) return
   
         if (nnew == 0) then
            if (allocated(arr)) deallocate(arr)
            return
         end if
         
         if (allocated(arr)) then
            nold = size(arr, dim=2)
         else
            nold = 0
         end if
   
         if (nnew == nold) return
         
         allocate(tmp(NDIM, nnew))
         if (nold > 0) then
            if (nnew > nold) then
               tmp(:,1:nold) = arr(:,1:nold)
               do i = nold+1, nnew
                  tmp(:,i) = init_val(:)
               end do
            else
               tmp(:,1:nnew) = arr(:,1:nnew)
            end if
         else
            do i = 1, nnew
               tmp(:, i) = init_val(:)
            end do
         end if
         call move_alloc(tmp, arr)
   
         return
   
         return
      end subroutine base_util_resize_arr_DPvec
   
   
      subroutine base_util_resize_arr_I4B(arr, nnew)
         !! author: David A. Minton 
         !! 
         !! Resizes an array component of integer type. Passing nnew = 0 will deallocate. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr  
            !! Array to resize 
         integer(I4B),                            intent(in)    :: nnew 
            !! New size 
         ! Internals
         integer(I4B), dimension(:), allocatable :: tmp 
            !! Temporary storage array in case the input array is already allocated 
         integer(I4B) :: nold 
            !! Old size 
         integer(I4B), parameter :: init_val = -1
   
         if (nnew < 0) return
   
         if (nnew == 0) then
            if (allocated(arr)) deallocate(arr)
            return
         end if
         
         if (allocated(arr)) then
            nold = size(arr)
         else
            nold = 0
         end if
   
         if (nnew == nold) return
         
         allocate(tmp(nnew))
         if (nold > 0) then
            if (nnew > nold) then
               tmp(1:nold) = arr(1:nold)
               tmp(nold+1:nnew) = init_val
            else
               tmp(1:nnew) = arr(1:nnew)
            end if
         else
            tmp(1:nnew) = init_val
         end if
         call move_alloc(tmp, arr)
   
         return
      end subroutine base_util_resize_arr_I4B
   

      subroutine base_util_resize_arr_logical(arr, nnew)
         !! author: David A. Minton 
         !! 
         !! Resizes an array component of logical type. Passing nnew = 0 will deallocate. 
         implicit none
         ! Arguments
         logical, dimension(:), allocatable, intent(inout) :: arr  
            !! Array to resize 
         integer(I4B),                       intent(in)    :: nnew 
            !! New size 
         ! Internals
         logical, dimension(:), allocatable :: tmp 
            !! Temporary storage array in case the input array is already allocated 
         integer(I4B) :: nold 
            !! Old size 
         logical, parameter :: init_val = .false.
 
         if (nnew < 0) return
   
         if (nnew == 0) then
            if (allocated(arr)) deallocate(arr)
            return
         end if
         
         if (allocated(arr)) then
            nold = size(arr)
         else
            nold = 0
         end if
   
         if (nnew == nold) return
         
         allocate(tmp(nnew))
         if (nold > 0) then
            if (nnew > nold) then
               tmp(1:nold) = arr(1:nold)
               tmp(nold+1:nnew) = init_val
            else
               tmp(1:nnew) = arr(1:nnew)
            end if
         else
            tmp = init_val
         end if
         call move_alloc(tmp, arr)
   
         return
      end subroutine base_util_resize_arr_logical
   
      
      subroutine base_util_resize_storage(self, nnew)
         !! author: David A. Minton 
         !! 
         !! Checks the current size of a Swiftest against the requested size and resizes it if it is too small. 
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self 
            !! Storage object 
         integer(I4B),        intent(in)    :: nnew 
            !! New size 
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
         class(base_storage), intent(inout) :: self 
            !! Storage object 
         integer(I4B),        intent(in)    :: n    
            !! New size 

         if (allocated(self%frame)) deallocate(self%frame)
         allocate(self%frame(n))
         self%nframes = n
   
         return
      end subroutine base_util_setup_storage 
  
      
      subroutine base_util_snapshot_save(self, snapshot)
         !! author: David A. Minton 
         !! 
         !! Checks the current size of the storage object against the required size and extends it by a factor of 2 more than  
         !! requested if it is too small.   
         !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing  
         !! every time you want to add an encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff  
         !! between performance (fewer resize calls) and memory managment. Memory usage grows by a factor of 2 each time it fills  
         !! up, but no more.  
         implicit none
         ! Arguments
         class(base_storage), intent(inout) :: self     
            !! Storage encounter storage object 
         class(*),            intent(in)    :: snapshot 
            !! Object to snapshot 
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


      subroutine base_util_spill_arr_char_string(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of type character strings 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: keeps        
            !! Array of values to keep  
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: discards     
            !! Array of discards 
         logical,               dimension(:),              intent(in)    :: lspill_list  
            !! Logical array of bodies to spill into the discards 
         logical,                                          intent(in)    :: ldestructive 
            !! Logical flag indicating whether or not this operation should alter the keeps array or not 
         ! Internals
         integer(I4B) :: nspill, nkeep, nlist
         character(len=STRMAX), dimension(:), allocatable                :: tmp          
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(nspill))
         else if (size(discards) /= nspill) then
            deallocate(discards)
            allocate(discards(nspill))
         end if

         discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(nkeep))
               tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_char_string
      

      subroutine base_util_spill_arr_DP(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of type DP 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         real(DP), dimension(:), allocatable, intent(inout) :: keeps        
            !! Array of values to keep  
         real(DP), dimension(:), allocatable, intent(inout) :: discards     
            !! Array of discards 
         logical,  dimension(:),              intent(in)    :: lspill_list  
            !! Logical array of bodies to spill into the discardss 
         logical,                             intent(in)    :: ldestructive 
            !! Logical flag indicating whether or not this operation  
                                                                            
            !!    should alter the keeps array or not 
         ! Internals
         integer(I4B) :: nspill, nkeep, nlist
         real(DP), dimension(:), allocatable                :: tmp          
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(nspill))
         else if (size(discards) /= nspill) then
            deallocate(discards)
            allocate(discards(nspill))
         end if

         discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(nkeep))
               tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_DP


      subroutine base_util_spill_arr_DPvec(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of DP vectors with shape (NDIM, n) 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         real(DP), dimension(:,:), allocatable, intent(inout) :: keeps        
            !! Array of values to keep  
         real(DP), dimension(:,:), allocatable, intent(inout) :: discards     
            !! Array discards 
         logical,  dimension(:),                intent(in)    :: lspill_list  
            !! Logical array of bodies to spill into the discards 
         logical,                               intent(in)    :: ldestructive 
            !! Logical flag indicating whether or not this operation should alter the keeps array or not 
         ! Internals
         integer(I4B) :: i, nspill, nkeep, nlist
         real(DP), dimension(:,:), allocatable                :: tmp          
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(NDIM, nspill))
         else if (size(discards, dim=2) /= nspill) then
            deallocate(discards)
            allocate(discards(NDIM, nspill))
         end if

         do i = 1, NDIM
            discards(i,:) = pack(keeps(i,1:nlist), lspill_list(1:nlist))
         end do
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(NDIM, nkeep))
               do i = 1, NDIM
                  tmp(i, :) = pack(keeps(i, 1:nlist), .not. lspill_list(1:nlist))
               end do
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_DPvec


      subroutine base_util_spill_arr_I4B(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of type I4B 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: keeps       
            !! Array of values to keep  
         integer(I4B), dimension(:), allocatable, intent(inout) :: discards    
            !! Array of discards 
         logical,      dimension(:),              intent(in)    :: lspill_list 
            !! Logical array of bodies to spill into the discards 
         logical,                                 intent(in)    :: ldestructive
            !! Logical flag indicating whether or not this  operation should alter the keeps array or not 
         ! Internals
         integer(I4B) :: nspill, nkeep, nlist
         integer(I4B), dimension(:), allocatable                :: tmp         
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(nspill))
         else if (size(discards) /= nspill) then
            deallocate(discards)
            allocate(discards(nspill))
         end if

         discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(nkeep))
               tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_I4B


      subroutine base_util_spill_arr_I8B(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of type I4B 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         integer(I8B), dimension(:), allocatable, intent(inout) :: keeps       
            !! Array of values to keep  
         integer(I8B), dimension(:), allocatable, intent(inout) :: discards    
            !! Array of discards 
         logical,      dimension(:),              intent(in)    :: lspill_list 
            !! Logical array of bodies to spill into the discards 
         logical,                                 intent(in)    :: ldestructive
            !! Logical flag indicating whether or not this  
                                                                               
            !!   operation should alter the keeps array or not 
         ! Internals
         integer(I4B) :: nspill, nkeep, nlist
         integer(I8B), dimension(:), allocatable                :: tmp          
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(nspill))
         else if (size(discards) /= nspill) then
            deallocate(discards)
            allocate(discards(nspill))
         end if

         discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(nkeep))
               tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_I8B


      subroutine base_util_spill_arr_logical(keeps, discards, lspill_list, ldestructive)
         !! author: David A. Minton 
         !! 
         !! Performs a spill operation on a single array of logicals 
         !! This is the inverse of a spill operation 
         implicit none
         ! Arguments
         logical, dimension(:), allocatable, intent(inout) :: keeps        
            !! Array of values to keep  
         logical, dimension(:), allocatable, intent(inout) :: discards     
            !! Array of discards 
         logical, dimension(:),              intent(in)    :: lspill_list  
            !! Logical array of bodies to spill into the discards 
         logical,                            intent(in)    :: ldestructive 
            !! Logical flag indicating whether or not this operation  
                                                                           
            !!    should alter the keeps array or no 
         ! Internals
         integer(I4B) :: nspill, nkeep, nlist
         logical, dimension(:), allocatable                :: tmp          
            !! Array of values to keep  

         nkeep = count(.not.lspill_list(:))
         nspill = count(lspill_list(:))
         nlist = size(lspill_list(:))

         if (.not.allocated(keeps) .or. nspill == 0) return
         if (size(keeps) < nkeep) return
         if (.not.allocated(discards)) then
            allocate(discards(nspill))
         else if (size(discards) /= nspill) then
            deallocate(discards)
            allocate(discards(nspill))
         end if

         discards(:) = pack(keeps(1:nlist), lspill_list(1:nlist))
         if (ldestructive) then
            if (nkeep > 0) then
               allocate(tmp(nkeep))
               tmp(:) = pack(keeps(1:nlist), .not. lspill_list(1:nlist))
               call move_alloc(tmp, keeps)
            else
               deallocate(keeps)
            end if
         end if

         return
      end subroutine base_util_spill_arr_logical


      pure subroutine base_util_sort_dp(arr)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array in place into ascending numerical order using quicksort. 
         implicit none
         ! Arguments
         real(DP), dimension(:), intent(inout) :: arr

         call base_util_sort_qsort_DP(arr)

         return
      end subroutine base_util_sort_dp


      pure subroutine base_util_sort_index_dp(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array by index in ascending numerical order using quick sort. 
         !! This algorithm works well for partially sorted arrays (which is usually the case here). 
         !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously 
         !! sorted array). If it is not allocated, this subroutine swiftest_allocates it. 
         implicit none
         ! Arguments
         real(DP),     dimension(:),              intent(in)    :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
         ! Internals
         integer(I4B) :: n, i
         real(DP), dimension(:), allocatable :: tmparr

         n = size(arr)
         if (.not.allocated(ind)) then
            allocate(ind(n))
            ind = [(i, i=1, n)]
         end if
         allocate(tmparr, mold=arr)
         tmparr(:) = arr(ind(:))
         call base_util_sort_qsort_DP(tmparr, ind)
      
         return
      end subroutine base_util_sort_index_dp


      recursive pure subroutine base_util_sort_qsort_DP(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array by index in ascending numerical order using quicksort sort. 
         implicit none
         ! Arguments
         real(DP), dimension(:), intent(inout)           :: arr
         integer(I4B),dimension(:),intent(out), optional :: ind
         
            !! Internals 
         integer :: iq

         if (size(arr) > 1) then
            if (present(ind)) then
               call base_util_sort_partition_DP(arr, iq, ind)
               call base_util_sort_qsort_DP(arr(:iq-1),ind(:iq-1))
               call base_util_sort_qsort_DP(arr(iq:),  ind(iq:))
            else
               call base_util_sort_partition_DP(arr, iq)
               call base_util_sort_qsort_DP(arr(:iq-1))
               call base_util_sort_qsort_DP(arr(iq:))
            end if
         end if

         return
      end subroutine base_util_sort_qsort_DP

   
      pure subroutine base_util_sort_partition_DP(arr, marker, ind)
         !! author: David A. Minton 
         !! 
         !! Partition function for quicksort on DP type 
         !! 
         implicit none
         ! Arguments
         real(DP),     intent(inout), dimension(:)           :: arr
         integer(I4B), intent(inout), dimension(:), optional :: ind
         integer(I4B), intent(out)                           :: marker
         ! Internals
         integer(I4B) :: i, j, itmp, narr, ipiv
         real(DP) :: temp
         real(DP) :: x   ! pivot point

         narr = size(arr)

         ! Get center as pivot, as this is likely partially sorted
         ipiv = narr / 2
         x = arr(ipiv)
         i = 0
         j = narr + 1
      
         do
            j = j - 1
            do
               if (arr(j) <= x) exit
               j = j - 1
            end do
            i = i + 1
            do
               if (arr(i) >= x) exit
               i = i + 1
            end do
            if (i < j) then
               ! exchange A(i) and A(j)
               temp = arr(i)
               arr(i) = arr(j)
               arr(j) = temp
               if (present(ind)) then
                  itmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = itmp
               end if
            else if (i == j) then
               marker = i + 1
               return
            else
               marker = i
               return
            endif
         end do
   
         return
      end subroutine base_util_sort_partition_DP
   

      pure subroutine base_util_sort_i4b(arr)
         !! author: David A. Minton 
         !! 
         !! Sort input integer array in place into ascending numerical order using quick sort. 
         !! This algorithm works well for partially sorted arrays (which is usually the case here) 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), intent(inout) :: arr

         call base_util_sort_qsort_I4B(arr)

         return
      end subroutine base_util_sort_i4b


      pure subroutine base_util_sort_index_I4B(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input integer array by index in ascending numerical order using quicksort. 
         !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously 
         !! sorted array). If it is not allocated, this subroutine swiftest_allocates it. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:),              intent(in)  :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
         ! Internals
         integer(I4B) :: n, i
         integer(I4B), dimension(:), allocatable :: tmparr

         n = size(arr)
         if (.not.allocated(ind)) then
            allocate(ind(n))
            ind = [(i, i=1, n)]
         end if
         allocate(tmparr, mold=arr)
         tmparr(:) = arr(ind(:))
         call base_util_sort_qsort_I4B(tmparr, ind)

         return
      end subroutine base_util_sort_index_I4B


      pure subroutine base_util_sort_index_I4B_I8Bind(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input integer array by index in ascending numerical order using quicksort. 
         !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously 
         !! sorted array). If it is not allocated, this subroutine swiftest_allocates it. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:),              intent(in)  :: arr
         integer(I8B), dimension(:), allocatable, intent(inout) :: ind
         ! Internals
         integer(I8B) :: n, i
         integer(I4B), dimension(:), allocatable :: tmparr

         n = size(arr)
         if (.not.allocated(ind)) then
            allocate(ind(n))
            ind = [(i, i=1_I8B, n)]
         end if
         allocate(tmparr, mold=arr)
         tmparr(:) = arr(ind(:))
         call base_util_sort_qsort_I4B_I8Bind(tmparr, ind)

         return
      end subroutine base_util_sort_index_I4B_I8Bind


      pure subroutine base_util_sort_index_I8B_I8Bind(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input integer array by index in ascending numerical order using quicksort. 
         !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously 
         !! sorted array). If it is not allocated, this subroutine swiftest_allocates it. 
         implicit none
         ! Arguments
         integer(I8B), dimension(:),              intent(in)  :: arr
         integer(I8B), dimension(:), allocatable, intent(inout) :: ind
         ! Internals
         integer(I8B) :: n, i
         integer(I8B), dimension(:), allocatable :: tmparr

         n = size(arr)
         if (.not.allocated(ind)) then
            allocate(ind(n))
            ind = [(i, i=1_I8B, n)]
         end if
         allocate(tmparr, mold=arr)
         tmparr(:) = arr(ind(:))
         call base_util_sort_qsort_I8B_I8Bind(tmparr, ind)

         return
      end subroutine base_util_sort_index_I8B_I8Bind


      recursive pure subroutine base_util_sort_qsort_I4B(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input I4B array by index in ascending numerical order using quicksort. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), intent(inout)          :: arr
         integer(I4B), dimension(:), intent(out),  optional :: ind
         ! Internals
         integer(I4B) :: iq

         if (size(arr) > 1) then
            if (present(ind)) then
               call base_util_sort_partition_I4B(arr, iq, ind)
               call base_util_sort_qsort_I4B(arr(:iq-1),ind(:iq-1))
               call base_util_sort_qsort_I4B(arr(iq:),  ind(iq:))
            else
               call base_util_sort_partition_I4B(arr, iq)
               call base_util_sort_qsort_I4B(arr(:iq-1))
               call base_util_sort_qsort_I4B(arr(iq:))
            end if
         end if

         return
      end subroutine base_util_sort_qsort_I4B


      recursive pure subroutine base_util_sort_qsort_I4B_I8Bind(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input I4B array by index in ascending numerical order using quicksort. 
         !! 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), intent(inout)          :: arr
         integer(I8B), dimension(:), intent(out),  optional :: ind
         ! Internals
         integer(I8B) :: iq

         if (size(arr) > 1_I8B) then
            if (present(ind)) then
               call base_util_sort_partition_I4B_I8Bind(arr, iq, ind)
               call base_util_sort_qsort_I4B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
               call base_util_sort_qsort_I4B_I8Bind(arr(iq:),  ind(iq:))
            else
               call base_util_sort_partition_I4B_I8Bind(arr, iq)
               call base_util_sort_qsort_I4B_I8Bind(arr(:iq-1_I8B))
               call base_util_sort_qsort_I4B_I8Bind(arr(iq:))
            end if
         end if

         return
      end subroutine base_util_sort_qsort_I4B_I8Bind


      recursive pure subroutine base_util_sort_qsort_I8B_I8Bind(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input I8B array by index in ascending numerical order using quicksort. 
         implicit none
         ! Arguments
         integer(I8B), dimension(:), intent(inout)          :: arr
         integer(I8B), dimension(:), intent(out),  optional :: ind
         ! Internals
         integer(I8B) :: iq

         if (size(arr) > 1_I8B) then
            if (present(ind)) then
               call base_util_sort_partition_I8B_I8Bind(arr, iq, ind)
               call base_util_sort_qsort_I8B_I8Bind(arr(:iq-1_I8B),ind(:iq-1_I8B))
               call base_util_sort_qsort_I8B_I8Bind(arr(iq:),  ind(iq:))
            else
               call base_util_sort_partition_I8B_I8Bind(arr, iq)
               call base_util_sort_qsort_I8B_I8Bind(arr(:iq-1_I8B))
               call base_util_sort_qsort_I8B_I8Bind(arr(iq:))
            end if
         end if

         return
      end subroutine base_util_sort_qsort_I8B_I8Bind

   
      pure subroutine base_util_sort_partition_I4B(arr, marker, ind)
         
            !! author: David A. Minton 
         
            !! 
         
            !! Partition function for quicksort on I4B type 
         
            !! 
         implicit none
         ! Arguments
         integer(I4B), intent(inout), dimension(:)           :: arr
         integer(I4B), intent(inout), dimension(:), optional :: ind
         integer(I4B), intent(out)                           :: marker
         ! Internals
         integer(I4B) :: i, j, itmp, narr, ipiv
         integer(I4B) :: temp
         integer(I4B) :: x   ! pivot point

         narr = size(arr)

         ! Get center as pivot, as this is likely partially sorted
         ipiv = narr / 2
         x = arr(ipiv)
         i = 0
         j = narr + 1
      
         do
            j = j - 1
            do
               if (arr(j) <= x) exit
               j = j - 1
            end do
            i = i + 1
            do
               if (arr(i) >= x) exit
               i = i + 1
            end do
            if (i < j) then
               ! exchange A(i) and A(j)
               temp = arr(i)
               arr(i) = arr(j)
               arr(j) = temp
               if (present(ind)) then
                  itmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = itmp
               end if
            else if (i == j) then
               marker = i + 1
               return
            else
               marker = i
               return
            endif
         end do
   
         return
      end subroutine base_util_sort_partition_I4B


      pure subroutine base_util_sort_partition_I4B_I8Bind(arr, marker, ind)
         !! author: David A. Minton 
         !! 
         !! Partition function for quicksort on I4B type 
         implicit none
         ! Arguments
         integer(I4B), intent(inout), dimension(:)           :: arr
         integer(I8B), intent(inout), dimension(:), optional :: ind
         integer(I8B), intent(out)                           :: marker
         ! Internals
         integer(I8B) :: i, j, itmp, narr, ipiv
         integer(I4B) :: temp
         integer(I8B) :: x   ! pivot point

         narr = size(arr)

         ! Get center as pivot, as this is likely partially sorted
         ipiv = narr / 2_I8B
         x = arr(ipiv)
         i = 0_I8B
         j = narr + 1_I8B
      
         do
            j = j - 1_I8B
            do
               if (arr(j) <= x) exit
               j = j - 1_I8B
            end do
            i = i + 1_I8B
            do
               if (arr(i) >= x) exit
               i = i + 1_I8B
            end do
            if (i < j) then
               ! exchange A(i) and A(j)
               temp = arr(i)
               arr(i) = arr(j)
               arr(j) = temp
               if (present(ind)) then
                  itmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = itmp
               end if
            else if (i == j) then
               marker = i + 1_I8B
               return
            else
               marker = i
               return
            endif
         end do
   
         return
      end subroutine base_util_sort_partition_I4B_I8Bind


      pure subroutine base_util_sort_partition_I8B_I8Bind(arr, marker, ind)
         !! author: David A. Minton 
         !! 
         !! Partition function for quicksort on I8B type with I8B index 
         implicit none
         ! Arguments
         integer(I8B), intent(inout), dimension(:)           :: arr
         integer(I8B), intent(inout), dimension(:), optional :: ind
         integer(I8B), intent(out)                           :: marker
         ! Internals
         integer(I8B) :: i, j, itmp, narr, ipiv
         integer(I8B) :: temp
         integer(I8B) :: x   ! pivot point

         narr = size(arr)

         ! Get center as pivot, as this is likely partially sorted
         ipiv = narr / 2_I8B
         x = arr(ipiv)
         i = 0_I8B
         j = narr + 1_I8B
      
         do
            j = j - 1_I8B
            do
               if (arr(j) <= x) exit
               j = j - 1_I8B
            end do
            i = i + 1_I8B
            do
               if (arr(i) >= x) exit
               i = i + 1_I8B
            end do
            if (i < j) then
               ! exchange A(i) and A(j)
               temp = arr(i)
               arr(i) = arr(j)
               arr(j) = temp
               if (present(ind)) then
                  itmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = itmp
               end if
            else if (i == j) then
               marker = i + 1_I8B
               return
            else
               marker = i
               return
            endif
         end do
   
         return
      end subroutine base_util_sort_partition_I8B_I8Bind


      pure subroutine base_util_sort_sp(arr)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array in place into ascending numerical order using quicksort. 
         !! 
         implicit none
         ! Arguments
         real(SP), dimension(:), intent(inout) :: arr

         call base_util_sort_qsort_SP(arr)

         return
      end subroutine base_util_sort_sp


      pure subroutine base_util_sort_index_sp(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array by index in ascending numerical order using quicksort. 
         !! If ind is supplied already allocated, we assume it is an existing index array (e.g. a previously 
         !! sorted array). If it is not allocated, this subroutine swiftest_allocates it. 
         implicit none
         ! Arguments
         real(SP),     dimension(:),              intent(in)    :: arr
         integer(I4B), dimension(:), allocatable, intent(inout) :: ind
         ! Internals
         integer(I4B) :: n, i
         real(SP), dimension(:), allocatable :: tmparr

         n = size(arr)
         if (.not.allocated(ind)) then
            allocate(ind(n))
            ind = [(i, i=1, n)]
         end if
         allocate(tmparr, mold=arr)
         tmparr(:) = arr(ind(:))
         call base_util_sort_qsort_SP(tmparr, ind)
      
         return
      end subroutine base_util_sort_index_sp


      recursive pure subroutine base_util_sort_qsort_SP(arr, ind)
         !! author: David A. Minton 
         !! 
         !! Sort input DP precision array by index in ascending numerical order using quicksort. 
         implicit none
         ! Arguments
         real(SP), dimension(:), intent(inout)           :: arr
         integer(I4B),dimension(:),intent(out), optional :: ind
         
         ! Internals 
         integer :: iq

         if (size(arr) > 1) then
            if (present(ind)) then
               call base_util_sort_partition_SP(arr, iq, ind)
               call base_util_sort_qsort_SP(arr(:iq-1),ind(:iq-1))
               call base_util_sort_qsort_SP(arr(iq:),  ind(iq:))
            else
               call base_util_sort_partition_SP(arr, iq)
               call base_util_sort_qsort_SP(arr(:iq-1))
               call base_util_sort_qsort_SP(arr(iq:))
            end if
         end if

         return
      end subroutine base_util_sort_qsort_SP


      pure subroutine base_util_sort_partition_SP(arr, marker, ind)
         !! author: David A. Minton 
         !! 
         !! Partition function for quicksort on SP type 
         implicit none
         ! Arguments
         real(SP),     intent(inout), dimension(:)           :: arr
         integer(I4B), intent(inout), dimension(:), optional :: ind
         integer(I4B), intent(out)                           :: marker
         ! Internals
         integer(I4B) :: i, j, itmp, narr, ipiv
         real(SP) :: temp
         real(SP) :: x   ! pivot point

         narr = size(arr)

         ! Get center as pivot, as this is likely partially sorted
         ipiv = narr / 2
         x = arr(ipiv)
         i = 0
         j = narr + 1
      
         do
            j = j - 1
            do
               if (arr(j) <= x) exit
               j = j - 1
            end do
            i = i + 1
            do
               if (arr(i) >= x) exit
               i = i + 1
            end do
            if (i < j) then
               ! exchange A(i) and A(j)
               temp = arr(i)
               arr(i) = arr(j)
               arr(j) = temp
               if (present(ind)) then
                  itmp = ind(i)
                  ind(i) = ind(j)
                  ind(j) = itmp
               end if
            else if (i == j) then
               marker = i + 1
               return
            else
               marker = i
               return
            endif
         end do
   
         return
      end subroutine base_util_sort_partition_SP


      pure subroutine base_util_sort_rearrange_arr_char_string(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of character string in-place from an index list. 
         implicit none
         ! Arguments
         character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I4B),          dimension(:),              intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I4B),                                     intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         character(len=STRMAX), dimension(:), allocatable  :: tmp 
            !! Temporary copy of arry used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_char_string


      pure subroutine base_util_sort_rearrange_arr_DP(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of DP type in-place from an index list. 
         implicit none
         ! Arguments
         real(DP),     dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I4B), dimension(:),              intent(in)  :: ind 
            !! Index to rearrange against 
         integer(I4B),                            intent(in)  :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         real(DP), dimension(:), allocatable :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_DP


      pure subroutine base_util_sort_rearrange_arr_DPvec(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of (NDIM,n) DP-type vectors in-place from an index list. 
         implicit none
         ! Arguments
         real(DP),     dimension(:,:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I4B), dimension(:),                intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I4B),                              intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         real(DP), dimension(:,:), allocatable :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(:,1:n) = arr(:, ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_DPvec


      pure subroutine base_util_sort_rearrange_arr_I4B(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of integers in-place from an index list. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I4B), dimension(:),              intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I4B),                             intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         integer(I4B), dimension(:), allocatable                :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_I4B

      pure subroutine base_util_sort_rearrange_arr_I4B_I8Bind(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of integers in-place from an index list. 
         implicit none
         ! Arguments
         integer(I4B), dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I8B), dimension(:),              intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I8B),                            intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         integer(I4B), dimension(:), allocatable                :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0_I8B) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_I4B_I8Bind


      pure subroutine base_util_sort_rearrange_arr_logical_I8Bind(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of logicals in-place from an index list. 
         implicit none
         ! Arguments
         logical,      dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I8B), dimension(:),              intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I8B),                            intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         logical, dimension(:), allocatable                :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_logical_I8Bind


      pure subroutine base_util_sort_rearrange_arr_logical(arr, ind, n)
         !! author: David A. Minton 
         !! 
         !! Rearrange a single array of logicals in-place from an index list. 
         implicit none
         ! Arguments
         logical,      dimension(:), allocatable, intent(inout) :: arr 
            !! Destination array  
         integer(I4B), dimension(:),              intent(in)    :: ind 
            !! Index to rearrange against 
         integer(I4B),                            intent(in)    :: n   
            !! Number of elements in arr and ind to rearrange 
         ! Internals
         logical, dimension(:), allocatable                :: tmp 
            !! Temporary copy of array used during rearrange operation 

         if (.not. allocated(arr) .or. n <= 0) return
         allocate(tmp, mold=arr)
         tmp(1:n) = arr(ind)
         call move_alloc(tmp, arr)

         return
      end subroutine base_util_sort_rearrange_arr_logical


      subroutine base_util_unique_DP(input_array, output_array, index_map)
         !! author: David A. Minton 
         !! 
         !! Takes an input unsorted integer array and returns a new array of sorted, unique values (DP version) 
         implicit none
         ! Arguments
         real(DP),     dimension(:),              intent(in)  :: input_array  
            !! Unsorted input array  
         real(DP),     dimension(:), allocatable, intent(out) :: output_array 
            !! Sorted array of unique values  
         integer(I4B), dimension(:), allocatable, intent(out) :: index_map    
            !! An array of the same size as input_array that such  that any for any index i,  
            !!    output_array(index_map(i)) = input_array(i) 
         ! Internals
         real(DP), dimension(:), allocatable :: unique_array
         integer(I4B) :: n
         real(DP) :: lo, hi
   
         allocate(unique_array, mold=input_array)
         allocate(index_map(size(input_array)))
         lo = minval(input_array) - 1
         hi = maxval(input_array)
   
         n = 0
         do 
            n = n + 1
            lo = minval(input_array(:), mask=input_array(:) > lo)
            unique_array(n) = lo
            where(abs(input_array(:) - lo) < epsilon(1.0_DP) * lo) index_map(:) = n
            if (lo >= hi) exit
         enddo
         allocate(output_array(n), source=unique_array(1:n)) 
   
         return
      end subroutine base_util_unique_DP
   

      subroutine base_util_unique_I4B(input_array, output_array, index_map)
         !! author: David A. Minton 
         !! 
         !! Takes an input unsorted integer array and returns a new array of sorted, unique values (I4B version)
         implicit none
         ! Arguments
         integer(I4B), dimension(:),              intent(in)  :: input_array  
            !! Unsorted input array  
         integer(I4B), dimension(:), allocatable, intent(out) :: output_array 
            !! Sorted array of unique values 
         integer(I4B), dimension(:), allocatable, intent(out) :: index_map    
            !! An array of the same size as input_array that such  that any for any index i,  
            !! output_array(index_map(i)) = input_array(i)      
         ! Internals
         integer(I4B), dimension(:), allocatable :: unique_array
         integer(I4B) :: n, lo, hi
   
         allocate(unique_array, mold=input_array)
         allocate(index_map, mold=input_array)
         lo = minval(input_array) - 1
         hi = maxval(input_array)
   
         n = 0
         do 
            n = n + 1
            lo = minval(input_array(:), mask=input_array(:) > lo)
            unique_array(n) = lo
            where(input_array(:) == lo) index_map(:) = n
            if (lo >= hi) exit
         enddo
         allocate(output_array(n), source=unique_array(1:n)) 
   
         return
      end subroutine base_util_unique_I4B


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
         class(base_parameters),intent(inout),codimension[*]  :: self  
            !! Collection of parameters 
         ! Internals

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
         call coclone(self%lenc_save_closest)
         call coclone(self%interaction_loops)
         call coclone(self%encounter_check_plpl)
         call coclone(self%encounter_check_pltp)
         call coclone(self%lflatten_interactions)
         call coclone(self%lencounter_sas_plpl)
         call coclone(self%lencounter_sas_pltp)
         call coclone(self%lrhill_present)
         call coclone(self%lextra_force)
         call coclone(self%lbig_discard)
         call coclone(self%lclose)
         call coclone(self%lenergy)
         call coclone(self%lnon_spherical_cb)
         call coclone(self%lrotation)
         call coclone(self%ltides)
         call coclone(self%E_orbit_orig)
         call coclone(self%GMtot_orig)
         call coclonevec(self%L_total_orig)
         call coclonevec(self%L_orbit_orig)
         call coclonevec(self%L_spin_orig)
         call coclonevec(self%L_escape)
         call coclone(self%GMescape)
         call coclone(self%E_collisions)
         call coclone(self%E_untracked)
         call coclone(self%lfirstenergy)
         call coclone(self%lfirstkick)
         call coclone(self%lrestart)
         call coclone(self%display_style)
         call coclone(self%display_unit)
         call coclone(self%log_output )
         call coclone(self%lgr)
         call coclone(self%lyarkovsky)
         call coclone(self%lyorp)
         call coclone(self%seed)
         call coclone(self%lcoarray)

         return
      end subroutine base_coclone_param 
#endif
    

end module base
