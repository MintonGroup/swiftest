module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   implicit none
   private

   !! An abstract superclass for a generic Swiftest particle. All particle types are derived from this class
   type, abstract, public :: swiftest_particle           
      !! Superclass that defines the generic elements of a Swiftest particle 
      !private
      integer(I4B)                              :: nbody  !! Number of bodies
      integer(I4B), dimension(:),   allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),   allocatable :: status !! Status
      logical                                   :: is_allocated = .false. !! Flag to indicate whether or not the components are allocated
   contains
      procedure, public :: superalloc => swiftest_particle_allocate  !! A base constructor that sets nbody and allocates the common components
      !procedure, public :: get_nbody => swiftest_get_nbody
      !procedure, public :: get_name => swiftest_get_name
      procedure(swiftest_read_particle_input_file), public, deferred :: set_from_file
      final :: swiftest_particle_dealloc
   end type swiftest_particle

   abstract interface
      !! An abstract interface that defines the generic file read method for all Swiftest particle types
      subroutine swiftest_read_particle_input_file(self,param) 
         use user
         implicit none
         class(swiftest_particle), intent(inout)  :: self  !! Generic swiftest body initial conditions
         type(user_input_parameters),intent(in) :: param   !! Input collection of user-defined parameters
      end subroutine swiftest_read_particle_input_file
   end interface

   !! Basic Swiftest test particle class
   type, public, extends(swiftest_particle) :: swiftest_tp 
      !private
      integer(I4B), dimension(:),   allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),   allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),   allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:), allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:), allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:), allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:), allocatable :: vb     !! Barycentric velocity
   contains
      procedure, public :: alloc => swiftest_tp_allocate
      procedure, public :: set_from_file => swiftest_read_tp_in 
      final :: swiftest_tp_deallocate
   end type swiftest_tp

   !! Basic Swiftest massive body particle class
   type,public,extends(swiftest_tp) :: swiftest_pl
      !private
      real(DP), dimension(:), allocatable :: mass   !! Mass
      real(DP), dimension(:), allocatable :: radius !! Radius
      real(DP), dimension(:), allocatable :: rhill  !! Hill's sphere radius
   contains
      procedure, public :: alloc => swiftest_pl_allocate
      procedure, public :: set_from_file => swiftest_read_pl_in 
      final :: swiftest_pl_deallocate
   end type swiftest_pl

   !> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the swiftest submodules
   interface
      !! Interfaces for all Swiftest particle methods that are implemented in separate submodules 
      module subroutine swiftest_read_pl_in(self,param) 
         use user
         implicit none
         class(swiftest_pl), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(user_input_parameters),intent(in) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_pl_in

      module subroutine swiftest_read_tp_in(self,param) 
         use user
         implicit none
         class(swiftest_tp), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(user_input_parameters),intent(in) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_tp_in
   end interface

   contains
      !! Swiftest constructor and desctructor methods
      subroutine swiftest_particle_allocate(self,n)
         !! Basic Swiftest generic particle constructor method
         implicit none

         class(swiftest_particle), intent(inout)   :: self !! Swiftest generic particle object
         integer, intent(in)                       :: n    !! Number of particles to allocate

         self%nbody = n
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest particle structure already alllocated'
            return
         end if
         write(*,*) 'Allocating the basic Swiftest particle'
         allocate(self%name(n))
         allocate(self%status(n))

         self%name(:) = 0
         self%status(:) = 0
         self%is_allocated = .true.

         return
      end subroutine swiftest_particle_allocate

      subroutine swiftest_particle_deallocate(self)
         !! Basic Swiftest generic particle destructor/finalizer
         implicit none

         type(swiftest_particle), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%name)
            deallocate(self%status)
            self%is_allocated = .false.
         end if
         return
      end subroutine swiftest_particle_deallocate

      subroutine swiftest_tp_allocate(self,n)
         !! Basic Swiftest test particle constructor method
         implicit none

         class(swiftest_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate
        
         call self%superalloc(n)
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest test particle structure already alllocated'
            return
         end if
         write(*,*) 'Allocating the Swiftest test particle'

         allocate(self%peri(n))
         allocate(self%atp(n))
         allocate(self%isperi(n))
         allocate(self%xh(NDIM,n))
         allocate(self%vh(NDIM,n))
         allocate(self%xb(NDIM,n))
         allocate(self%vb(NDIM,n))

         self%peri(:) = 0.0_DP
         self%atp(:) = 0.0_DP
         self%isperi(:) = 0.0_DP
         self%xh(:,:) = 0.0_DP
         self%vh(:,:) = 0.0_DP
         self%xb(:,:) = 0.0_DP
         self%vb(:,:) = 0.0_DP

         return
      end subroutine swiftest_tp_allocate

      subroutine swiftest_tp_deallocate(self)
         !! Basic Swiftest test particle destructor/finalizer
         implicit none

         type(swiftest_tp), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%isperi)
            deallocate(self%peri)
            deallocate(self%atp)
            deallocate(self%xh)
            deallocate(self%vh)
            deallocate(self%xb)
            deallocate(self%vb)
         end if
         return
      end subroutine swiftest_tp_deallocate

      subroutine swiftest_pl_allocate(self,n)
         !! Basic Swiftest massive body constructor method
         implicit none

         class(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
         integer, intent(in)                  :: n    !! Number of massive bodies to allocate

         call self%superalloc(n)
         if (n <= 0) return

         if (self%is_allocated) then
            write(*,*) 'Swiftest massive body structure already alllocated'
            return
         end if

         call self%swiftest_tp%alloc(n)

         allocate(self%mass(n))
         allocate(self%radius(n))
         allocate(self%rhill(n))

         self%mass(:) = 0.0_DP
         self%radius(:) = 0.0_DP
         self%rhill(:) = 0.0_DP
         return
      end subroutine swiftest_pl_allocate

      subroutine swiftest_pl_deallocate(self)
         !! Basic Swiftest massive body destructor/finalizer
         implicit none

         type(swiftest_pl), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%mass)
            deallocate(self%radius)
            deallocate(self%rhill)
         end if
         return
      end subroutine swiftest_pl_deallocate

end module swiftest_data_structures
