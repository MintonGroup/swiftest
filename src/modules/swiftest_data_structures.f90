module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   use user
   implicit none


   type,public :: swiftest_tp
      integer(I4B)                                :: nbody  !! Number of bodies
      integer(I4B), dimension(:),     allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      integer(I4B), dimension(:),     allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),     allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),     allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
   contains
      procedure :: alloc => swiftest_tp_allocate
      procedure :: dealloc => swiftest_tp_deallocate
      procedure :: read_from_file => swiftest_read_tp_in 
   end type swiftest_tp

   type,public,extends(swiftest_tp) :: swiftest_pl
      real(DP),     dimension(:),     allocatable :: mass   !! Mass
      real(DP),     dimension(:),     allocatable :: radius !! Radius
      real(DP),     dimension(:),     allocatable :: rhill  !! Hill's sphere radius
   contains
      procedure :: alloc => swiftest_pl_allocate
      procedure :: dealloc => swiftest_pl_deallocate
      procedure :: read_from_file => swiftest_read_pl_in 
   end type swiftest_pl

   interface

      module subroutine swiftest_read_pl_in(self,param) 
         class(swiftest_pl), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(user_input_parameters),intent(in) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_pl_in

      module subroutine swiftest_read_tp_in(self,param) 
         class(swiftest_tp), intent(inout)  :: self  !! Swiftest data structure to store massive body initial conditions
         type(user_input_parameters),intent(in) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_tp_in
   end interface


   contains

      subroutine swiftest_tp_allocate(self,n)
         implicit none

         class(swiftest_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate

         self%nbody = n
         if (n <= 0) return
         allocate(self%name(n))
         allocate(self%status(n))
         allocate(self%peri(n))
         allocate(self%atp(n))
         allocate(self%isperi(n))
         allocate(self%xh(NDIM,n))
         allocate(self%vh(NDIM,n))
         allocate(self%xb(NDIM,n))
         allocate(self%vb(NDIM,n))

         self%name = 0
         self%status = 0
         self%peri = 0.0_DP
         self%atp = 0.0_DP
         self%isperi = 0.0_DP
         self%xh = 0.0_DP
         self%vh = 0.0_DP
         self%xb = 0.0_DP
         self%vb = 0.0_DP
         return
      end subroutine swiftest_tp_allocate

      subroutine swiftest_pl_allocate(self,n)
         implicit none

         class(swiftest_pl), intent(inout)    :: self !! Swiftest massive body object
         integer, intent(in)                  :: n    !! Number of massive bodies to allocate

         self%nbody = n
         if (n <= 0) return
         call self%swiftest_tp%alloc(n)

         allocate(self%mass(n))
         allocate(self%radius(n))
         allocate(self%rhill(n))

         self%mass = 0.0_DP
         self%radius = 0.0_DP
         self%rhill = 0.0_DP
         return
      end subroutine swiftest_pl_allocate

      subroutine swiftest_tp_deallocate(self)
         implicit none

         class(swiftest_tp), intent(inout)    :: self

         if (allocated(self%name)) deallocate(self%name)
         if (allocated(self%status)) deallocate(self%status)
         if (allocated(self%isperi)) deallocate(self%isperi)
         if (allocated(self%peri)) deallocate(self%peri)
         if (allocated(self%atp)) deallocate(self%atp)
         if (allocated(self%xh)) deallocate(self%xh)
         if (allocated(self%vh)) deallocate(self%vh)
         if (allocated(self%xb)) deallocate(self%xb)
         if (allocated(self%vb)) deallocate(self%vb)
      return
      end subroutine swiftest_tp_deallocate

      subroutine swiftest_pl_deallocate(self)
         implicit none

         class(swiftest_pl), intent(inout)    :: self
         call self%swiftest_tp%dealloc()

         if (allocated(self%mass)) deallocate(self%mass)
         if (allocated(self%radius)) deallocate(self%radius)
         if (allocated(self%rhill)) deallocate(self%rhill)
      return
      end subroutine swiftest_pl_deallocate


end module swiftest_data_structures
