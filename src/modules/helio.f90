module helio
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter modules: helio.f90
   use swiftest_globals
   use swiftest_data_structures
   implicit none

   !! Helio test particle class
   type, public, extends(swiftest_tp) :: helio_tp
      real(DP), dimension(:,:), allocatable :: ah  !! Total heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ahi !! Teliocentric acceleration due to interactions
   contains
      procedure :: alloc => helio_tp_allocate
      final :: helio_tp_deallocate
      procedure, public :: get_acch => helio_getacch_tp !! Compute heliocentric accelerations of test particles
   end type helio_tp

   !! Helio massive body particle class
   type, public, extends(swiftest_pl) :: helio_pl
      real(DP), dimension(:,:), allocatable :: ah  !! Total heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ahi !! Heliocentric acceleration due to interactions
   contains
      procedure :: alloc => helio_pl_allocate   !! Constructor method - Allocates space for number of particles
      final :: helio_pl_deallocate              !! Finalizer method - Deallocates all allocatables 
      procedure, public :: get_acch => helio_getacch !! Compute heliocentric accelerations of plAnetss
   end type helio_pl



!> Interfaces for all helio particle methods that are implemented in separate submodules 
interface
   subroutine helio_drift(helio_plA, dt)
      implicit none
      class(helio_pl), intent(inout) :: helio_plA
      real(DP), intent(in)             :: dt
   end subroutine helio_drift

   subroutine helio_drift_tp(helio_tpA, mu, dt)
      implicit none
      class(helio_tp), intent(inout) :: helio_tpA
      real(DP), intent(in)             :: mu, dt
   end subroutine helio_drift_tp

   subroutine helio_getacch(helio_plA, param, t, lflag)
      implicit none
      class(helio_pl), intent(inout)  :: helio_plA
      type(user_input_parameters),intent(in) :: param   !! Input collection of user-defined parameter
      real(DP), intent(in)           :: t
      logical, intent(in)       :: lflag
   end subroutine helio_getacch

   subroutine helio_getacch_int(helio_plA)
      implicit none
      class(helio_pl), intent(inout)  :: helio_plA
      integer(I4B), intent(in)       :: npl
   end subroutine helio_getacch_int

   subroutine helio_getacch_int_tp(npl, ntp, helio_plA, helio_tpA, xh)
      implicit none
      integer(I4B), intent(in)                   :: npl, ntp
      real(DP), dimension(NDIM, npl), intent(in) :: xh
      type(helio_pl), intent(inout)           :: helio_plA
      type(helio_tp), intent(inout)              :: helio_tpA
   end subroutine helio_getacch_int_tp

   subroutine helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xh, j2rp2, j4rp4)
      implicit none
      logical, intent(in)                   :: lflag, lextra_force
      integer(I4B), intent(in)                   :: npl, nplmax, ntp, ntpmax
      real(DP), intent(in)                       :: t, j2rp2, j4rp4
      real(DP), dimension(NDIM, npl), intent(in) :: xh
      type(helio_pl), intent(inout)              :: helio_plA
      type(helio_tp), intent(inout)              :: helio_tpA
   end subroutine helio_getacch_tp

   subroutine helio_kickvb(npl, helio_plA, dt)
      implicit none
      integer(I4B), intent(in)       :: npl
      real(DP), intent(in)           :: dt
      type(helio_pl), intent(inout)  :: helio_plA
   end subroutine helio_kickvb

   subroutine helio_kickvb_tp(ntp, helio_tpA, dt)
      implicit none
      integer(I4B), intent(in)       :: ntp
      real(DP), intent(in)           :: dt
      type(helio_tp), intent(inout)  :: helio_tpA
   end subroutine helio_kickvb_tp

   subroutine helio_lindrift(npl, swiftest_plA, dt, pt)
      implicit none
      integer(I4B), intent(in)               :: npl
      real(DP), intent(in)                   :: dt
      real(DP), dimension(NDIM), intent(out) :: pt
      type(swiftest_pl), intent(inout)       :: swiftest_plA
   end subroutine helio_lindrift

   subroutine helio_lindrift_tp(ntp, swiftest_tpA, dt, pt)
      implicit none
      integer(I4B), intent(in)              :: ntp
      real(DP), intent(in)                  :: dt
      real(DP), dimension(NDIM), intent(in) :: pt
      type(swiftest_tp), intent(inout)      :: swiftest_tpA
   end subroutine helio_lindrift_tp

   subroutine helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt)
      implicit none
      logical, intent(in)      :: lextra_force
      logical, intent(inout)   :: lfirst
      integer(I4B), intent(in)      :: npl, nplmax, ntp, ntpmax
      real(DP), intent(in)          :: t, j2rp2, j4rp4, dt
      type(helio_pl), intent(inout) :: helio_plA
      type(helio_tp), intent(inout) :: helio_tpA
   end subroutine helio_step

   subroutine helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
      implicit none
      logical, intent(in)                    :: lextra_force
      logical, intent(inout)                 :: lfirst
      integer(I4B), intent(in)                    :: npl, nplmax
      real(DP), intent(in)                        :: t, j2rp2, j4rp4, dt
      real(DP), dimension(NDIM), intent(out)      :: ptb, pte
      real(DP), dimension(NDIM, npl), intent(out) :: xbeg,xend
      type(helio_pl), intent(inout)               :: helio_plA
   end subroutine helio_step_pl

   subroutine helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
      implicit none
      logical, intent(in)                   :: lextra_force
      logical, intent(inout)                :: lfirsttp
      integer(I4B), intent(in)                   :: npl, nplmax, ntp, ntpmax
      real(DP), intent(in)                       :: t, j2rp2, j4rp4, dt
      real(DP), dimension(NDIM), intent(in)      :: ptb, pte
      real(DP), dimension(NDIM, npl), intent(in) :: xbeg, xend
      type(helio_pl), intent(inout)              :: helio_plA
      type(helio_tp), intent(inout)              :: helio_tpA
   end subroutine helio_step_tp

   subroutine helio_user_getacch(t, npl, helio_plA)
      implicit none
      integer(I4B), intent(in)      :: npl
      real(DP), intent(in)          :: t
      type(helio_pl), intent(inout) :: helio_plA
   end subroutine helio_user_getacch

   subroutine helio_user_getacch_tp(t, ntp, helio_tpA)
      implicit none
      integer(I4B), intent(in)      :: ntp
      real(DP), intent(in)          :: t
      type(helio_tp), intent(inout) :: helio_tpA
   end subroutine helio_user_getacch_tp
end interface

!> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the helio submodules.
   contains
      !! Helio constructor and desctructor methods
      subroutine helio_tp_allocate(self,n)
         !! Helio test particle constructor method
         implicit none

         class(helio_tp), intent(inout) :: self !! Swiftest test particle object
         integer, intent(in)            :: n    !! Number of test particles to allocate

         if (self%is_allocated) then
            write(*,*) 'Helio test particle structure already alllocated'
            return
         end if
         write(*,*) 'Allocating the Swiftest test particle structure'
         
         call self%swiftest_tp%alloc(n)
         if (n <= 0) return
         allocate(self%ah(NDIM,n))
         allocate(self%ahi(NDIM,n))

         self%ah(:,:) = 0.0_DP
         self%ahi(:,:) = 0.0_DP
         return
      end subroutine helio_tp_allocate

      subroutine helio_tp_deallocate(self)
         !! Helio test particle destructor/finalizer
         implicit none

         type(helio_tp), intent(inout)    :: self

         if (self%is_allocated) then
            deallocate(self%ah)
            deallocate(self%ahi)
         end if
         return
      end subroutine helio_tp_deallocate

      subroutine helio_pl_allocate(self,n)
         !! Helio massive body constructor method
         implicit none

         class(helio_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate

         if (self%is_allocated) then
            write(*,*) 'Helio plAnet structure already alllocated'
            return
         end if
         call self%swiftest_pl%alloc(n)

         if (n <= 0) return
         allocate(self%ah(NDIM,n))
         allocate(self%ahi(NDIM,n))

         self%ah(:,:) = 0.0_DP
         self%ahi(:,:) = 0.0_DP
         return
      end subroutine helio_pl_allocate

      subroutine helio_pl_deallocate(self)
         !! Helio massive body destructor/finalizer
         implicit none

         type(helio_pl), intent(inout)    :: self
         if (self%is_allocated) then
            deallocate(self%ah)
            deallocate(self%ahi)
         end if
         return
      end subroutine helio_pl_deallocate

      

end module helio