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
      real(DP),     dimension(:,:),   allocatable :: ah  ! total heliocentric acceleration
      real(DP),     dimension(:,:),   allocatable :: ahi ! heliocentric acceleration due to interactions
   contains
      procedure :: alloc => helio_tp_allocate
      final :: helio_tp_deallocate
   end type helio_tp

   !! Helio massive body particle class
   type, public, extends(swiftest_pl) :: helio_pl
      real(DP),     dimension(:,:),   allocatable :: ah  ! total heliocentric acceleration
      real(DP),     dimension(:,:),   allocatable :: ahi ! heliocentric acceleration due to interactions
   contains
      procedure :: alloc => helio_pl_allocate
      final :: helio_pl_deallocate
   end type helio_pl



!> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the helio submodules.
!interface
!! Interfaces for all helio particle methods that are implemented in separate submodules 
!end interface
   contains
      !! Helio constructor and desctructor methods
      subroutine helio_tp_allocate(self,n)
         !! Helio test particle constructor method
         implicit none

         class(helio_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)                  :: n    !! Number of test particles to allocate

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
            write(*,*) 'Helio planet structure already alllocated'
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