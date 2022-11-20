module io_progress_bar
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use swiftest_globals
   use swiftest_classes
   implicit none
   public

   type :: progress_bar
      !! author: David A. Minton
      !! 
      !! Implements a class for a simple progress bar that can print on the screen.
      integer(I4B)                  :: PBARSIZE = 80 !! Number of characters acros for a whole progress bar
      integer(I8B)                  :: nloops        !! The total number of loops that the progrees bar is executing
      character(len=:), allocatable :: barstr        !! The string that prints out as the progress bar
      integer(I4B)                  :: spinner       !! Position of the "spinner" that indicates that progress is being made
      character(len=1)              :: barchar = "=" !! The progress bar character
      character(len=32)             :: fmt           !! The format string that is used to define the progress bar itself
      integer(I4B)                  :: pos           !! The current position of the progress bar
      character(len=32)             :: message       !! The current message displayed at the end of the progress bar
   contains
      procedure :: reset => io_pbar_reset   !! Resets the progress bar to the beginning
      procedure :: update => io_pbar_update !! Updates the progress bar with new values and causes the "spinner" to flip.
   end type progress_bar

contains

   subroutine io_pbar_reset(self, nloops)
      !! author: David A. Minton
      !! 
      !! Resets the progress bar to the beginning
      implicit none
      ! Arguments
      class(progress_bar),intent(inout) :: self
      integer(I8B),       intent(in)    :: nloops
      ! Internals
      character(len=2) :: numchar
      integer(I4B) :: k

      if (.not.allocated(self%barstr)) then
         allocate(character(self%PBARSIZE) :: self%barstr)
      end if
      do k = 1, self%PBARSIZE
         self%barstr(k:k) = " "
      end do
      write(numchar,'(I2)') self%PBARSIZE
      self%fmt  = '(A1,"[",A' // numchar // ',"]",$)'
      self%nloops = nloops
      self%spinner = 0
      self%pos = 0
      self%message = ""

      write(*,fmt=self%fmt) char(13),self%barstr

      return
   end subroutine io_pbar_reset


   subroutine io_pbar_update(self,i,message)
      !! author: David A. Minton
      !! 
      !! Updates the progress bar with new values and causes the "spinner" to flip.
      implicit none
      ! Arguments
      class(progress_bar), intent(inout)        :: self    !! Progres bar object
      integer(I8B),        intent(in)           :: i       !! The current loop index of the progress loop
      character(len=*),    intent(in), optional :: message !! An optional message to display to the right of the progress bar
      ! Internals
      real(DP)     :: frac
      integer(I4B) :: pos           !! The current integer position of the progress bar
      character(len=1), dimension(4), parameter :: spinstr = ["/","-","\","|"]



      ! Compute the current position
      frac = real(i,kind=DP) / real(self%nloops,kind=DP)
      pos = min(int(ceiling(frac * self%PBARSIZE),kind=I4B),self%PBARSIZE)

      if (pos /= self%pos) then
         self%pos = pos
         ! Fill in the bar character up to the current position
         self%barstr(pos:pos) = self%barchar
      end if

      ! Compute the current value of the spinner and set the spinner character
      self%spinner = self%spinner + 1
      if (self%spinner > size(spinstr)) self%spinner = 1

      self%barstr(pos+1:pos+1) = spinstr(self%spinner)

      write(*,fmt=self%fmt) char(13),self%barstr


      return
   end subroutine io_pbar_update


end module io_progress_bar
