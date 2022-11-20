module io_progress_bar
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use swiftest_globals
   use swiftest_classes
   implicit none
   public

   character(len=1), dimension(4), parameter, private :: spinstr = ["\","-","|","/"] !! The progress spinner sequency
   character(len=1),parameter, private                :: barchar = "#" !! The progress bar character

   type :: progress_bar
      !! author: David A. Minton
      !! 
      !! Implements a class for a simple progress bar that can print on the screen.
      integer(I4B)                  :: PBARSIZE = 80 !! Number of characters acros for a whole progress bar
      integer(I8B)                  :: loop_length   !! The total number of loops that the progrees bar is executing
      integer(I8B)                  :: spinner_cycle_length  !! The number of loop iterations it takes to execute a full spinner cycle
      integer(I8B)                  :: spinner_skip  !! The number of loop iterations to execute before updating the spinner
      character(len=:), allocatable :: barstr        !! The string that prints out as the progress bar
      integer(I4B)                  :: bar_pos       !! The current position of the progress bar
      integer(I4B)                  :: spin_pos      !! Position of the "spinner" that indicates that progress is being made
      character(len=32)             :: fmt           !! The format string that is used to define the progress bar itself
      character(len=64)             :: message       !! The current message displayed at the end of the progress bar
   contains
      procedure :: reset => io_pbar_reset   !! Resets the progress bar to the beginning
      procedure :: update => io_pbar_update !! Updates the progress bar with new values and causes the "spinner" to flip.
   end type progress_bar

contains

   subroutine io_pbar_reset(self, loop_length, spinner_cycle_length)
      !! author: David A. Minton
      !! 
      !! Resets the progress bar to the beginning
      implicit none
      ! Arguments
      class(progress_bar),intent(inout)        :: self         !! The progress bar object
      integer(I8B),       intent(in)           :: loop_length  !! The length of the loop that the progress bar is attached to
      integer(I8B),       intent(in), optional :: spinner_cycle_length !! The number of iterations between a complete spinner cycle. If not passd, the default is to update the spinner each call
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
      self%fmt  = '(A1,"[",A' // numchar // ',"] ",A,$)'
      self%loop_length = loop_length
      self%bar_pos = 0
      self%message = ""
      if (present(spinner_cycle_length)) then
         self%spinner_cycle_length = spinner_cycle_length
         self%spinner_skip = self%spinner_cycle_length / size(spinstr) 
      else
         self%spinner_cycle_length = size(spinstr)
         self%spinner_skip = 1
      end if
      self%spin_pos = 1

      write(*,fmt=self%fmt) char(13),self%barstr,trim(adjustl(self%message))

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
      integer(I4B) :: bar_pos, spin_pos  !! The current integer position of the progress bar and spinner
      logical :: update = .false.

      ! Compute the current position
      frac = real(i,kind=DP) / real(self%loop_length,kind=DP)
      bar_pos = min(int(ceiling(frac * self%PBARSIZE),kind=I4B),self%PBARSIZE)

      if (bar_pos /= self%bar_pos) then
         ! Fill in the bar character up to the current position
         self%barstr(bar_pos:bar_pos) = barchar
         update = .true.
         self%bar_pos = bar_pos
      end if

      ! Compute the current value of the spinner and set the spinner character
      spin_pos = mod(i / self%spinner_skip, size(spinstr)) + 1
      if (spin_pos /= self%spin_pos) then
         self%barstr(bar_pos+1:bar_pos+1) = spinstr(spin_pos)
         self%spin_pos = spin_pos
         update = .true.
      end if

      if (present(message)) then
         if (message /= self%message) then
            update = .true.
            self%message = message 
         end if
      end if

      if (update) write(*,fmt=self%fmt) char(13),self%barstr,trim(adjustl(self%message))


      return
   end subroutine io_pbar_update


end module io_progress_bar
