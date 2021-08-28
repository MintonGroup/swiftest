submodule (swiftest_classes) s_util_resize
   use swiftest
contains

   module subroutine util_resize_arr_char_string(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                     intent(in)    :: nnew !! New size
      ! Internals
      character(len=STRMAX), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_char_string


   module subroutine util_resize_arr_DP(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision type. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                        intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_DP


   module subroutine util_resize_arr_DPvec(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of double precision vectors of size (NDIM, n). Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                          intent(in)    :: nnew !! New size
      ! Internals
      real(DP), dimension(:,:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr, dim=2)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(NDIM, nnew))
      if (nnew > nold) then
         tmp(:, 1:nold) = arr(:, 1:nold)
      else
         tmp(:, 1:nnew) = arr(:, 1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_DPvec


   module subroutine util_resize_arr_I4B(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of integer type. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                            intent(in)    :: nnew !! New size
      ! Internals
      integer(I4B), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_I4B


   module subroutine util_resize_arr_info(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(swiftest_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                         intent(in)    :: nnew !! New size
      ! Internals
      type(swiftest_particle_info), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size
      logical :: is_symba

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         call util_copy_particle_info_arr(arr(1:nold), tmp(1:nold))
      else
         call util_copy_particle_info_arr(arr(1:nnew), tmp(1:nnew))
      end if

      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_info


   module subroutine util_resize_arr_logical(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of logical type. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                       intent(in)    :: nnew !! New size
      ! Internals
      logical, dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
      integer(I4B) :: nold !! Old size

      if (.not. allocated(arr) .or. nnew < 0) return

      nold = size(arr)
      if (nnew == nold) return

      if (nnew == 0) then
         deallocate(arr)
         return
      end if
      
      allocate(tmp(nnew))
      if (nnew > nold) then
         tmp(1:nold) = arr(1:nold)
      else
         tmp(1:nnew) = arr(1:nnew)
      end if
      call move_alloc(tmp, arr)

      return
   end subroutine util_resize_arr_logical


   module subroutine util_resize_body(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  !! Swiftest body object
      integer(I4B),               intent(in)    :: nnew  !! New size neded

      call util_resize(self%info, nnew)
      call util_resize(self%id, nnew)
      call util_resize(self%status, nnew)
      call util_resize(self%ldiscard, nnew)
      call util_resize(self%lmask, nnew)
      call util_resize(self%mu, nnew)
      call util_resize(self%xh, nnew)
      call util_resize(self%vh, nnew)
      call util_resize(self%xb, nnew)
      call util_resize(self%vb, nnew)
      call util_resize(self%ah, nnew)
      call util_resize(self%aobl, nnew)
      call util_resize(self%atide, nnew)
      call util_resize(self%agr, nnew)
      call util_resize(self%ir3h, nnew)
      call util_resize(self%a, nnew)
      call util_resize(self%e, nnew)
      call util_resize(self%inc, nnew)
      call util_resize(self%capom, nnew)
      call util_resize(self%omega, nnew)
      call util_resize(self%capm, nnew)
      self%nbody = count(self%status(1:nnew) /= INACTIVE)

      return
   end subroutine util_resize_body


   module subroutine util_resize_encounter(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      class(swiftest_encounter), intent(inout) :: self !! Swiftest encounter list 
      integer(I4B),              intent(in)    :: nnew !! New size of list needed
      ! Internals
      class(swiftest_encounter), allocatable :: enc_temp
      integer(I4B)                           :: nold
      logical                                :: lmalloc

      lmalloc = allocated(self%status)
      if (lmalloc) then
         nold = size(self%status)
      else
         nold = 0
      end if
      if (nnew > nold) then
         if (lmalloc) allocate(enc_temp, source=self)
         call self%setup(2 * nnew)
         if (lmalloc) then
            call self%copy(enc_temp)
            deallocate(enc_temp)
         end if
      else
         self%status(nnew+1:nold) = INACTIVE
      end if
      self%nenc = nnew

      return
   end subroutine util_resize_encounter


   module subroutine util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_pl), intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call util_resize_body(self, nnew)

      call util_resize(self%mass, nnew)
      call util_resize(self%Gmass, nnew)
      call util_resize(self%rhill, nnew)
      call util_resize(self%radius, nnew)
      call util_resize(self%xbeg, nnew)
      call util_resize(self%xend, nnew)
      call util_resize(self%vbeg, nnew)
      call util_resize(self%density, nnew)
      call util_resize(self%Ip, nnew)
      call util_resize(self%rot, nnew)
      call util_resize(self%k2, nnew)
      call util_resize(self%Q, nnew)
      call util_resize(self%tlag, nnew)

      return
   end subroutine util_resize_pl


   module subroutine util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a Swiftest body against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(swiftest_tp), intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),       intent(in)    :: nnew  !! New size neded

      call util_resize_body(self, nnew)

      call util_resize(self%isperi, nnew)
      call util_resize(self%peri, nnew)
      call util_resize(self%atp, nnew)

      return
   end subroutine util_resize_tp


end submodule s_util_resize