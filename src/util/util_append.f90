submodule (swiftest_classes) s_util_append
   use swiftest
contains

   module subroutine util_append_arr_char_string(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr    !! Destination array 
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: source !! Array to append 
      logical,               dimension(:), optional,    intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask)
      else
         nsrc = size(source)
      end if

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr+nsrc)

      arr(narr+1:nsrc) = source(:)

      return
   end subroutine util_append_arr_char_string


   module subroutine util_append_arr_DP(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr    !! Destination array 
      real(DP), dimension(:), allocatable, intent(inout) :: source !! Array to append 
      logical,  dimension(:), optional,    intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask)
      else
         nsrc = size(source)
      end if

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr+nsrc)

      arr(narr+1:nsrc) = source(:)

      return
   end subroutine util_append_arr_DP


   module subroutine util_append_arr_DPvec(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr    !! Destination array 
      real(DP), dimension(:,:), allocatable, intent(inout) :: source !! Array to append 
      logical,  dimension(:),   optional,    intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask)
      else
         nsrc = size(source, dim=2)
      end if

      if (allocated(arr)) then
         narr = size(arr, dim=2)
      else
         allocate(arr(NDIM,nsrc))
         narr = 0
      end if

      call util_resize(arr, narr+nsrc)

      arr(:,narr+1:nsrc) = source(:,:)

      return
   end subroutine util_append_arr_DPvec


   module subroutine util_append_arr_I4B(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr    !! Destination array 
      integer(I4B), dimension(:), allocatable, intent(inout) :: source !! Array to append 
      logical,      dimension(:), optional,    intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask)
      else
         nsrc = size(source)
      end if

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr+nsrc)

      arr(narr+1:nsrc) = source(:)

      return
   end subroutine util_append_arr_I4B


   module subroutine util_append_arr_logical(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr    !! Destination array 
      logical, dimension(:), allocatable, intent(inout) :: source !! Array to append 
      logical, dimension(:), optional,    intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (present(lsource_mask)) then
         nsrc = count(lsource_mask)
      else
         nsrc = size(source)
      end if

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr+nsrc)

      arr(narr+1:nsrc) = source(:)

      return
   end subroutine util_append_arr_logical


   module subroutine util_append_body(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self   !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source !! Source object to append
      logical, dimension(:), optional, intent(in)    :: lsource_mask  !! Logical mask indicating which elements to append to

      associate(nold => self%nbody, nnew => source%nbody)

      end associate
      return
   end subroutine util_append_body

end submodule s_util_append