submodule (swiftest_classes) s_util_append
   use swiftest
contains

   module subroutine util_append_arr_char_string(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical,               dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      nsrc = count(lsource_mask)

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr + nsrc)

      arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))

      return
   end subroutine util_append_arr_char_string


   module subroutine util_append_arr_DP(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical,  dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      nsrc = count(lsource_mask)

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr + nsrc)

      arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))

      return
   end subroutine util_append_arr_DP


   module subroutine util_append_arr_DPvec(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
      logical,  dimension(:),                intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      nsrc = count(lsource_mask)

      if (allocated(arr)) then
         narr = size(arr, dim=2)
      else
         allocate(arr(NDIM, nsrc))
         narr = 0
      end if

      call util_resize(arr, narr + nsrc)

      arr(1, narr + 1:narr + nsrc) = pack(source(1,:), lsource_mask(:))
      arr(2, narr + 1:narr + nsrc) = pack(source(2,:), lsource_mask(:))
      arr(3, narr + 1:narr + nsrc) = pack(source(3,:), lsource_mask(:))

      return
   end subroutine util_append_arr_DPvec


   module subroutine util_append_arr_I4B(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical,      dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      nsrc = count(lsource_mask)

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      call util_resize(arr, narr + nsrc)

      arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))

      return
   end subroutine util_append_arr_I4B


   module subroutine util_append_arr_logical(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical, dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: narr, nsrc

      if (.not. allocated(source)) return

      if (allocated(arr)) then
         narr = size(arr)
      else
         allocate(arr(nsrc))
         narr = 0
      end if

      nsrc = count(lsource_mask)

      call util_resize(arr, narr + nsrc)

      arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))

      return
   end subroutine util_append_arr_logical


   module subroutine util_append_body(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_body),            intent(inout) :: self         !! Swiftest body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      call util_append(self%name, source%name, lsource_mask)
      call util_append(self%id, source%id, lsource_mask)
      call util_append(self%status, source%status, lsource_mask)
      call util_append(self%ldiscard, source%ldiscard, lsource_mask)
      call util_append(self%lmask, source%lmask, lsource_mask)
      call util_append(self%mu, source%mu, lsource_mask)
      call util_append(self%xh, source%xh, lsource_mask)
      call util_append(self%vh, source%vh, lsource_mask)
      call util_append(self%xb, source%xb, lsource_mask)
      call util_append(self%vb, source%vb, lsource_mask)
      call util_append(self%ah, source%ah, lsource_mask)
      call util_append(self%aobl, source%aobl, lsource_mask)
      call util_append(self%atide, source%atide, lsource_mask)
      call util_append(self%agr, source%agr, lsource_mask)
      call util_append(self%ir3h, source%ir3h, lsource_mask)
      call util_append(self%a, source%a, lsource_mask)
      call util_append(self%e, source%e, lsource_mask)
      call util_append(self%inc, source%inc, lsource_mask)
      call util_append(self%capom, source%capom, lsource_mask)
      call util_append(self%omega, source%omega, lsource_mask)
      call util_append(self%capm, source%capm, lsource_mask)

      self%nbody = count(self%status(:) /= INACTIVE)

      return
   end subroutine util_append_body


   module subroutine util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_pl),              intent(inout) :: self         !! Swiftest massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to


      select type(source)
      class is (swiftest_pl)
         call util_append_body(self, source, lsource_mask)

         call util_append(self%mass, source%mass, lsource_mask)
         call util_append(self%Gmass, source%Gmass, lsource_mask)
         call util_append(self%rhill, source%rhill, lsource_mask)
         call util_append(self%radius, source%radius, lsource_mask)
         call util_append(self%xbeg, source%xbeg, lsource_mask)
         call util_append(self%xend, source%xend, lsource_mask)
         call util_append(self%vbeg, source%vbeg, lsource_mask)
         call util_append(self%density, source%density, lsource_mask)
         call util_append(self%Ip, source%Ip, lsource_mask)
         call util_append(self%rot, source%rot, lsource_mask)
         call util_append(self%k2, source%k2, lsource_mask)
         call util_append(self%Q, source%Q, lsource_mask)
         call util_append(self%tlag, source%tlag, lsource_mask)

         call self%eucl_index()
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_pl or its descendents"
         call util_exit(FAILURE)
      end select
   
      return
   end subroutine util_append_pl


   module subroutine util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_tp),              intent(inout) :: self         !! Swiftest test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (swiftest_tp)
         call util_append_body(self, source, lsource_mask)

         call util_append(self%isperi, source%isperi, lsource_mask)
         call util_append(self%peri, source%peri, lsource_mask)
         call util_append(self%atp, source%atp, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_tp or its descendents"
         call util_exit(FAILURE)
      end select

      return
   end subroutine util_append_tp

end submodule s_util_append