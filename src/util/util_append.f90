submodule (swiftest_classes) s_util_append
   use swiftest
contains

   module subroutine util_append_arr_char_string(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of character string type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      character(len=STRMAX), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      character(len=STRMAX), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                                     intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,               dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_char_string


   module subroutine util_append_arr_DP(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                        intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_DP


   module subroutine util_append_arr_DPvec(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of double precision vector type of size (NDIM, n) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      real(DP), dimension(:,:), allocatable, intent(inout) :: arr          !! Destination array 
      real(DP), dimension(:,:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                          intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,  dimension(:),                intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(NDIM,nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(1, nold + 1:nold + nnew) = pack(source(1,1:nsrc), lsource_mask(1:nsrc))
      arr(2, nold + 1:nold + nnew) = pack(source(2,1:nsrc), lsource_mask(1:nsrc))
      arr(3, nold + 1:nold + nnew) = pack(source(3,1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_DPvec


   module subroutine util_append_arr_I4B(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of integer(I4B) onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      integer(I4B), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      integer(I4B), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                            intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical,      dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

      return
   end subroutine util_append_arr_I4B


   module subroutine util_append_arr_logical(arr, source, nold, nsrc, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of logical type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      logical, dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      logical, dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      integer(I4B),                       intent(in)    :: nold, nsrc   !! Extend of the old array and the source array, respectively
      logical, dimension(:),              intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nnew

      if (.not. allocated(source)) return

      nnew = count(lsource_mask(1:nsrc))
      if (.not.allocated(arr)) then
         allocate(arr(nold+nnew))
      else
         call util_resize(arr, nold + nnew)
      end if

      arr(nold + 1:nold + nnew) = pack(source(1:nsrc), lsource_mask(1:nsrc))

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

      associate(nold => self%nbody, nsrc => source%nbody)
         call util_append(self%name, source%name, nold, nsrc, lsource_mask)
         call util_append(self%id, source%id, nold, nsrc, lsource_mask)
         call util_append(self%status, source%status, nold, nsrc, lsource_mask)
         call util_append(self%ldiscard, source%ldiscard, nold, nsrc, lsource_mask)
         call util_append(self%lmask, source%lmask, nold, nsrc, lsource_mask)
         call util_append(self%mu, source%mu, nold, nsrc, lsource_mask)
         call util_append(self%xh, source%xh, nold, nsrc, lsource_mask)
         call util_append(self%vh, source%vh, nold, nsrc, lsource_mask)
         call util_append(self%xb, source%xb, nold, nsrc, lsource_mask)
         call util_append(self%vb, source%vb, nold, nsrc, lsource_mask)
         call util_append(self%ah, source%ah, nold, nsrc, lsource_mask)
         call util_append(self%aobl, source%aobl, nold, nsrc, lsource_mask)
         call util_append(self%atide, source%atide, nold, nsrc, lsource_mask)
         call util_append(self%agr, source%agr, nold, nsrc, lsource_mask)
         call util_append(self%ir3h, source%ir3h, nold, nsrc, lsource_mask)
         call util_append(self%a, source%a, nold, nsrc, lsource_mask)
         call util_append(self%e, source%e, nold, nsrc, lsource_mask)
         call util_append(self%inc, source%inc, nold, nsrc, lsource_mask)
         call util_append(self%capom, source%capom, nold, nsrc, lsource_mask)
         call util_append(self%omega, source%omega, nold, nsrc, lsource_mask)
         call util_append(self%capm, source%capm, nold, nsrc, lsource_mask)
         self%nbody = nold + count(lsource_mask(:))
      end associate

      return
   end subroutine util_append_body


   module subroutine util_append_encounter(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(swiftest_encounter), intent(inout) :: self         !! Swiftest encounter list object
      class(swiftest_encounter), intent(in)    :: source       !! Source object to append
      logical, dimension(:),     intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      associate(nold => self%nenc, nsrc => source%nenc)
         call util_append(self%lvdotr, source%lvdotr, nold, nsrc, lsource_mask)
         call util_append(self%status, source%status, nold, nsrc, lsource_mask)
         call util_append(self%index1, source%index1, nold, nsrc, lsource_mask)
         call util_append(self%index2, source%index2, nold, nsrc, lsource_mask)
         call util_append(self%x1, source%x1, nold, nsrc, lsource_mask)
         call util_append(self%x2, source%x2, nold, nsrc, lsource_mask)
         call util_append(self%v1, source%v1, nold, nsrc, lsource_mask)
         call util_append(self%v2, source%v2, nold, nsrc, lsource_mask)
         call util_append(self%t, source%t, nold, nsrc, lsource_mask)
         self%nenc = nold + count(lsource_mask(:))
      end associate

      return
   end subroutine util_append_encounter


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
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%mass, source%mass, nold, nsrc, lsource_mask)
            call util_append(self%Gmass, source%Gmass, nold, nsrc, lsource_mask)
            call util_append(self%rhill, source%rhill, nold, nsrc, lsource_mask)
            call util_append(self%radius, source%radius, nold, nsrc, lsource_mask)
            call util_append(self%xbeg, source%xbeg, nold, nsrc, lsource_mask)
            call util_append(self%xend, source%xend, nold, nsrc, lsource_mask)
            call util_append(self%vbeg, source%vbeg, nold, nsrc, lsource_mask)
            call util_append(self%density, source%density, nold, nsrc, lsource_mask)
            call util_append(self%Ip, source%Ip, nold, nsrc, lsource_mask)
            call util_append(self%rot, source%rot, nold, nsrc, lsource_mask)
            call util_append(self%k2, source%k2, nold, nsrc, lsource_mask)
            call util_append(self%Q, source%Q, nold, nsrc, lsource_mask)
            call util_append(self%tlag, source%tlag, nold, nsrc, lsource_mask)

            call util_append_body(self, source, lsource_mask)
         end associate

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
         associate(nold => self%nbody, nsrc => source%nbody)
            call util_append(self%isperi, source%isperi, nold, nsrc, lsource_mask)
            call util_append(self%peri, source%peri, nold, nsrc, lsource_mask)
            call util_append(self%atp, source%atp, nold, nsrc, lsource_mask)

            call util_append_body(self, source, lsource_mask)
         end associate
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class swiftest_tp or its descendents"
         call util_exit(FAILURE)
      end select

      return
   end subroutine util_append_tp

end submodule s_util_append