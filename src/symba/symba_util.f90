submodule(symba_classes) s_symba_util
   use swiftest
contains

   module subroutine symba_util_append_arr_info(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of particle information type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(symba_particle_info), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      type(symba_particle_info), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical,                   dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
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

      call util_resize(arr, narr + nsrc)

      if (present(lsource_mask)) then
         arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))
      else
         arr(narr + 1:narr + nsrc) = source(:)
      end if

      return
   end subroutine symba_util_append_arr_info


   module subroutine symba_util_append_arr_kin(arr, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append a single array of kinship type onto another. If the destination array is not allocated, or is not big enough, this will allocate space for it.
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr          !! Destination array 
      type(symba_kinship), dimension(:), allocatable, intent(in)    :: source       !! Array to append 
      logical,             dimension(:), optional,    intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
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

      call util_resize(arr, narr + nsrc)

      if (present(lsource_mask)) then
         arr(narr + 1:narr + nsrc) = pack(source(:), lsource_mask(:))
      else
         arr(narr + 1:narr + nsrc) = source(:)
      end if

      return
   end subroutine symba_util_append_arr_kin


   module subroutine symba_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_pl),                 intent(inout) :: self         !! SyMBA massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_pl)
         call util_append_pl(self, source, lsource_mask) ! Note: helio_pl does not have its own append method, so we skip back to the base class

         call util_append(self%lcollision, source%lcollision, lsource_mask)
         call util_append(self%lencounter, source%lencounter, lsource_mask)
         call util_append(self%lmtiny, source%lmtiny, lsource_mask)
         call util_append(self%nplenc, source%nplenc, lsource_mask)
         call util_append(self%ntpenc, source%ntpenc, lsource_mask)
         call util_append(self%levelg, source%levelg, lsource_mask)
         call util_append(self%levelm, source%levelm, lsource_mask)
         call util_append(self%isperi, source%isperi, lsource_mask)
         call util_append(self%peri, source%peri, lsource_mask)
         call util_append(self%atp, source%atp, lsource_mask)
         call util_append(self%kin, source%kin, lsource_mask)
         call util_append(self%info, source%info, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_pl or its descendents!"
         call util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_pl


   module subroutine symba_util_append_merger(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(symba_merger),             intent(inout) :: self         !! SyMBA massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B), dimension(:), allocatable        :: ncomp_tmp    !! Temporary placeholder for ncomp incase we are appending a symba_pl object to a symba_merger

      select type(source)
      class is (symba_merger)
         call symba_util_append_pl(self, source, lsource_mask) 
         call util_append(self%ncomp, source%ncomp, lsource_mask)
      class is (symba_pl)
         call symba_util_append_pl(self, source, lsource_mask) 
         allocate(ncomp_tmp, mold=source%id)
         ncomp_tmp(:) = 0
         call util_append(self%ncomp, ncomp_tmp, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_pl or its descendents!"
         call util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_merger


   module subroutine symba_util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from test particle object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_tp),                 intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:), optional, intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_tp)
         call util_append_tp(self, source, lsource_mask) ! Note: helio_tp does not have its own append method, so we skip back to the base class

         call util_append(self%nplenc, source%nplenc, lsource_mask)
         call util_append(self%levelg, source%levelg, lsource_mask)
         call util_append(self%levelm, source%levelm, lsource_mask)
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_tp or its descendents!"
         call util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_tp


   module subroutine symba_util_fill_arr_info(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(symba_particle_info), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(symba_particle_info), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,                   dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
      return
   end subroutine symba_util_fill_arr_info


   module subroutine symba_util_fill_arr_kin(keeps, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Performs a fill operation on a single array of particle kinship types
      !! This is the inverse of a spill operation   
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps      !! Array of values to keep 
      type(symba_kinship), dimension(:), allocatable, intent(in)    :: inserts    !! Array of values to insert into keep
      logical,             dimension(:),              intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      if (.not.allocated(keeps) .or. .not.allocated(inserts)) return

      keeps(:) = unpack(keeps(:),   .not.lfill_list(:), keeps(:))
      keeps(:) = unpack(inserts(:),      lfill_list(:), keeps(:))
   
      return
   end subroutine symba_util_fill_arr_kin


   module subroutine symba_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self       !! SyMBA masive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_pl)
            call util_fill(keeps%lcollision, inserts%lcollision, lfill_list)
            call util_fill(keeps%lencounter, inserts%lencounter, lfill_list)
            call util_fill(keeps%lmtiny, inserts%lmtiny, lfill_list)
            call util_fill(keeps%nplenc, inserts%nplenc, lfill_list)
            call util_fill(keeps%ntpenc, inserts%ntpenc, lfill_list)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)
            call util_fill(keeps%isperi, inserts%isperi, lfill_list)
            call util_fill(keeps%peri, inserts%peri, lfill_list)
            call util_fill(keeps%atp, inserts%atp, lfill_list)
            call util_fill(keeps%kin, inserts%kin, lfill_list)
            call util_fill(keeps%info, inserts%info, lfill_list)
            
            call util_fill_pl(keeps, inserts, lfill_list)  ! Note: helio_pl does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_pl


   module subroutine symba_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self       !! SyMBA test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_tp)
            call util_fill(keeps%nplenc, inserts%nplenc, lfill_list)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)
            
            call util_fill_tp(keeps, inserts, lfill_list) ! Note: helio_tp does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_tp or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_tp


   module subroutine symba_util_peri_pl(self, system, param)
      !! author: David A. Minton
      !!
      !! Determine system pericenter passages for planets in SyMBA
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_peri.f90
      !! Adapted from Hal Levison's Swift routine util_mass_peri.f
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)       :: i
      real(DP)           :: vdotr, e

      associate(pl => self, npl => self%nbody)
         if (pl%lfirst) then
            if (param%qmin_coord == "HELIO") then
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xh(:,i), pl%vh(:,i))
                     if (vdotr > 0.0_DP) then
                        pl%isperi(i) = 1
                     else
                        pl%isperi(i) = -1
                     end if
                  end if
               end do
            else
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xb(:,i), pl%vb(:,i))
                     if (vdotr > 0.0_DP) then
                        pl%isperi(i) = 1
                     else
                        pl%isperi(i) = -1
                     end if
                  end if
               end do
            end if
         else
            if (param%qmin_coord == "HELIO") then
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xh(:,i), pl%vh(:,i))
                     if (pl%isperi(i) == -1) then
                        if (vdotr >= 0.0_DP) then
                           pl%isperi(i) = 0
                           CALL orbel_xv2aeq(pl%mu(i), pl%xh(:,i), pl%vh(:,i), pl%atp(i), e, pl%peri(i))
                        end if
                     else
                        if (vdotr > 0.0_DP) then
                           pl%isperi(i) = 1
                        else
                           pl%isperi(i) = -1
                        end if
                     end if
                  end if
               end do
            else
               do i = 1, npl
                  if (pl%status(i) == ACTIVE) then
                     vdotr = dot_product(pl%xb(:,i), pl%vb(:,i))
                     if (pl%isperi(i) == -1) then
                        if (vdotr >= 0.0_DP) then
                           pl%isperi(i) = 0
                           CALL orbel_xv2aeq(system%Gmtot, pl%xb(:,i), pl%vb(:,i), pl%atp(i), e, pl%peri(i))
                        end if
                     else
                        if (vdotr > 0.0_DP) then
                           pl%isperi(i) = 1
                        else
                           pl%isperi(i) = -1
                        end if
                     end if
                  end if
               end do
            end if
         end if
      end associate
 
     return
   end subroutine symba_util_peri_pl


   module subroutine symba_util_rearray_pl(self, system, param)
      !! Author: the Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Clean up the massive body structures to remove discarded bodies and add new bodies
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: self   !! SyMBA massive body object
      class(symba_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(symba_parameters),   intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      class(symba_pl), allocatable :: tmp !! The discarded body list.

      associate(pl => self, pl_adds => system%pl_adds)
         allocate(tmp, mold=pl)
         ! Remove the discards and destroy the list, as the system already tracks pl_discards elsewhere
         call pl%spill(tmp, lspill_list=(pl%ldiscard(:) .or. pl%status(:) == INACTIVE), ldestructive=.true.)
         call tmp%setup(0,param)
         deallocate(tmp)

         ! Add in any new bodies
         call pl%append(pl_adds)

         ! If there are still bodies in the system, sort by mass in descending order and re-index
         if (pl%nbody > 0) then
            call pl%sort("mass", ascending=.false.)
            pl%lmtiny(:) = pl%Gmass(:) > param%MTINY
            pl%nplm = count(pl%lmtiny(:))
            call pl%eucl_index()
         end if

      end associate

      return
   end subroutine symba_util_rearray_pl


   module subroutine symba_util_resize_arr_info(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(symba_particle_info), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                         intent(in)    :: nnew !! New size
      ! Internals
      type(symba_particle_info), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
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
   end subroutine symba_util_resize_arr_info


   module subroutine symba_util_resize_arr_kin(arr, nnew)
      !! author: David A. Minton
      !!
      !! Resizes an array component of type character string. Array will only be resized if has previously been allocated. Passing nnew = 0 will deallocate.
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: arr  !! Array to resize
      integer(I4B),                                   intent(in)    :: nnew !! New size
      ! Internals
      type(symba_kinship), dimension(:), allocatable :: tmp !! Temporary storage array in case the input array is already allocated
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
   end subroutine symba_util_resize_arr_kin


   module subroutine symba_util_resize_merger(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a SyMBA merger list against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_merger), intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),        intent(in)    :: nnew  !! New size neded

      call symba_util_resize_pl(self, nnew)

      call util_resize(self%ncomp, nnew)

      return
   end subroutine symba_util_resize_merger


   module subroutine symba_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a SyMBA massive body object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize_pl(self, nnew)

      call util_resize(self%lcollision, nnew)
      call util_resize(self%lencounter, nnew)
      call util_resize(self%lmtiny, nnew)
      call util_resize(self%nplenc, nnew)
      call util_resize(self%ntpenc, nnew)
      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)
      call util_resize(self%isperi, nnew)
      call util_resize(self%peri, nnew)
      call util_resize(self%atp, nnew)
      call util_resize(self%kin, nnew)
      call util_resize(self%info, nnew)

      return
   end subroutine symba_util_resize_pl


   module subroutine symba_util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a test particle object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize_tp(self, nnew)

      call util_resize(self%nplenc, nnew)
      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)

      return
   end subroutine symba_util_resize_tp


   module subroutine symba_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self      !! SyMBA massive body object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * pl%nplenc(1:npl), ind(1:npl))
         case("ntpenc")
            call util_sort(direction * pl%ntpenc(1:npl), ind(1:npl))
         case("levelg")
            call util_sort(direction * pl%levelg(1:npl), ind(1:npl))
         case("levelm")
            call util_sort(direction * pl%levelm(1:npl), ind(1:npl))
         case("peri")
            call util_sort(direction * pl%peri(1:npl), ind(1:npl))
         case("atp")
            call util_sort(direction * pl%atp(1:npl), ind(1:npl))
         case("lcollision", "lencounter", "lmtiny", "nplm", "nplplm", "kin", "info")
            write(*,*) 'Cannot sort by ' // trim(adjustl(sortby)) // '. Component not sortable!'
         case default ! Look for components in the parent class
            call util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate
      return
   end subroutine symba_util_sort_pl


   module subroutine symba_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA test particle object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self      !! SyMBA test particle object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(self%nbody) :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * tp%nplenc(1:ntp), ind(1:ntp))
         case("levelg")
            call util_sort(direction * tp%levelg(1:ntp), ind(1:ntp))
         case("levelm")
            call util_sort(direction * tp%levelm(1:ntp), ind(1:ntp))
         case default ! Look for components in the parent class
            call util_sort_tp(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)
      end associate

      return
   end subroutine symba_util_sort_tp


   module subroutine symba_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_pl),               intent(inout) :: self !! SyMBA massive body object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      ! Internals
      class(symba_pl), allocatable :: pl_sorted  !! Temporary holder for sorted body
      integer(I4B) :: i, j

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange_pl(pl,ind)
         allocate(pl_sorted, source=self)
         if (allocated(pl%lcollision)) pl%lcollision(1:npl) = pl_sorted%lcollision(ind(1:npl))
         if (allocated(pl%lencounter)) pl%lencounter(1:npl) = pl_sorted%lencounter(ind(1:npl))
         if (allocated(pl%lmtiny))     pl%lmtiny(1:npl) = pl_sorted%lmtiny(ind(1:npl))
         if (allocated(pl%nplenc))     pl%nplenc(1:npl) = pl_sorted%nplenc(ind(1:npl))
         if (allocated(pl%ntpenc))     pl%ntpenc(1:npl) = pl_sorted%ntpenc(ind(1:npl))
         if (allocated(pl%levelg))     pl%levelg(1:npl) = pl_sorted%levelg(ind(1:npl))
         if (allocated(pl%levelm))     pl%levelm(1:npl) = pl_sorted%levelm(ind(1:npl))
         if (allocated(pl%isperi))     pl%isperi(1:npl) = pl_sorted%isperi(ind(1:npl))
         if (allocated(pl%peri))       pl%peri(1:npl) = pl_sorted%peri(ind(1:npl))
         if (allocated(pl%atp))        pl%atp(1:npl) = pl_sorted%atp(ind(1:npl))
         if (allocated(pl%info))       pl%info(1:npl) = pl_sorted%info(ind(1:npl))
         if (allocated(pl%kin)) then
            pl%kin(1:npl) = pl_sorted%kin(ind(1:npl))
            do i = 1, npl
               do j = 1, pl%kin(i)%nchild
                  pl%kin(i)%child(j) = ind(pl%kin(i)%child(j))
               end do
            end do
         end if
         deallocate(pl_sorted)
      end associate

      return
   end subroutine symba_util_sort_rearrange_pl


   module subroutine symba_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA test particle object in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_tp),               intent(inout) :: self !! SyMBA test particle object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)
      ! Internals
      class(symba_tp), allocatable :: tp_sorted  !! Temporary holder for sorted body

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange_tp(tp,ind)
         allocate(tp_sorted, source=self)
         if (allocated(tp%nplenc)) tp%nplenc(1:ntp) = tp_sorted%nplenc(ind(1:ntp))
         if (allocated(tp%levelg)) tp%levelg(1:ntp) = tp_sorted%levelg(ind(1:ntp))
         if (allocated(tp%levelm)) tp%levelm(1:ntp) = tp_sorted%levelm(ind(1:ntp))
         deallocate(tp_sorted)
      end associate
      
      return
   end subroutine symba_util_sort_rearrange_tp


   module subroutine symba_util_spill_arr_info(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle origin information types
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(symba_particle_info), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(symba_particle_info), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,                   dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                                              intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not

      if (.not.allocated(keeps) .or. count(lspill_list(:)) == 0) return
      if (.not.allocated(discards)) allocate(discards(count(lspill_list(:))))

      discards(:) = pack(keeps(:), lspill_list(:))
      if (ldestructive) then
         if (count(.not.lspill_list(:)) > 0) then
            keeps(:) = pack(keeps(:), .not. lspill_list(:))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine symba_util_spill_arr_info


   module subroutine symba_util_spill_arr_kin(keeps, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Performs a spill operation on a single array of particle kinships
      !! This is the inverse of a spill operation
      implicit none
      ! Arguments
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: keeps        !! Array of values to keep 
      type(symba_kinship), dimension(:), allocatable, intent(inout) :: discards     !! Array of discards
      logical,             dimension(:),              intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discardss
      logical,                                        intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter the keeps array or not

      if (.not.allocated(keeps) .or. count(lspill_list(:)) == 0) return
      if (.not.allocated(discards)) allocate(discards(count(lspill_list(:))))

      discards(:) = pack(keeps(:), lspill_list(:))
      if (ldestructive) then
         if (count(.not.lspill_list(:)) > 0) then
            keeps(:) = pack(keeps(:), .not. lspill_list(:))
         else
            deallocate(keeps)
         end if
      end if

      return
   end subroutine symba_util_spill_arr_kin


   module subroutine symba_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA massive body particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self        !! SyMBA massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_pl)
            call util_spill(keeps%lcollision, discards%lcollision, lspill_list, ldestructive)
            call util_spill(keeps%lencounter, discards%lencounter, lspill_list, ldestructive)
            call util_spill(keeps%lmtiny, discards%lmtiny, lspill_list, ldestructive)
            call util_spill(keeps%nplenc, discards%nplenc, lspill_list, ldestructive)
            call util_spill(keeps%ntpenc, discards%ntpenc, lspill_list, ldestructive)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)
            call util_spill(keeps%isperi, discards%isperi, lspill_list, ldestructive)
            call util_spill(keeps%peri, discards%peri, lspill_list, ldestructive)
            call util_spill(keeps%atp, discards%atp, lspill_list, ldestructive)
            call util_spill(keeps%info, discards%info, lspill_list, ldestructive)
            call util_spill(keeps%kin, discards%kin, lspill_list, ldestructive)

            call util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_pl or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_pl


   module subroutine symba_util_spill_pltpenc(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA encounter structure from active list to discard list
      !! Note: Because the symba_plplenc currently does not contain any additional variable components, this method can recieve it as an input as well.
      implicit none
      ! Arguments
      class(symba_pltpenc),      intent(inout) :: self         !! SyMBA pl-tp encounter list 
      class(swiftest_encounter), intent(inout) :: discards     !! Discarded object 
      logical, dimension(:),     intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                   intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I4B) :: i
  
      associate(keeps => self)
         select type(discards)
         class is (symba_pltpenc)
            call util_spill(keeps%level, discards%level, lspill_list, ldestructive)
            call util_spill_encounter(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_pltpenc or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
   
      return
   end subroutine symba_util_spill_pltpenc


   module subroutine symba_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I4B) :: i

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_tp)
            call util_spill(keeps%nplenc, discards%nplenc, lspill_list, ldestructive)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)

            call util_spill_tp(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_tp or its descendents!"
            call util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_tp

end submodule s_symba_util