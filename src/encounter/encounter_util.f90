submodule (encounter_classes) s_encounter_util
   use swiftest
contains

   module subroutine encounter_util_append_list(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one Swiftest body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self         !! Swiftest encounter list object
      class(encounter_list), intent(in)    :: source       !! Source object to append
      logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
      ! Internals
      integer(I4B) :: nold, nsrc

      nold = self%nenc
      nsrc = source%nenc
      call util_append(self%lvdotr, source%lvdotr, nold, nsrc, lsource_mask)
      call util_append(self%status, source%status, nold, nsrc, lsource_mask)
      call util_append(self%index1, source%index1, nold, nsrc, lsource_mask)
      call util_append(self%index2, source%index2, nold, nsrc, lsource_mask)
      call util_append(self%id1, source%id1, nold, nsrc, lsource_mask)
      call util_append(self%id2, source%id2, nold, nsrc, lsource_mask)
      call util_append(self%x1, source%x1, nold, nsrc, lsource_mask)
      call util_append(self%x2, source%x2, nold, nsrc, lsource_mask)
      call util_append(self%v1, source%v1, nold, nsrc, lsource_mask)
      call util_append(self%v2, source%v2, nold, nsrc, lsource_mask)
      call util_append(self%t, source%t, nold, nsrc, lsource_mask)
      self%nenc = nold + count(lsource_mask(1:nsrc))

      return
   end subroutine encounter_util_append_list


   module subroutine encounter_util_copy_list(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self   !! Encounter list 
      class(encounter_list), intent(in)    :: source !! Source object to copy into

      associate(n => source%nenc)
         self%nenc = n
         self%lvdotr(1:n) = source%lvdotr(1:n) 
         self%status(1:n) = source%status(1:n) 
         self%index1(1:n) = source%index1(1:n)
         self%index2(1:n) = source%index2(1:n)
         self%id1(1:n) = source%id1(1:n)
         self%id2(1:n) = source%id2(1:n)
         self%x1(:,1:n) = source%x1(:,1:n)
         self%x2(:,1:n) = source%x2(:,1:n)
         self%v1(:,1:n) = source%v1(:,1:n)
         self%v2(:,1:n) = source%v2(:,1:n)
         self%t(1:n) = source%t(1:n)
      end associate

      return
   end subroutine encounter_util_copy_list


   module subroutine encounter_util_collapse_ragged_list(ragged_list, n1, nenc, index1, index2, lvdotr)
      !! author: David A. Minton
      !!    
      !! Collapses a ragged index list (one encounter list per body) into a pair of index arrays and a vdotr logical array (optional)
      implicit none
      ! Arguments
      type(encounter_list), dimension(:),              intent(in)            :: ragged_list !! The ragged encounter list
      integer(I4B),                                    intent(in)            :: n1          !! Number of bodies 1
      integer(I4B),                                    intent(out)           :: nenc        !! Total number of encountersj 
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
      integer(I4B),         dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
      integer(I4B),         dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      ! Internals
      integer(I4B) :: i, j0, j1, nenci

      associate(nenc_arr => ragged_list(:)%nenc)
         nenc = sum(nenc_arr(:))
      end associate
      if (nenc == 0) return

      allocate(index1(nenc))
      allocate(index2(nenc))
      j0 = 1
      do i = 1, n1
         nenci = ragged_list(i)%nenc
         if (nenci > 0) then
            j1 = j0 + nenci - 1
            index1(j0:j1) = i
            index2(j0:j1) = ragged_list(i)%index2(:)
            if (present(lvdotr)) lvdotr(j0:j1) = ragged_list(i)%lvdotr(:)
            j0 = j1 + 1
         end if
      end do

      return
   end subroutine encounter_util_collapse_ragged_list


   module subroutine encounter_util_resize_list(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Note: The reason to extend it by a factor of 2 is for performance. When there are many enounters per step, resizing every time you want to add an 
      !! encounter takes significant computational effort. Resizing by a factor of 2 is a tradeoff between performance (fewer resize calls) and memory managment
      !! Memory usage grows by a factor of 2 each time it fills up, but no more. 
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self !! Swiftest encounter list 
      integer(I4B),              intent(in)    :: nnew !! New size of list needed
      ! Internals
      class(encounter_list), allocatable :: enc_temp
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
   end subroutine encounter_util_resize_list


   module subroutine encounter_util_sort_aabb_1D(self, n, extent_arr)
      !! author: David A. Minton
      !!
      !! Sorts the bounding box extents along a single dimension prior to the sweep phase. 
      !! This subroutine sets the sorted index array (ind) and the beginning/ending index list (beg & end)
      implicit none
      ! Arguments
      class(encounter_bounding_box_1D), intent(inout) :: self       !! Bounding box structure along a single dimension
      integer(I4B),                     intent(in)    :: n          !! Number of bodies with extents
      real(DP), dimension(:),           intent(in)    :: extent_arr !! Array of extents of size 2*n
      ! Internals
      integer(I4B) :: i, j, ibox, jbox
      logical, dimension(:), allocatable :: lfresh

      call util_sort(extent_arr, self%ind)
      allocate(lfresh(n))

      ! Determine the interval starting points and sizes
      lfresh(:) = .true. ! This will prevent double-counting of pairs
      do ibox = 1, 2*n
         i = self%ind(ibox)
         if (i > n) i = i - n ! If this is an endpoint index, shift it to the correct range
         if (.not.lfresh(i)) cycle
         do jbox = ibox + 1, 2*n
            j = self%ind(jbox)
            if (j > n) j = j - n ! If this is an endpoint index, shift it to the correct range
            if (j == i) then
               lfresh(i) = .false.
               self%ibeg(i) = ibox
               self%iend(i) = jbox
               exit ! We've reached the end of this interval 
            end if
         end do
      end do

      return
   end subroutine encounter_util_sort_aabb_1D


   module subroutine encounter_util_spill_list(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest encounter structure from active list to discard list
      implicit none
      ! Arguments
      class(encounter_list), intent(inout) :: self         !! Swiftest encounter list 
      class(encounter_list), intent(inout) :: discards     !! Discarded object 
      logical, dimension(:),     intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,                   intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      ! Internals
      integer(I4B) :: nenc_old
  
      associate(keeps => self)
         call util_spill(keeps%lvdotr, discards%lvdotr, lspill_list, ldestructive)
         call util_spill(keeps%status, discards%status, lspill_list, ldestructive)
         call util_spill(keeps%index1, discards%index1, lspill_list, ldestructive)
         call util_spill(keeps%index2, discards%index2, lspill_list, ldestructive)
         call util_spill(keeps%id1, discards%id1, lspill_list, ldestructive)
         call util_spill(keeps%id2, discards%id2, lspill_list, ldestructive)
         call util_spill(keeps%x1, discards%x1, lspill_list, ldestructive)
         call util_spill(keeps%x2, discards%x2, lspill_list, ldestructive)
         call util_spill(keeps%v1, discards%v1, lspill_list, ldestructive)
         call util_spill(keeps%v2, discards%v2, lspill_list, ldestructive)
         call util_spill(keeps%t, discards%t, lspill_list, ldestructive)

         nenc_old = keeps%nenc

         ! This is the base class, so will be the last to be called in the cascade. 
         ! Therefore we need to set the nenc values for both the keeps and discareds
         discards%nenc = count(lspill_list(1:nenc_old))
         if (ldestructive) keeps%nenc = nenc_old - discards%nenc
      end associate
   
      return
   end subroutine encounter_util_spill_list

   
   module subroutine encounter_util_sweep_aabb_double_list(self, n1, n2, ind_arr2, nenc, index1, index2)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the encounter candidates
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self     !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n1       !! Number of bodies 1
      integer(I4B),                            intent(in)    :: n2       !! Number of bodies 2
      integer(I4B), dimension(:),              intent(in)    :: ind_arr2 !! index array for mapping body 2 indexes
      integer(I4B),                            intent(out)   :: nenc     !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1   !! List of indices for body 1 in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2   !! List of indices for body 2 in each encounter candidate pair
      !Internals
      Integer(I4B) :: i, n
      type(encounter_list), dimension(n1) :: lenc         !! Array of encounter lists (one encounter list per body)

      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(self, lenc, ind_arr2) &
      !$omp firstprivate(n1, n2)
      do i = 1, n1
         call encounter_util_sweep_aabb_one_double_list(i, n1, n2, self%aabb(1)%ind(:), self%aabb(1)%ibeg(:), self%aabb(1)%iend(:), self%aabb(2)%ibeg(:), self%aabb(2)%iend(:), ind_arr2, lenc(i))
      end do
      !$omp end parallel do 

      call encounter_util_collapse_ragged_list(lenc, n1, nenc, index1, index2)

      return
   end subroutine encounter_util_sweep_aabb_double_list


   module subroutine encounter_util_sweep_aabb_single_list(self, n, ind_arr, nenc, index1, index2)
      !! author: David A. Minton
      !!
      !! Sweeps the sorted bounding box extents and returns the encounter candidates. Mutual encounters
      !! allowed. That is, all bodies are from the same list
      implicit none
      ! Arguments
      class(encounter_bounding_box),           intent(inout) :: self    !! Multi-dimensional bounding box structure
      integer(I4B),                            intent(in)    :: n       !! Number of bodies 1
      integer(I4B), dimension(:),              intent(in)    :: ind_arr !! index array for mapping body 2 indexes
      integer(I4B),                            intent(out)   :: nenc    !! Total number of encounter candidates
      integer(I4B), dimension(:), allocatable, intent(out)   :: index1  !! List of indices for body 1 in each encounter candidate pair
      integer(I4B), dimension(:), allocatable, intent(out)   :: index2  !! List of indices for body 2 in each encounter candidate pair
      !Internals
      Integer(I4B) :: i
      type(encounter_list), dimension(n) :: lenc         !! Array of encounter lists (one encounter list per body)

      ! Sweep the intervals for each of the massive bodies along one dimension
      ! This will build a ragged pair of index lists inside of the lenc data structure
      !$omp parallel do default(private) schedule(static)&
      !$omp shared(self, lenc, ind_arr) &
      !$omp firstprivate(n)
      do i = 1, n
         call encounter_util_sweep_aabb_one_single_list(i, n, self%aabb(1)%ind(:), self%aabb(1)%ibeg(:), self%aabb(1)%iend(:), self%aabb(2)%ibeg(:), self%aabb(2)%iend(:), ind_arr, lenc(i))
      end do
      !$omp end parallel do 

      call encounter_util_collapse_ragged_list(lenc, n, nenc, index1, index2)

      return
   end subroutine encounter_util_sweep_aabb_single_list


   subroutine encounter_util_sweep_aabb_one_double_list(i, n1, n2, ext_ind, ibegx, iendx, ibegy, iendy, ind_arr2, lenc)
      !! author: David A. Minton
      !!
      !! Performs a sweep operation on a single body. Encounters from the same lists not allowed (e.g. pl-tp encounters only)
      implicit none
      ! Arguments
      integer(I4B),               intent(in)    :: i            !! The current index of the ith body
      integer(I4B),               intent(in)    :: n1           !! Number of bodies 1
      integer(I4B),               intent(in)    :: n2           !! Number of bodies 2
      integer(I4B), dimension(:), intent(in)    :: ext_ind      !! Sorted index array of extents
      integer(I4B), dimension(:), intent(in)    :: ibegx, iendx !! Beginning and ending index lists in the x-dimension
      integer(I4B), dimension(:), intent(in)    :: ibegy, iendy !! Beginning and ending index lists in the y-dimension
      integer(I4B), dimension(:), intent(in)    :: ind_arr2     !! index array for mapping body 2 indexes
      type(encounter_list),       intent(inout) :: lenc         !! Encounter list for the ith body
      ! Internals
      integer(I4B) :: ibox, jbox, nbox, j, ybegi, yendi, ntot
      logical, dimension(n1) :: lencounteri

      ntot = n1 + n2
      ibox = ibegx(i) + 1
      nbox = iendx(i) - 1
      ybegi = ibegy(i) 
      yendi = iendy(i)
      lencounteri(:) = .false.
      do concurrent(jbox = ibox:nbox) ! Sweep forward until the end of the interval
         j = ext_ind(jbox)
         if (j > ntot) j = j - ntot ! If this is an endpoint index, shift it to the correct range
         if (((i <= n1) .and. (j <= n1)) .or. ((i > n1) .and. (j > n1))) cycle  ! only pairs
         ! Check the y-dimension
         lencounteri(j) = (iendy(j) > ybegi) .and. (ibegy(j) < yendi)
      end do

      lenc%nenc = count(lencounteri(:))
      if (lenc%nenc > 0) then
         allocate(lenc%index2(lenc%nenc))
         lenc%index2(:) = pack(ind_arr2(:), lencounteri(:)) 
      end if

      return
   end subroutine encounter_util_sweep_aabb_one_double_list
 

   subroutine encounter_util_sweep_aabb_one_single_list(i, n, ext_ind, ibegx, iendx, ibegy, iendy, ind_arr, lenc)
      !! author: David A. Minton
      !!
      !! Performs a sweep operation on a single body. Mutual encounters allowed (e.g. pl-pl)
      implicit none
      ! Arguments
      integer(I4B),               intent(in)    :: i            !! The current index of the ith body
      integer(I4B),               intent(in)    :: n            !! Number of bodies
      integer(I4B), dimension(:), intent(in)    :: ext_ind      !! Sorted index array of extents
      integer(I4B), dimension(:), intent(in)    :: ibegx, iendx !! Beginning and ending index lists in the x-dimension
      integer(I4B), dimension(:), intent(in)    :: ibegy, iendy !! Beginning and ending index lists in the y-dimension
      integer(I4B), dimension(:), intent(in)    :: ind_arr      !! index array for mapping body 2 indexes
      type(encounter_list),       intent(inout) :: lenc         !! Encounter list for the ith body
      ! Internals
      integer(I4B) :: ibox, jbox, nbox, j, ybegi, yendi
      logical, dimension(n) :: lencounteri

      ibox = ibegx(i) + 1
      nbox = iendx(i) - 1
      ybegi = ibegy(i) 
      yendi = iendy(i)
      lencounteri(:) = .false.
      do concurrent(jbox = ibox:nbox) ! Sweep forward until the end of the interval
         j = ext_ind(jbox)
         if (j > n) j = j - n ! If this is an endpoint index, shift it to the correct range
         ! Check the y-dimension
         lencounteri(j) = (iendy(j) > ybegi) .and. (ibegy(j) < yendi)
      end do

      lenc%nenc = count(lencounteri(:))
      if (lenc%nenc > 0) then
         allocate(lenc%index2(lenc%nenc))
         lenc%index2(:) = pack(ind_arr(:), lencounteri(:)) 
      end if

      return
   end subroutine encounter_util_sweep_aabb_one_single_list


end submodule s_encounter_util