submodule(symba_classes) s_symba_util
   use swiftest
contains
   module subroutine symba_util_copy_pltpenc(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(symba_pltpenc), intent(inout) :: self   !! SyMBA pl-tp encounter list 
      class(symba_pltpenc), intent(in)    :: source !! Source object to copy into

      associate(n => source%nenc)
         self%nenc = n
         self%lvdotr(1:n) = source%lvdotr(1:n) 
         self%status(1:n) = source%status(1:n) 
         self%level(1:n)  = source%level(1:n)
         self%index1(1:n) = source%index1(1:n)
         self%index2(1:n) = source%index2(1:n)
      end associate
   end subroutine symba_util_copy_pltpenc

   module subroutine symba_util_copy_plplenc(self, source)
      !! author: David A. Minton
      !!
      !! Copies elements from the source encounter list into self.
      implicit none
      ! Arguments
      class(symba_plplenc), intent(inout) :: self   !! SyMBA pl-pl encounter list 
      class(symba_pltpenc), intent(in)    :: source !! Source object to copy into

      call symba_util_copy_pltpenc(self, source)
      associate(n => source%nenc)
         select type(source)
         class is (symba_plplenc)
            self%xh1(:,1:n) = source%xh1(:,1:n) 
            self%xh2(:,1:n) = source%xh2(:,1:n) 
            self%vb1(:,1:n) = source%vb1(:,1:n) 
            self%vb2(:,1:n) = source%vb2(:,1:n) 
         end select
      end associate
   end subroutine symba_util_copy_plplenc

   module subroutine symba_util_resize_pltpenc(self, nrequested)
      !! author: David A. Minton
      !!
      !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
      !! Polymorphic method works on both symba_pltpenc and symba_plplenc types
      implicit none
      ! Arguments
      class(symba_pltpenc), intent(inout) :: self       !! SyMBA pl-tp encounter list 
      integer(I4B),         intent(in)    :: nrequested !! New size of list needed
      ! Internals
      class(symba_pltpenc), allocatable   :: enc_temp
      integer(I4B)                        :: nold

      nold = size(self%status)
      if (nrequested > nold) then
         allocate(enc_temp, source=self)
         call self%setup(2 * nrequested)
         call self%copy(enc_temp)
         deallocate(enc_temp)
      else
         self%status(nrequested+1:nold) = INACTIVE
      end if
      self%nenc = nrequested
      return
   end subroutine symba_util_resize_pltpenc

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
         pl%lcollision(1:npl) = pl_sorted%lcollision(ind(1:npl))
         pl%lencounter(1:npl) = pl_sorted%lencounter(ind(1:npl))
         pl%nplenc(1:npl) = pl_sorted%nplenc(ind(1:npl))
         pl%ntpenc(1:npl) = pl_sorted%ntpenc(ind(1:npl))
         pl%levelg(1:npl) = pl_sorted%levelg(ind(1:npl))
         pl%levelm(1:npl) = pl_sorted%levelm(ind(1:npl))
         pl%isperi(1:npl) = pl_sorted%isperi(ind(1:npl))
         pl%peri(1:npl) = pl_sorted%peri(ind(1:npl))
         pl%atp(1:npl) = pl_sorted%atp(ind(1:npl))
         pl%info(1:npl) = pl_sorted%info(ind(1:npl))
         pl%kin(1:npl) = pl_sorted%kin(ind(1:npl))
         do i = 1, npl
            do j = 1, pl%kin(i)%nchild
               pl%kin(i)%child(j) = ind(pl%kin(i)%child(j))
            end do
         end do
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
         tp%nplenc(1:ntp) = tp_sorted%nplenc(ind(1:ntp))
         tp%levelg(1:ntp) = tp_sorted%levelg(ind(1:ntp))
         tp%levelm(1:ntp) = tp_sorted%levelm(ind(1:ntp))
         deallocate(tp_sorted)
      end associate
      return
   end subroutine symba_util_sort_rearrange_tp

end submodule s_symba_util