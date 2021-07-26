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
      end if
      self%nenc = nrequested
      return
   end subroutine symba_util_resize_pltpenc


end submodule s_symba_util