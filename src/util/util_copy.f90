submodule(swiftest_classes) s_util_copy
   use swiftest
contains

module subroutine util_copy_encounter(self, source)
   !! author: David A. Minton
   !!
   !! Copies elements from the source encounter list into self.
   implicit none
   ! Arguments
   class(swiftest_encounter), intent(inout) :: self   !! Encounter list 
   class(swiftest_encounter), intent(in)    :: source !! Source object to copy into

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
end subroutine util_copy_encounter

end submodule s_util_copy
