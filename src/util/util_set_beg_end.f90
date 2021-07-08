submodule(swiftest_classes) s_util_set_beg_end
   use swiftest
contains
   
   module subroutine util_set_beg_end(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
      implicit none
      ! Arguments
      class(swiftest_pl),       intent(inout)          :: self !! Swiftest massive body object
      real(DP), dimension(:,:), intent(in),   optional :: xbeg, xend, vbeg

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end subroutine util_set_beg_end
end submodule s_util_set_beg_end