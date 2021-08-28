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
         self%kidx(1:n) = source%kidx(1:n)
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


   module subroutine util_copy_particle_info(self, source)
      !! author: David A. Minton
      !!
      !! Copies one set of information object components into another, component-by-component
      implicit none
      class(swiftest_particle_info),  intent(inout) :: self
      class(swiftest_particle_info),  intent(in)    :: source

      self%name = source%name 
      self%particle_type = source%particle_type
      self%origin_type = source%origin_type
      self%origin_time = source%origin_time
      self%origin_xh(:) = source%origin_xh(:)
      self%origin_vh(:) = source%origin_vh(:)

      return
   end subroutine util_copy_particle_info


   module subroutine util_copy_particle_info_arr(source, dest, idx)
      !! author: David A. Minton
      !!
      !! Copies contents from an array of one particle information objects to another.
      implicit none
      class(swiftest_particle_info), dimension(:), intent(in)             :: source !! Source object to copy into
      class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata information object
      integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source object
      ! Internals
      integer(I4B) :: i, j, n, nsource, ndest

      if (size(source) == 0) return

      if (present(idx)) then
         n = size(idx)
      else
         n = size(source)
      end if

      nsource = size(source)
      ndest = size(dest)

      if ((n == 0) .or. (n > ndest) .or. (n > nsource)) then
         write(*,*) 'Particle info copy operation failed. n, nsource, ndest: ',n, nsource, ndest
         return
      end if

      do i = 1, n
         if (present(idx)) then
            j = idx(i)
         else
            j = i
         end if 
         call dest(i)%copy(source(j))
      end do
      
      return
   end subroutine util_copy_particle_info_arr



end submodule s_util_copy
