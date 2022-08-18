submodule(helio_classes) s_helio_util
   use swiftest
contains

   module subroutine helio_util_final_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio massive body object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_pl),  intent(inout) :: self !! Helio massive body object

      call self%dealloc()

      return
   end subroutine helio_util_final_pl


   module subroutine helio_util_final_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio nbody system object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_nbody_system),  intent(inout) :: self !! Helio nbody system object

      call self%dealloc()

      return
   end subroutine helio_util_final_system


   module subroutine helio_util_final_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the Helio test particle object - deallocates all allocatables
      implicit none
      ! Arguments
      type(helio_tp),  intent(inout) :: self !! Helio test particle object

      call self%dealloc()

      return
   end subroutine helio_util_final_tp

end submodule s_helio_util