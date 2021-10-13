submodule (swiftest_classes) s_util_dealloc
   use swiftest
contains

   module subroutine util_dealloc_body(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_body),  intent(inout) :: self

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%lmask)) deallocate(self%lmask)
      if (allocated(self%mu)) deallocate(self%mu)
      if (allocated(self%xh)) deallocate(self%xh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%xb)) deallocate(self%xb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%ah)) deallocate(self%ah)
      if (allocated(self%aobl)) deallocate(self%aobl)
      if (allocated(self%agr)) deallocate(self%lmask)
      if (allocated(self%atide)) deallocate(self%lmask)
      if (allocated(self%ir3h)) deallocate(self%ir3h)
      if (allocated(self%a)) deallocate(self%a)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%inc)) deallocate(self%inc)
      if (allocated(self%capom)) deallocate(self%capom)
      if (allocated(self%omega)) deallocate(self%omega)
      if (allocated(self%capm)) deallocate(self%capm)

      return
   end subroutine util_dealloc_body


   module subroutine util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest massive body object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_pl),  intent(inout) :: self !! Swiftest massive body object

      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%Gmass)) deallocate(self%Gmass)
      if (allocated(self%rhill)) deallocate(self%rhill)
      if (allocated(self%renc)) deallocate(self%renc)
      if (allocated(self%radius)) deallocate(self%radius)
      if (allocated(self%density)) deallocate(self%density)
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%Ip)) deallocate(self%Ip)
      if (allocated(self%k2)) deallocate(self%k2)
      if (allocated(self%Q)) deallocate(self%Q)
      if (allocated(self%tlag)) deallocate(self%tlag)
      if (allocated(self%k_plpl)) deallocate(self%k_plpl)

      call util_dealloc_body(self)

      return
   end subroutine util_dealloc_pl


   module subroutine util_dealloc_system(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest nbody system object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_nbody_system),  intent(inout) :: self !! Swiftest nbody system object

      call self%pl%dealloc()
      call self%tp%dealloc()
      call self%tp_discards%dealloc()
      call self%pl_discards%dealloc()

      if (allocated(self%cb)) deallocate(self%cb)
      if (allocated(self%pl)) deallocate(self%pl)
      if (allocated(self%tp)) deallocate(self%tp)
      if (allocated(self%tp_discards)) deallocate(self%tp_discards)
      if (allocated(self%pl_discards)) deallocate(self%pl_discards)

      return
   end subroutine util_dealloc_system


   module subroutine util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Finalize the swiftest test particle object - deallocates all allocatables
      implicit none
      ! Argument
      class(swiftest_tp),  intent(inout) :: self !! Swiftest test particle object

      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)
      if (allocated(self%k_pltp)) deallocate(self%k_pltp)

      call util_dealloc_body(self)

      return
   end subroutine util_dealloc_tp

end submodule s_util_dealloc