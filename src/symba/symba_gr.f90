submodule(symba_classes) s_symba_gr
   use swiftest
contains

   module pure subroutine symba_gr_p4_pl(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to massive bodies due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_symba_p4.f90
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B)                                 :: i

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         select type(system)
         class is (symba_nbody_system)
            do concurrent(i = 1:npl, pl%lmask(i) .and. pl%levelg(i) == system%irec )
               call gr_p4_pos_kick(param, pl%xh(:, i), pl%vb(:, i), dt)
            end do
         end select
      end associate

   return
   end subroutine symba_gr_p4_pl


   module pure subroutine symba_gr_p4_tp(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Position kick to test particles due to p**4 term in the post-Newtonian correction
      !!    Based on Saha & Tremaine (1994) Eq. 28
      !!
      !! Adapted from David A. Minton's Swifter routine routine gr_symba_p4.f90
      implicit none
      ! Arguments
      class(symba_tp),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Step size
      ! Internals
      integer(I4B)                              :: i

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody)
         select type(system)
         class is (symba_nbody_system)
            do concurrent(i = 1:ntp, tp%lmask(i) .and. tp%levelg(i) == system%irec)
               call gr_p4_pos_kick(param, tp%xh(:, i), tp%vb(:, i), dt)
            end do
         end select
      end associate

   return
   end subroutine symba_gr_p4_tp

end submodule s_symba_gr
