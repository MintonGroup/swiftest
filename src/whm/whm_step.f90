submodule(whm_classes) s_whm_step
contains
   module procedure whm_step
      !! author: David A. Minton
      !!
      !! Step planets and active test particles ahead in heliocentric coordinates
      !!
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      use swiftest
      implicit none
      logical                       :: lfirsttp
      logical, save                 :: lmalloc = .true.
      integer(I4B)                     :: i
      real(DP), dimension(:, :), allocatable, save :: xbeg, xend
      type(whm_pl), pointer                :: whm_plp
   
      lfirsttp = lfirst
      associate(ntp => self%tp%nbody, npl => self%pl%nbody)
         if (ntp > 0) xbeg(:, :) = self%pl%xh(:, 1:npl)
         call self%pl%step(config, t, dt)
         if (ntp > 0) then 
            xend(:, :) = self%pl%xh(:, 1:npl)
            call self%tp%step(config, t, dt, xbeg, xend)
         end if
      end associate
      self%lintegrate = .true.
      return
   end procedure whm_step
end submodule s_whm_step
