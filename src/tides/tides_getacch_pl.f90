submodule(swiftest_classes) s_tides_getacch
   use swiftest
contains
   module subroutine tides_getacch_pl(self, system)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! Calculated tidal torques from central body to any planet and from any planet to central body
      !! planet - planet interactions are considered negligeable
      !! Adapted from Mercury-T code from Bolmont et al. 2015
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      ! Internals
      integer(I4B)                          :: i
      real(DP)                              :: rmag
      real(DP), dimension(NDIM)             :: ej, vj, F_central
      real(DP), dimension(:,:), allocatable :: F_tot

      associate(pl => self, npl => self%nbody, cb => system%cb)
         allocate(F_tot, mold=pl%ah)
         do i = 1, npl
            ! Placeholders until model is implemented
            ! ***************************************
            F_tot(:,i) = 0.0_DP
            F_central(:) = 0.0_DP
            ! ***************************************
            !rmag = norm2(pl%xh(:,i))
            !ej = pl%xh(:,i) / rmag
            !vj = pl%vh(:, i)
            !Ftr = 
            !Pto = 
            !Pto_central =  !Eq 5 Bolmont et al. 2015 
            !F_tot(:,i) = (Ftr + (Pto + Pto_central) * dot_product(vj, ej) / rmag * ej + Pto * cross_product((rotj - theta_j), ej) + Pto_central * cross_product((rot_central - theta_j), ej) !Eq 6 Bolmont et al. 2015
            !F_central = F_central + F_tot(:,i)
         end do 

         do i = 1, npl
            pl%ah(:,i) = pl%ah(:,i) +  F_tot(:,i) / pl%Gmass(i) + F_central(:) / cb%Gmass
         end do 
      end associate

      return

   end subroutine tides_getacch_pl
end submodule s_tides_getacch