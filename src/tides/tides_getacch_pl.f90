submodule(swiftest_classes) s_tides_getacch
   use swiftest
contains
   module subroutine tides_getacch_pl(self, system)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! Calculated tidal torques from central body to any planet and from any planet to central body
      !! planet - planet interactions are considered negligable.
      !! This is a constant time lag model.
      !! 
      !! Adapted from Mercury-T code from Bolmont et al. (2015)
      !!
      !! Reference:
      !!    Bolmont, E., Raymond, S.N., Leconte, J., Hersant, F., Correia, A.C.M., 2015. 
      !!       Mercury-T : A new code to study tidally evolving multi-planet systems. 
      !!       Applications to Kepler-62. A&A 583, A116. https://doi.org/10.1051/0004-6361/201525909
      implicit none
      ! Arguments
      class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      ! Internals
      integer(I4B)                          :: i
      real(DP)                              :: rmag, vmag
      real(DP), dimension(NDIM)             :: r_unit, v_unit, h_unit, vj, F_central
      real(DP), dimension(:,:), allocatable :: F_tot

      associate(pl => self, npl => self%nbody, cb => system%cb)
         allocate(F_tot, mold=pl%ah)
         do i = 1, npl
            ! Placeholders until model is implemented
            ! ***************************************
            F_tot(:,i) = 0.0_DP
            F_central(:) = 0.0_DP
            ! ***************************************
            rmag = norm2(pl%xh(:,i))
            vmag = norm2(pl%vh(:,i))
            r_unit(:) = pl%xh(:,i) / rmag
            v_unit(:) = pl%vh(:,i) / vmag
            h_unit(:) = r_unit(:) .cross. v_unit(:)

             
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