submodule(tides) s_tides_kick_getacch
   use swiftest
contains

   module subroutine tides_kick_getacch_pl(self, nbody_system)
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
      class(base_object),           intent(inout) :: self   !! Swiftest massive body object
      class(base_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      ! Internals
      integer(I4B)                          :: i
      real(DP)                              :: rmag, vmag
      real(DP), dimension(NDIM)             :: r_unit, v_unit, h_unit, theta_unit, theta_dot, F_T
      real(DP)                              :: Ftr, Ptopl, Ptocb, r5cbterm, r5plterm

      select type(pl => self)
      class is (swiftest_pl)
      select type(nbody_system)
      class is (swiftest_nbody_system)
         associate(npl => pl%nbody, cb => nbody_system%cb)
            pl%atide(:,:) = 0.0_DP
            cb%atide(:) = 0.0_DP
            do i = 1, npl
               rmag = norm2(pl%rh(:,i))
               vmag = norm2(pl%vh(:,i))
               r_unit(:) = pl%rh(:,i) / rmag
               v_unit(:) = pl%vh(:,i) / vmag
               h_unit(:) = r_unit(:) .cross. v_unit(:)
               theta_unit(:) = h_unit(:) .cross. r_unit(:)
               theta_dot = dot_product(pl%vh(:,i), theta_unit(:))

               ! First calculate the tangential component of the force vector (eq. 5 & 6 of Bolmont et al. 2015)
               ! The radial component is already computed in the obl_acc methods
               r5cbterm = pl%Gmass(i)**2 * cb%k2 * cb%radius**5
               r5plterm = cb%Gmass**2 * pl%k2(i) * pl%radius(i)**5

               Ptopl = 3 * r5plterm * pl%tlag(i) / rmag**7
               Ptocb = 3 * r5cbterm * cb%tlag / rmag**7

               Ftr = -3 / rmag**7 * (r5cbterm + r5plterm) - 3 * vmag / rmag * (Ptocb + Ptopl)

               F_T(:) = (Ftr + (Ptocb + Ptopl) * dot_product(v_unit, r_unit) / rmag) * r_unit(:)  &
                        + Ptopl * ((pl%rot(:,i) - theta_dot(:)) .cross. r_unit(:))  &
                        + Ptocb * ((cb%rot(:)   - theta_dot(:)) .cross. r_unit(:))
               cb%atide(:) = cb%atide(:) + F_T(:) / cb%Gmass
               pl%atide(:,i) = F_T(:) / pl%Gmass(i) 
            end do 

            do i = 1, npl
               pl%ah(:,i) = pl%ah(:,i) + pl%atide(:,i) + cb%atide(:)
            end do 
         end associate
      end select
      end select

      return
   end subroutine tides_kick_getacch_pl
   
end submodule s_tides_kick_getacch