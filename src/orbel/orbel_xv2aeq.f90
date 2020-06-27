submodule (swiftest_classes) s_orbel_xv2aeq


contains

   module procedure orbel_xv2aeq_vec
      !! author: David A. Minton
      !!
      !! A wrapper method that converts all of the cartesian position and velocity vectors of a Swiftest body object to orbital elements.
      use swiftest
      implicit none
   
      call self%set_vec(cb)
      associate(n => self%nbody)
         call orbel_xv2aeq(self%mu_vec(1:n), self%xh(1:n, 1),  self%xh(1:n, 2),  self%xh(1:n, 3), &
                                          self%vh(1:n, 1),  self%vh(1:n, 2),  self%vh(1:n, 3), &
                                          self%a(1:n),     self%e(1:n),     self%inc(1:n),  &
                                          self%capom(1:n), self%omega(1:n), self%capm(1:n))
      end associate
      return
   end procedure orbel_xv2aeq_vec 

   module procedure orbel_xv2aeq
      !! author: David A. Minton
      !!
      !! Compute semimajor axis, eccentricity, and pericentric distance from relative Cartesian position and velocity
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2aeq.f90
      !! Adapted from Luke Dones' Swift routine orbel_xv2aeq.f
      use swiftest
      implicit none
      integer(I4B) :: iorbit_type
      real(DP)   :: r, v2, h2, energy, fac
      real(DP), dimension(NDIM) :: x, v, hvec

      a = 0.0_DP
      e = 0.0_DP
      q = 0.0_DP
      x = (/px, py, pz/)
      v = (/vx, vy, vz/)
      r = sqrt(dot_product(x(:), x(:)))
      v2 = dot_product(v(:), v(:))
      hvec(:) = x(:) .cross. v(:)
      h2 = dot_product(hvec(:), hvec(:))
      if (h2 == 0.0_DP) return
      energy = 0.5_DP * v2 - mu / r
      if (abs(energy * r / mu) < sqrt(VSMALL)) then
         iorbit_type = PARABOLA
      else
         a = -0.5_DP * mu / energy
         if (a < 0.0_DP) then
            fac = -h2 / (mu * a)
            if (fac > VSMALL) then
               iorbit_type = HYPERBOLA
            else
               iorbit_type = PARABOLA
            end if
         else
            iorbit_type = ELLIPSE
         end if
      end if
      select case (iorbit_type)
      case (ELLIPSE)
         fac = 1.0_DP - h2 / (mu * a)
         if (fac > VSMALL) e = sqrt(fac)
         q = a * (1.0_DP - e)
      case (PARABOLA)
         a = 0.5_DP * h2 / mu
         e = 1.0_DP
         q = a
      case (HYPERBOLA)
         e = sqrt(1.0_DP + fac)
         q = a * (1.0_DP - e)
      end select

      return

   end procedure orbel_xv2aeq
end submodule s_orbel_xv2aeq
