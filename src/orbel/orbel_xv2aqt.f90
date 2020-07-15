submodule (swiftest_classes) s_orbel_xv2aqt
contains
   module procedure orbel_xv2aqt ! (mu, px, py, pz, vx, vy, vz, a, q, capm, tperi
      !! author: David A. Minton
      !!
      !! Compute semimajor axis, pericentric distance, mean anomaly, and time to nearest pericenter passage from
      !! relative Cartesian position and velocity
      !!      tperi > 0 means nearest pericenter passage is in the future
      !!      tperi < 0 means nearest pericenter passage is in the past
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: orbel_xv2aqt.f90
      use swiftest
      implicit none
      integer(I4B) :: iorbit_type
      real(DP)   :: r, v2, h2, rdotv, energy, fac, w, face, cape, e, tmpf, capf, mm
      real(DP), dimension(NDIM) :: hvec

      a = 0.0_DP
      q = 0.0_DP
      capm = 0.0_DP
      tperi = 0.0_DP
      r = norm2(x(:))
      v2 = dot_product(v(:), v(:))
      hvec(:) = x(:) .cross. v(:)
      h2 = dot_product(hvec(:), hvec(:))
      if (h2 == 0.0_DP) return
      rdotv = dot_product(x(:), v(:))
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
         if (fac > VSMALL) then
            e = sqrt(fac)
            cape = 0.0_DP
            face = (a - r) / (a * e)
            if (face < -1.0_DP) then
               cape = PI
            else if (face < 1.0_DP) then
               cape = acos(face)
            end if
            if (rdotv < 0.0_DP) cape = TWOPI - cape
         else
            e = 0.0_DP
            cape = 0.0_DP
         end if
         capm = cape - e * sin(cape)
         q = a * (1.0_DP - e)
         mm = sqrt(mu / a**3)
         if (capm < PI) then
            tperi = -1.0_DP * capm / mm
         else
            tperi = -1.0_DP * (capm - TWOPI) / mm
         end if
      case (PARABOLA)
         a = 0.5_DP * h2 / mu
         e = 1.0_DP
         w = 0.0_DP
         fac = 2 * a / r - 1.0_DP
         if (fac < -1.0_DP) then
            w = PI
         else if (fac < 1.0_DP) then
            w = acos(fac)
         end if
         if (rdotv < 0.0_DP) w = TWOPI - w
         tmpf = tan(0.5_DP * w)
         capm = tmpf*(1.0_DP + tmpf * tmpf / 3.0_DP)
         q = a
         mm = sqrt(0.5_DP * mu / q**3)
         tperi = -1.0_DP * capm / mm
      case (HYPERBOLA)
         e = sqrt(1.0_DP + fac)
         tmpf = (a - r) / (a * e)
         if (tmpf < 1.0_DP) tmpf = 1.0_DP
         capf = log(tmpf + sqrt(tmpf * tmpf - 1.0_DP))
         if (rdotv < 0.0_DP) capf = -capf
         capm = e * sinh(capf) - capf
         q = a * (1.0_DP - e)
         mm = sqrt(-mu / a**3)
         tperi = -1.0_DP * capm / mm
      end select

      return

   end procedure orbel_xv2aqt
end submodule s_orbel_xv2aqt
