submodule (swiftest_data_structures) s_drift_dan
contains
   module procedure drift_dan
   !! author: David A. Minton
   !!
   !! Perform Kepler drift, solving Kepler's equation in appropriate variables
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: drift_dan.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift routine drift_dan.f
   use swiftest
   implicit none
   real(DP)            :: dt, f, g, fdot, gdot, c1, c2, c3, u, alpha, fp, r0
   real(DP)            :: v0s, a, asq, en, dm, ec, es, esq, xkep, fchk, s, c
   real(DP), dimension(ndim) :: x, v

! executable code
   iflag = 0
   dt = dt0
   r0 = sqrt(dot_product(x0(:), x0(:)))
   v0s = dot_product(v0(:), v0(:))
   u = dot_product(x0(:), v0(:))
   alpha = 2.0_DP*mu/r0 - v0s
   if (alpha > 0.0_DP) then
      a = mu/alpha
      asq = a*a
      en = sqrt(mu/(a*asq))
      ec = 1.0_DP - r0/a
      es = u/(en*asq)
      esq = ec*ec + es*es
      dm = dt*en - int(dt*en/twopi)*twopi
      dt = dm/en
      if ((esq < e2max) .and. (dm*dm < dm2max) .and. (esq*dm*dm < e2dm2max)) then
         call drift_kepmd(dm, es, ec, xkep, s, c)
         fchk = (xkep - ec*s + es*(1.0_DP - c) - dm)
! dek - original code compared fchk*fchk with danbyb, but i think it should
! dek - be compared with danbyb*danbyb, and i changed it accordingly - please
! dek - check with hal and/or martin about this
         if (fchk*fchk > danbyb*danbyb) then
            iflag = 1
            return
         end if
         fp = 1.0_DP - ec*c + es*s
         f = a/r0*(c - 1.0_DP) + 1.0_DP
         g = dt + (s - xkep)/en
         fdot = -(a/(r0*fp))*en*s
         gdot = (c - 1.0_DP)/fp + 1.0_DP
         x(:) = x0(:)*f + v0(:)*g
         v(:) = x0(:)*fdot + v0(:)*gdot
         x0(:) = x(:)
         v0(:) = v(:)
         iflag = 0
         return
      end if
   end if
   call drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
   if (iflag == 0) then
      f = 1.0_DP - mu/r0*c2
      g = dt - mu*c3
      fdot = -mu/(fp*r0)*c1
      gdot = 1.0_DP - mu/fp*c2
      x(:) = x0(:)*f + v0(:)*g
      v(:) = x0(:)*fdot + v0(:)*gdot
      x0(:) = x(:)
      v0(:) = v(:)
   end if

   return

   end procedure drift_dan
end submodule s_drift_dan
