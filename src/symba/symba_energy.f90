submodule (symba) s_symba_energy
contains
   module procedure symba_energy
   !! author: David A. Minton
   !!
   !! Compute total system angular momentum vector and kinetic, potential and total system energy
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_energy.f90
   !! Adapted from Martin Duncan's Swift routine anal_energy.f
use swiftest
implicit none
   logical(lgt)                     :: lmalloc = .true.
   integer(I4B)                     :: i, j
   real(DP)                       :: mass, msys, r2, v2, oblpot
   real(DP), dimension(ndim)            :: h, x, v, dx
   real(DP), dimension(npl)             :: irh


! executable code

   call coord_h2b(npl, swiftest_plA, msys)
   htot = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
   ke = 0.0_DP
   pe = 0.0_DP

!$omp parallel do default(none) &
!$omp shared (swiftest_plA, npl) &
!$omp private (i, x, v, mass, h, htot, v2, dx, r2) &
!$omp reduction (+:ke, pe)
   do i = 1, npl - 1
      x(:) = swiftest_plA%xb(:,i)
      v(:) = swiftest_plA%vb(:,i)
      mass = swiftest_plA%mass(i)
      h(1) = mass*(x(2)*v(3) - x(3)*v(2))
      h(2) = mass*(x(3)*v(1) - x(1)*v(3))
      h(3) = mass*(x(1)*v(2) - x(2)*v(1))
      htot(:) = htot(:) + h(:)
      v2 = dot_product(v(:), v(:))
      ke = ke + 0.5_DP*mass*v2
      do j = i + 1, npl
         dx(:) = swiftest_plA%xb(:,j) - x(:) !this is 0 for the removed ps because swiftest_plA%xb(:,j) = swiftest_plA%xb(:,i)
         r2 = dot_product(dx(:), dx(:)) !this is 0 for the removed ps
         if (r2 /= 0) then
            pe = pe - mass*swiftest_plA%mass(j)/sqrt(r2) !division !!!!!!r2 is 0 for ps 12 which is the removed ps
         end if
      end do
   end do
!$omp end parallel do
   i = npl ! needed to account for the parllelization above
   x(:) = swiftest_plA%xb(:,i)
   v(:) = swiftest_plA%vb(:,i)
   mass = swiftest_plA%mass(i)
   h(1) = mass*(x(2)*v(3) - x(3)*v(2))
   h(2) = mass*(x(3)*v(1) - x(1)*v(3))
   h(3) = mass*(x(1)*v(2) - x(2)*v(1))
   htot(:) = htot(:) + h(:)
   v2 = dot_product(v(:), v(:))
   ke = ke + 0.5_DP*mass*v2
   if (config%j2rp2 /= 0.0_DP) then
      do i = 2, npl
         r2 = dot_product(swiftest_plA%xh(:,i),swiftest_plA%xh(:,i))
         irh(i) = 1.0_DP/sqrt(r2)
      end do
      call obl_pot(swiftest_plA, config%j2rp2, config%j4rp4, swiftest_plA%xh(:,:), irh, oblpot)
      pe = pe + oblpot
   end if
   te = ke + pe

   return

   end procedure symba_energy
end submodule s_symba_energy
