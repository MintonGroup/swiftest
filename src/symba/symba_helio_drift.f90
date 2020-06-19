submodule (symba) s_symba_helio_drift
contains
   module procedure symba_helio_drift
   !! author: David A. Minton
   !!
   !! Loop through planets and call Danby drift routine
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_helio_drift.f90
   !! Adapted from Hal Levison's Swift routine symba5_helio_drift.f
use swiftest
implicit none
   integer(I4B)          :: i, iflag
   real(DP)            :: mu

! executable code
   mu = symba_plA%mass(1)
!$omp parallel do default(none) &
!$omp shared (symba_plA, npl, mu, dt, irec) &
!$omp private (i, iflag)
   do i = 2, npl
      if ((symba_plA%levelg(i) == irec) .and. (symba_plA%status(i) == ACTIVE)) then
         call drift_one(mu, symba_plA%xh(:,i), symba_plA%vb(:,i), dt, iflag)
         if (iflag /= 0) then
            write(*, *) " massive body ", symba_plA%name(i), " is lost!!!!!!!!!!"
            write(*, *) mu, dt
            write(*, *) symba_plA%xh(:,i)
            write(*, *) symba_plA%vb(:,i)
            write(*, *) " stopping "
            call util_exit(failure)
         end if
      end if
   end do
!$omp end parallel do

   return

   end procedure symba_helio_drift
end submodule s_symba_helio_drift
