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
   mu = symba_pla%helio%swiftest%mass(1)
!$omp parallel do default(none) &
!$omp shared (symba_pla, npl, mu, dt, irec) &
!$omp private (i, iflag)
   do i = 2, npl
      if ((symba_pla%levelg(i) == irec) .and. (symba_pla%helio%swiftest%status(i) == active)) then
         call drift_one(mu, symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%vb(:,i), dt, iflag)
         if (iflag /= 0) then
            write(*, *) " massive body ", symba_pla%helio%swiftest%name(i), " is lost!!!!!!!!!!"
            write(*, *) mu, dt
            write(*, *) symba_pla%helio%swiftest%xh(:,i)
            write(*, *) symba_pla%helio%swiftest%vb(:,i)
            write(*, *) " stopping "
            call util_exit(failure)
         end if
      end if
   end do
!$omp end parallel do

   return

   end procedure symba_helio_drift
end submodule s_symba_helio_drift
