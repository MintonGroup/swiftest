submodule (symba) s_symba_discard_peri_pl
contains
   module procedure symba_discard_peri_pl
   !! author: David A. Minton
   !!
   !! Check to see if planets should be discarded based on their pericenter distances
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: symba_discard_peri_pl.f90
   !! Adapted from Hal Levison's Swift routine discard_mass_peri.f
use swiftest
implicit none
   logical , save      :: lfirst = .true.
   integer(I4B)          :: i, j, ih
   real(DP)            :: r2
   real(DP), dimension(NDIM) :: dx


! executable code
   if (lfirst) then
      call symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
      lfirst = .false.
   else
      call symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
      do i = 2, npl
         if (symba_plA%status(i) == ACTIVE) then
            if ((symba_plA%isperi(i) == 0) .and. (symba_plA%nplenc(i)== 0)) then
               if ((symba_plA%atp(i) >= qmin_alo) .and. (symba_plA%atp(i) <= qmin_ahi) &
                .and. (symba_plA%peri(i) <= qmin)) then
                  ldiscards = .true.
                  symba_plA%status(i) = DISCARDED_PERI
                  write(*, *) "particle ", symba_plA%name(i), &
                   " perihelion distance too small at t = ", t
               end if
            end if
         end if
      end do
   end if

   return

   end procedure symba_discard_peri_pl
end submodule s_symba_discard_peri_pl
