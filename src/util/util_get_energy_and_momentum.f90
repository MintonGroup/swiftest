submodule (swiftest_classes) s_util_get_energy_and_momentum
contains
   module procedure util_get_energy_and_momentum
      !! author: David A. Minton
      !!
      !! Compute total system angular momentum vector and kinetic, potential and total system energy
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_energy.f90
      !! Adapted from Martin Duncan's Swift routine anal_energy.f 
      use swiftest
      implicit none

      integer(I4B)               :: i, j
      real(DP)                   :: mass, msys, r2, v2, oblpot, ke, pe
      real(DP), dimension(NDIM)  :: h, x, v, dx, htot
      real(DP), dimension(:), allocatable   :: irh

      select type(pl => self%pl)
      class is (swiftest_pl)
         call pl%h2b(self%cb)
         htot(:) = 0.0_DP
         ke = 0.0_DP
         pe = 0.0_DP

         !$omp parallel do default(private) &
         !$omp shared (self) &
         !$omp reduction (+:ke, pe, htot)
         do i = 1, pl%nbody 
            x(:) = pl%xb(:, i)
            v(:) = pl%vb(:, i)
            mass = pl%Gmass(i)
            h(:) = x(:) .cross. v(:)
            htot(:) = htot(:) + mass * h(:)
            v2 = dot_product(v(:), v(:))
            ke = ke + 0.5_DP * mass * v2
            do j = i + 1, pl%nbody
               dx(:) = pl%xb(:, j) - x(:)       
               r2 = dot_product(dx(:), dx(:))
               if (r2 /= 0) then
                  pe = pe - mass * pl%Gmass(j) / sqrt(r2) 
               end if
            end do
         end do
         !$omp end parallel do
      end select

      select type(cb => self%cb)
      class is (swiftest_cb)
         ! Add in the central body 
         x(:) = cb%xb(:)
         v(:) = cb%vb(:)
         mass = cb%Gmass
         h(:) = x(:) .cross. v(:)
         htot(:) = htot(:) + mass * h(:)
         v2 = dot_product(v(:), v(:))
         ke = ke + 0.5_DP * mass * v2
         if (cb%j2rp2 /= 0.0_DP) then
            allocate(irh(self%pl%nbody))
            do i = 1, self%pl%nbody
               r2 = dot_product(self%pl%xh(:, i), self%pl%xh(:, i)) 
               irh(i) = 1.0_DP / sqrt(r2)
            end do
            oblpot = self%pl%obl_pot(cb, irh)
            deallocate(irh)
            pe = pe + oblpot
         end if
      end select
   
      self%pe = pe
      self%ke = ke
      self%te = ke + pe
      self%htot(:) = htot(:) 

      return

   end procedure util_get_energy_and_momentum
end submodule s_util_get_energy_and_momentum
