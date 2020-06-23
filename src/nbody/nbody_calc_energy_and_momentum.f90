submodule (swiftest_classes) s_nbody_calc_energy_and_momentum
contains
   module procedure nbody_calc_energy_and_momentum
      !! author: David A. Minton
      !!
      !! Compute total system angular momentum vector and kinetic, potential and total system energy
      !!
      !! Adapted from David E. Kaufmann's Swifter modules: symba_energy.f90
      !! Adapted from Martin Duncan's Swift routine anal_energy.f 
      use swiftest
      implicit none

      integer(I4B)               :: i, j
      real(DP)                   :: mass, msys, r2, v2, oblpot, ke, pe
      real(DP), dimension(NDIM)  :: h, x, v, dx, htot
      real(DP), dimension(:), allocatable   :: irh

      call self%h2b()
      htot(:) = 0.0_DP
      ke = 0.0_DP
      pe = 0.0_DP

      !$omp parallel do default(private) &
      !$omp shared (self) &
      !$omp reduction (+:ke, pe, htot)
      do i = 1, self%pl%nbody 
         x(:) = self%pl%xb(:,i)
         v(:) = self%pl%vb(:,i)
         mass = self%pl%mass(i)
         call util_crossproduct(x,v,h)
         htot(:) = htot(:) + mass * h(:)
         v2 = dot_product(v(:), v(:))
         ke = ke + 0.5_DP * mass * v2
         do j = i + 1, self%pl%nbody
            dx(:) = self%xb(:,j) - x(:)       
            r2 = dot_product(dx(:), dx(:))    
            if (r2 /= 0) then
               pe = pe - mass * self%pl%mass(j) / sqrt(r2) 
            end if
         end do
      end do
      !$omp end parallel do

      ! Add in the central body 
      x(:) = self%cb%xb(:,i)
      v(:) = self%cb%vb(:,i)
      mass = self%cb%mass(i)
      call util_crossproduct(x,v,h)
      htot(:) = htot(:) + mass * h(:)
      v2 = dot_product(v(:), v(:))
      ke = ke + 0.5_DP * mass * v2
      if (self%cb%j2rp2 /= 0.0_DP) then
         allocate(irh(self%pl%nbody))
         do i = 1, self%pl%nbody
            r2 = dot_product(self%pl%xh(:,i),self%pl%xh(:,i))
            irh(i) = 1.0_DP / sqrt(r2)
         end do
         call obl_pot(self%pl, self%cb%j2rp2, self%cb%j4rp4, self%pl%xh(:,:), irh, oblpot)
         deallocate(irh)
         pe = pe + oblpot
      end if

      self%pe = pe
      self%ke = ke
      self%te = ke + pe
      self%htot(:) = htot(:) 

      return

   end procedure nbody_calc_energy_and_momentum
end submodule s_nbody_calc_energy_and_momentum
