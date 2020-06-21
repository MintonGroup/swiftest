submodule (swiftest_classes) s_nbody_calc_conserved
contains
   module procedure nbody_calc_conserved
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
   do i = 1, self%npl - 1
      x(:) = self%xb(:,i)
      v(:) = self%vb(:,i)
      mass = self%mass(i)
      call util_crossproduct(x,v,h)
      htot(:) = htot(:) + mass * h(:)
      v2 = dot_product(v(:), v(:))
      ke = ke + 0.5_DP * mass * v2
      do j = i + 1, self%npl
         dx(:) = self%xb(:,j) - x(:)       
         r2 = dot_product(dx(:), dx(:))    
         if (r2 /= 0) then
            pe = pe - mass * self%mass(j) / sqrt(r2) 
         end if
      end do
   end do
   !$omp end parallel do
   i = self%npl ! needed to account for the parllelization above
   x(:) = self%xb(:,i)
   v(:) = self%vb(:,i)
   mass = self%mass(i)
   call util_crossproduct(x,v,h)
   htot(:) = htot(:) + mass * h(:)
   v2 = dot_product(v(:), v(:))
   ke = ke + 0.5_DP * mass * v2
   if (config%j2rp2 /= 0.0_DP) then
      allocate(irh(self%npl))
      do i = 2, self%npl
         r2 = dot_product(self%xh(:,i),self%xh(:,i))
         irh(i) = 1.0_DP / sqrt(r2)
      end do
      call obl_pot(self, config%j2rp2, config%j4rp4, self%xh(:,:), irh, oblpot)
      deallocate(irh)
      pe = pe + oblpot
   end if
  
   self%pe = pe
   self%ke = ke
   self%te = ke + pe
   self%htot(:) = htot(:) 

   return

   end procedure nbody_calc_conserved
end submodule s_nbody_calc_conserved
