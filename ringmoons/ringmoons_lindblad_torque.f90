!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_lindblad
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the lindblad torques between each ring element and a given satellite. Function returns total torque on the
!                satellite, and stores the torques acting on each ring element in the ring
!
!  Input
!    Arguments : 
!                
!    Terminal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : Torque = ringmoons_lindblad_torque(swifter_pl1P,ring,Gm,a,e,inc)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_lindblad_torque(swifter_pl1P,ring,Gm,as,e,inc) result(Torque)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_lindblad_torque
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(in)       :: ring
   real(DP),intent(in)                    :: Gm, as, e, inc
   real(DP),dimension(0:ring%N+1)         :: Torque
   

! Internals
   integer(I4B)                           :: i,j, m, il,w,w1,w2,js, mshep
   real(DP)                               :: a, dTorque, beta, Amk, nw,lap,dlap,da3,Xs, Xlo,Xhi,Gfac,lind_factor,Xw2
   real(DP), parameter                    :: g = 2.24_DP
   real(DP),dimension(2:M_MAX)            :: Xr,ar
   logical(lgt),dimension(0:ring%N+1)     :: T_mask
   integer(I4B),dimension(2:M_MAX)        :: w1_arr,w2_arr


! Executable code
   Gfac = (Gm / swifter_pl1P%mass)

   ! Mask out any ring bins that don't have enough mass in them
   where ((ring%Gsigma(0:ring%N+1) * MU2GM / DU2CM**2 / GU) > 1e-3_DP)
      T_mask(0:ring%N+1) = .true.
   elsewhere
      T_mask(0:ring%N+1) = .false. 
   end where
   Xs = 2 * sqrt(as)
   Xw2 = 0.5_DP * Xs**2 * sqrt(Gfac)
   Xlo = ring%X_I + ring%deltaX * ring%inside
   Xhi = ring%X_F
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque(:) = 0.0_DP

   ! Calculate the number of modes that are separated by at least 1 bin width in X space
   mshep = max(2,min(M_MAX - 1,ceiling(0.5_DP * (sqrt(1._DP + 4._DP / 3._DP * Xs / ring%deltaX) - 1._DP))))
   
   ! Inner then outer Lindblads
   do il = -1,1,2
      ! Calculate resonance location in X space
      Xr(2:mshep) = Xs * marr(2:mshep,il)

      ! Calculate bin boundaries of resonance using its width in X space
      where((Xr(2:mshep) > Xlo).and.(Xr(2:mshep) < Xhi))
         w1_arr(2:mshep) = min(max(ceiling((sqrt(Xr(2:mshep)**2 - Xw2) - ring%X_I) / ring%deltaX),0),ring%N+1)
         w2_arr(2:mshep) = min(max(ceiling((sqrt(Xr(2:mshep)**2 + Xw2) - ring%X_I) / ring%deltaX),0),ring%N+1)
      elsewhere ! The resonance is outside the bins
         w1_arr(2:mshep) = 0
         w2_arr(2:mshep) = 0
      end where 
      do m  = 2, mshep
         if ((Xr(m) > Xlo).and.(Xr(m) < Xhi)) then
            beta = (Xs / Xr(m))**(il * 2)
            a = 0.25_DP * (Xr(m))**2
            lap  =  lapm(m,il)
            dlap = dlapm(m,il)
            Amk = (lap + dlap)
            w1 = w1_arr(m)
            w2 = w2_arr(m)
            nw = real(w2 - w1 + 1,kind=DP)
            ! Calculate the 1st order Lindblad torques and distribute them over the bins that include the resonance
            lind_factor = il * mfac(m) / nw * a**4 * (beta * Gfac  * Amk)**2 
            where(T_mask(w1:w2)) Torque(w1:w2) = Torque(w1:w2) + lind_factor * ring%Gsigma(w1:w2) * (ring%w(w1:w2))**2
         end if
      end do

      ! Add in shepherding torque
      a = as * marr(mshep+1,il)**2
      j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
      if ((j > ring%inside).and.(j < ring%N + 1)) then
         da3 = il * max(abs((a - as)**3),epsilon(a))
         if (T_mask(j)) Torque(j) = Torque(j) + g**2 / 6._DP * a**3 / da3 * (Gfac)**2 * ring%Gsigma(j) * (ring%w(j))**2 * a**4
      end if
   end do

   return
end function ringmoons_lindblad_torque
