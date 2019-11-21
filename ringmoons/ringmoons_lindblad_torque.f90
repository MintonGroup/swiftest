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
   real(DP)                               :: a, dTorque, beta, Amk, width, nw,lap,dlap,da3,Xs, Xlo,Xhi,Gfac
   real(DP), parameter                    :: g = 2.24_DP
   logical(lgt),save                      :: firstrun = .true.
   real(DP),dimension(M_MAX,-1:1),save    :: marr
   real(DP),dimension(M_MAX),save         :: mfac
   real(DP),dimension(M_MAX,-1:1)         :: Xr,Xw
   logical(lgt),dimension(0:ring%N+1)      :: T_mask


! Executable code
   if (firstrun) then
      do m = 2,M_MAX
         marr(m,-1) = (1._DP - 1._DP / real(m, kind=DP))**(1._DP / 3._DP)
         marr(m, 1) = (1._DP + 1._DP / real(m, kind=DP))**(1._DP / 3._DP)
         mfac(m) = 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) 
      end do
      firstrun = .false.
   end if 
   Gfac = (Gm / swifter_pl1P%mass)

   ! Mask out any ring bins that don't have enough mass in them
   where (ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk)
      T_mask(:) = .true.
   elsewhere
      T_mask(:) = .false. 
   end where
   Xs = 2 * sqrt(as)
   Xlo = ring%X_I + ring%deltaX * ring%inside
   Xhi = ring%X_F
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque(:) = 0.0_DP
   mshep = min(M_MAX,ceiling(0.5_DP * (sqrt(1._DP + 8._DP / 3._DP * sqrt(as) / ring%deltaX) - 1._DP)))
   
   ! Inner then outer lindblads
   do il = -1,1,2
      Xr(2:mshep,il) = Xs * marr(2:mshep,il)
      Xw(2:mshep,il) = Xr(2:mshep,il) * (Gfac)**(0.25_DP)
   
      do m  = 2, mshep
         if ((Xr(m,il) > Xlo).and.(Xr(m,il) < Xhi)) then
            beta = (Xs / Xr(m,il))**(il * 2)
            a = 0.25_DP * (Xr(m,il))**2
            lap  =  lapm(m,il)
            dlap = dlapm(m,il)
            Amk = (lap + dlap)
            w1 = min(max(ceiling((Xr(m,il) - Xw(m,il) - ring%X_I) / ring%deltaX),0),ring%N+1)
            w2 = min(max(ceiling((Xr(m,il) + Xw(m,il) - ring%X_I) / ring%deltaX),0),ring%N+1)
            nw = real(w2 - w1 + 1,kind=DP)
            where(T_mask(w1:w2)) Torque(w1:w2) = Torque(w1:w2) + il * mfac(m) / nw * ring%Gsigma(w1:w2) * (a**2 * beta * ring%w(w1:w2) * Gfac  * Amk)**2 
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
