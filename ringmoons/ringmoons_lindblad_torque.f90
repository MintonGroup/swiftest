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
   real(DP),dimension(0:ring%N+1)         :: Torque,ring_Gsigma
   

! Internals
   integer(I4B)                           :: i,j, m, inner_outer_sign,w,w1,w2,js, mshep
   real(DP)                               :: a, dTorque, beta, Amk, width, nw,lap,dlap,da3,Xs,Xlo,Xhi
   real(DP), parameter                    :: g = 2.24_DP
   real(DP)                               :: const_fac
   real(DP),dimension(M_MAX)              :: X_I,X_O,Xw_I,Xw_O,nw_I,nw_O,beta_I,beta_O,Amk_I,Amk_O
   integer(I4B),dimension(M_MAX)          :: w1_I,w2_I,w1_O,w2_O
   logical(lgt),save                      :: firstrun = .true.
   real(DP),dimension(M_MAX),save         :: marr


! Executable code

   if (firstrun) then
      do m = 2, M_MAX
         marr(m) = real(m, kind=DP)
      end do
      firstrun = .false.
   end if
      
     
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque(:) = 0.0_DP
   Xs = 2 * sqrt(as)
   Xlo = ring%X_I + ring%deltaX * ring%inside  
   Xhi = ring%X_F
   mshep = min(M_MAX,ceiling(0.5_DP * (sqrt(1._DP + (4._DP / 3._DP) * (Xs / ring%deltaX)) - 1._DP)))
   
   where (ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk)
      ring_Gsigma(:) = ring%Gsigma(:)
   elsewhere
      ring_Gsigma(:) = 0.0_DP
   end where 
   const_fac = 0.25_DP *  PI**2 / 3._DP * (Gm / swifter_pl1P%mass)**2

   ! Resonance locations
   X_I(2:mshep) = Xs * (1._DP - 1._DP / marr(2:mshep))**(1._DP / 3._DP)
   X_O(2:mshep) = Xs * (1._DP + 1._DP / marr(2:mshep))**(1._DP / 3._DP)

   ! Resonance widths
   Xw_I(2:mshep) = X_I(2:mshep) * (Gm / swifter_pl1P%mass)**(0.25_DP)
   Xw_O(2:mshep) = X_O(2:mshep) * (Gm / swifter_pl1P%mass)**(0.25_DP)


   ! Resonance betas
   beta_I(2:mshep) = (X_I(2:mshep) / Xs)**2
   beta_O(2:mshep) = (Xs / X_O(2:mshep))**2

   ! Laplace coefficient terms
   Amk_I(2:mshep) = lapm(-1,2:mshep) + dlapm(-1,2:mshep)
   Amk_O(2:mshep) = lapm( 1,2:mshep) + dlapm( 1,2:mshep)

   ! Bins where the resonances occur
   w1_I(2:mshep) = min(max(ceiling((X_I(2:mshep) - Xw_I(2:mshep) - ring%X_I) / ring%deltaX),0),ring%N+1)
   w2_I(2:mshep) = min(max(ceiling((X_I(2:mshep) + Xw_I(2:mshep) - ring%X_I) / ring%deltaX),0),ring%N+1)
   w1_O(2:mshep) = min(max(ceiling((X_O(2:mshep) - Xw_O(2:mshep) - ring%X_I) / ring%deltaX),0),ring%N+1)
   w2_O(2:mshep) = min(max(ceiling((X_O(2:mshep) + Xw_O(2:mshep) - ring%X_I) / ring%deltaX),0),ring%N+1)

   ! Number of bins for each resonance (used to distribute the torque over each bin)
   nw_I(2:mshep) = real(w2_I(2:mshep) - w2_I(2:mshep) + 1, kind=DP)
   nw_O(2:mshep) = real(w2_O(2:mshep) - w2_O(2:mshep) + 1, kind=DP)

   do m = 2,mshep

      ! Inner Lindblad  
      if ((X_I(m) < Xhi) .and. (X_I(m) > Xlo)) then
         w1 = w1_I(m)
         w2 = w2_I(m)
         Torque(w1:w2) = Torque(w1:w2) - const_fac * (m / real(m  - 1,kind=DP)) / nw_I(m)  &
                                      * ring_Gsigma(w1:w2) * (X_I(m)**4 * beta_I(m) * ring%w(w1:w2) * Amk_I(m))**2  
      end if
      
      ! Outer Lindblad
      if ((X_O(m) < Xhi) .and. (X_O(m) > Xlo)) then
         w1 = w1_O(m)
         w2 = w2_O(m)
         Torque(w1:w2) = Torque(w1:w2) + const_fac * (m / real(m  - 1,kind=DP))  / nw_O(m)  &
                                   * ring_Gsigma(w1:w2) * (X_O(m)**4 * beta_O(m) * ring%w(w1:w2) * Amk_I(m))**2  
      end if
   end do



   ! Add in shepherding torque
   do inner_outer_sign = -1,1,2
      !write(*,*) inner_outer_sign,'shep'
      a = (1._DP + inner_outer_sign * 1._DP / real(mshep + 1, kind=DP))**(2._DP / 3._DP) * as   
      j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
      if ((j > ring%inside).and.(j < ring%N + 1)) then 
         da3 = inner_outer_sign * max(abs((a - as)**3),epsilon(a))
         dTorque = g**2 / 6._DP * a**3 / da3 * (Gm / swifter_pl1P%mass)**2 * ring_Gsigma(j) * (ring%w(j))**2 * a**4
         Torque(j) = Torque(j) + dTorque
      end if
   end do

   return
end function ringmoons_lindblad_torque
