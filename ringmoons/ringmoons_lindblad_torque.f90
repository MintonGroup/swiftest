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
   integer(I4B)                           :: i,j, m, inner_outer_sign,w,w1,w2,js, mshep
   real(DP)                               :: a, dTorque, beta, Amk, width, nw,lap,dlap,da3
   real(DP), parameter                    :: g = 2.24_DP


! Executable code



     
   js = ringmoons_ring_bin_finder(ring,as) 
  
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque(:) = 0.0_DP
   mshep = min(M_MAX,ceiling(0.5_DP * (sqrt(1._DP + 8._DP / 3._DP * sqrt(as) / ring%deltaX) - 1._DP)))
   do m  = 2, mshep
      ! Go through modes up until resonance overlap occurs
      do inner_outer_sign = -1,1,2
         a = (1._DP + inner_outer_sign * 1._DP / real(m, kind=DP))**(2._DP / 3._DP) * as   !resonance location for first order resonances
         j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
         if ((j > 0).and.(j < ring%N + 1)) then !Resonance is in the bin, and don't consider overlapping
            beta = (as / a)**(inner_outer_sign)
            lap  =  lapm(inner_outer_sign,m)
            dlap = dlapm(inner_outer_sign,m)

            Amk = (lap + dlap)
            width = sqrt(Gm / swifter_pl1P%mass) * a 
            w1 = ringmoons_ring_bin_finder(ring,a - width)
            w2 = ringmoons_ring_bin_finder(ring,a + width)
            nw = real(w2 - w1 + 1,kind=DP)
            do w = w1,w2 
               dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) / nw * &
                      ring%Gsigma(w) * (a**2 * beta * ring%w(i) * (Gm / swifter_pl1P%mass) * Amk)**2 
               
               if ((abs(dTorque) > huge(dTorque)).or.(dTorque /= dTorque)) then
                  write(*,*) 'Bad Lindblad torque'
                  write(*,*) m,j,w,beta,a,as
                  write(*,*) ring%Gsigma(w),dTorque
                  call util_exit(FAILURE)
               end if
               Torque(w) = Torque(w) + dTorque
            end do

         end if
      end do
   end do   
   ! Add in shepherding torque
   do inner_outer_sign = -1,1,2
      a = (1._DP + inner_outer_sign * 1._DP / real(mshep + 1, kind=DP))**(2._DP / 3._DP) * as   
      j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
      if ((j > 0).and.(j < ring%N + 1)) then 
         da3 = inner_outer_sign * max(abs((a - as)**3),epsilon(a))
         dTorque = g**2 / 6._DP * a**3 / da3 * (Gm / swifter_pl1P%mass)**2 * ring%Gsigma(j) * (ring%w(j))**2 * a**4
         if ((abs(dTorque) > huge(dTorque)).or.(dTorque /= dTorque)) then
            write(*,*) 'Bad sheperding torque: m = ',mshep
            write(*,*) 'mcalc = ',(0.5_DP * (sqrt(1._DP + 8._DP / 3._DP * sqrt(as) / ring%deltaX) - 1._DP))
            write(*,*) a,as,da3
            write(*,*) ring%Gsigma(j),dTorque
            call util_exit(FAILURE)
         end if
         Torque(j) = Torque(j) + dTorque
      end if
   end do

   return
end function ringmoons_lindblad_torque
