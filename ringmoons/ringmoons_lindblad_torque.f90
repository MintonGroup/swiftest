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
   integer(I4B)                           :: i,j,j1, m, inner_outer_sign,w,w1,w2,js
   integer(I4B), parameter                :: m_max = 100 ! Maximum mode number 
   real(DP)                               :: a, dTorque, beta, Amk, width, nw,lap,dlap,da,a1
   logical(lgt), save                     :: first_run = .true.
   real(DP), dimension(-1:1,2:m_max), save :: lapm,dlapm
   real(DP), parameter                    :: g = 2.24_DP
   logical(lgt)                           :: shepflag


! Executable code


   ! For performance reasons, we compute a table of Laplace coefficient terms the first time through and then interpolate 
   if (first_run) then
      do m = 2, m_max
         do inner_outer_sign = -1,1,2
            beta =  (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(-inner_outer_sign * 2._DP / 3._DP)
            lapm(inner_outer_sign,m)  = m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) 
            dlapm(inner_outer_sign,m) = 0.5_DP * beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1) 
         end do
      end do
      first_run  = .false.
   end if
     
   js = ringmoons_ring_bin_finder(ring,as) 
  
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque(:) = 0.0_DP
   shepflag = .false.
   do m  = 2, m_max
      ! Go through modes up until resonance overlap occurs
      if (shepflag) exit
      do inner_outer_sign = -1,1,2
         a = (1._DP + inner_outer_sign * 1._DP / real(m, kind=DP))**(2._DP / 3._DP) * as   !resonance location for first order resonances
         j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
         if ((j > 0).and.(j < ring%N + 1)) then !Resonance is in the bin, and don't consider overlapping
            a1 = (1._DP + inner_outer_sign * 1._DP / real(m + 1, kind=DP))**(2._DP / 3._DP) * as ! Check for overlapping
            j1 = ringmoons_ring_bin_finder(ring, a1) 
            if ((j1 /= j).and.(.not.shepflag)) then ! Not overlapping within a bin
               select case(inner_outer_sign)
               case(-1) 
                  beta = a / as
               case(1)
                  beta = as / a
               end select
               lap  =  lapm(inner_outer_sign,m)
               dlap = dlapm(inner_outer_sign,m)

               Amk = (lap + dlap)
               width = sqrt(Gm / swifter_pl1P%mass) * a 
               w1 = ringmoons_ring_bin_finder(ring,a - width)
               w2 = ringmoons_ring_bin_finder(ring,a + width)
               nw = real(w2 - w1 + 1,kind=DP)
               do w = w1,w2 
                  dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) / nw * &
                         ring%Gsigma(w) * (a**2 * beta * sqrt(swifter_pl1P%mass / a**3) * Gm / swifter_pl1P%mass * Amk)**2 
                  Torque(w) = Torque(w) + dTorque
               end do
            else ! Resonances overlap inside a bin. Model this as shepherding torque
               shepflag = .true.
               da = a - as
               dTorque = g**2 / 6._DP * (a / da)**3 * (Gm / swifter_pl1P%mass)**2 * ring%Gsigma(j) * ring%w(j)**2 * a**4
               Torque(j) = Torque(j) + dTorque
            end if
         end if
      end do
   end do 

   return
end function ringmoons_lindblad_torque
