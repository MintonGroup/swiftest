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
   type(ringmoons_ring), intent(inout)    :: ring
   real(DP),intent(in)                    :: Gm, as, e, inc
   real(DP)                               :: Torque
   

! Internals
   integer(I4B)                           :: i,j, m, inner_outer_sign,w,w1,w2, j0
   integer(I4B), parameter                :: m_max = 6 ! Maximum mode number 
   real(DP)                               :: a, dTorque, beta, Amk, width, nw,lap,dlap,sigavg
   logical(lgt), save                     :: first_run = .true.
   real(DP), dimension(-1:1,2:m_max), save :: lapm,dlapm


! Executable code


   ! For performance reasons, we compute a table of Laplace coefficient terms the first time through and then interpolate 
   if (first_run) then
      do m = 2, m_max
         do inner_outer_sign = -1,1,2
         !do j = 1, NLAP
            beta =  (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(-inner_outer_sign * 2._DP / 3._DP)
            lapm(inner_outer_sign,m)  = m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) 
            dlapm(inner_outer_sign,m) = 0.5_DP * beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1) 
         end do
      end do
      first_run  = .false.
   end if
      
  
   ! Just do the first order resonances for now. The full suite of resonances will come later
   Torque = 0.0_DP
   do m  = 2, m_max
      !do inner Lindblad first
      j0 = -1
      do inner_outer_sign = -1,1,2
         a = (1._DP + inner_outer_sign * 1._DP / real(m, kind=DP))**(2._DP / 3._DP) * as   !resonance location for first order resonances
         j = ringmoons_ring_bin_finder(ring, a) !disk location of resonance
         if ((j > 0).and.(j < ring%N + 1)) then !.and.(j /= j0)) then ! Resonance is in the bin, and don't consider overlapping
            select case(inner_outer_sign)
            case(-1) 
               beta = a / as
            case(1)
               beta = as / a
            end select
            lap  =  lapm(inner_outer_sign,m)
            dlap = dlapm(inner_outer_sign,m)

            Amk = (lap + dlap)
            width = sqrt(Gm / swifter_pl1P%mass) * a !ring%r(j)
            w1 = ringmoons_ring_bin_finder(ring,a - width)
            w2 = ringmoons_ring_bin_finder(ring,a + width)
            nw = real(w2 - w1 + 1,kind=DP)
            sigavg = sum(ring%Gsigma(w1:w2)) / nw
            do w = w1,w2 
               dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) / nw * &
                      ring%Gsigma(w) * (a**2 * beta * sqrt(swifter_pl1P%mass / a**3) * Gm / swifter_pl1P%mass * Amk)**2 
               ring%Torque(w) = ring%Torque(w) + dTorque 
               Torque = Torque - dTorque
            end do
         end if
         j0 = j
      end do
   end do

   return
end function ringmoons_lindblad_torque
