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
function ringmoons_lindblad_torque(swifter_pl1P,ring,Gm,a,e,inc) result(Torque)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_lindblad_torque
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(inout)    :: ring
      real(DP),intent(in)                    :: Gm, a, e, inc
      real(DP)                               :: Torque
      

! Internals
      integer(I4B)                           :: j, m, inner_outer_sign
      integer(I4B), parameter                :: m_max = 3 ! Maximum mode number 
      real(DP)                               :: y, dTorque, beta, Amk

! Executable code

      ! Just do the first order resonances for now. The full suite of resonances will come later
     
      Torque = 0.0_DP
      do m  = 2, m_max
         !do inner Lindblad first
         do inner_outer_sign = -1,1,2
            y = (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(2._DP / 3._DP) * a   !resonance location for first order resonances
            j = ringmoons_ring_bin_finder(ring, y) !disk location of resonance
            if ((j == 0).or.(j == ring%N + 1)) cycle
            select case(inner_outer_sign)
            case(-1) 
               beta = ring%r(j) / a
            case(1)
               beta = a / ring%r(j)
            end select
            Amk = 0.5_DP * (2 * m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) + &
                           beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1))
            dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) * &
                      ring%Gsigma(j) * (ring%r(j)**2 * beta * ring%w(j) * Gm / swifter_pl1P%mass * Amk)**2
            ring%Torque(j) = ring%Torque(j) + dTorque
            Torque = Torque - dTorque
         end do
      end do

      return
end function ringmoons_lindblad_torque
