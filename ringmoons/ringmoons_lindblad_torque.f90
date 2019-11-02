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
   integer(I4B)                           :: i,j, m, inner_outer_sign,w,w1,w2
   integer(I4B), parameter                :: m_max = 3 ! Maximum mode number 
   real(DP)                               :: y, dTorque, beta, Amk, width, nw,lap,dlap,rem, bdb
   logical(lgt), save                     :: first_run = .true.
   integer(I4B), parameter                :: NLAP = 10000 ! Number of laplace coefficient distance ratios to pre-compute
   real(DP), dimension(-1:1,2:m_max), save :: lapm,dlapm
   !real(DP), parameter                    :: dbeta = 0.999999_DP / real(NLAP - 1, kind = DP)  


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
      do inner_outer_sign = -1,1,2
         y = (1._DP + inner_outer_sign * 1.0_DP / real(m, kind=DP))**(2._DP / 3._DP) * a   !resonance location for first order resonances
         j = ringmoons_ring_bin_finder(ring, y) !disk location of resonance
         if ((j > 0).and.(j < ring%N + 1)) then 
            select case(inner_outer_sign)
            case(-1) 
               beta = ring%r(j) / a
            case(1)
               beta = a / ring%r(j)
            end select
            !lap =  m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) 
            !dlap = 0.5_DP * beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1)
            !Amk = (lap + dlap)
            !write(*,*) 'full Laplace: ',Amk
 
            lap  =  lapm(inner_outer_sign,m)
            dlap = dlapm(inner_outer_sign,m)

            Amk = (lap + dlap)
            !write(*,*) 'table Laplace: ',Amk
            !read(*,*)
            dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) * &
                      ring%Gsigma(j) * (ring%r(j)**2 * beta * ring%w(j) * Gm / swifter_pl1P%mass * Amk)**2
            width = sqrt(Gm / swifter_pl1P%mass) * ring%r(j)
            w1 = ringmoons_ring_bin_finder(ring,ring%r(j) - width)
            w2 = ringmoons_ring_bin_finder(ring,ring%r(j) + width)
            nw = real(w2 - w1 + 1,kind=DP)
            do w = w1,w2 
               ring%Torque(w) = ring%Torque(w) + dTorque / nw
            end do
            Torque = Torque - dTorque
         end if
      end do
   end do

   return
end function ringmoons_lindblad_torque
