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
      integer(I4B)                           :: j, m, j0, inner_outer_sign
      integer(I4B), parameter                :: m_max = 10 ! Maximum mode number 
      real(DP)                               :: y, dTorque, beta, Amk

! Executable code

      ! Just do the first order resonances for now. The full suite of resonances will come later
     
      j0 = -1 
      Torque = 0.0_DP
      do m  = 2, m_max
         y = (1._DP + 1._DP / real(m, kind=DP))**(2._DP / 3._DP) * a   !resonance location for first order resonances
         j = ringmoons_ring_bin_finder(ring, y) !disk location of resonance
         if (j == j0) exit ! Resonances overlap in bin space, so our assumptions are no longer valid
         j0 = j
         if ((j == 0).or.(j == ring%N + 1)) cycle
         if (a > ring%r(j)) then ! Inner Lindblad
            beta = ring%r(j) / a
            inner_outer_sign = -1
         else
            beta = a / ring%r(j)
            inner_outer_sign = 1
         end if
         Amk = 0.5_DP * (2 * m * ringmoons_laplace_coefficient(beta,m,0.5_DP,0) + &
                        beta * ringmoons_laplace_coefficient(beta,m,0.5_DP,1))
         dTorque = inner_outer_sign * 4 * PI**2 / (3._DP) * m / real(m - 1, kind=DP) * &
                   ring%Gsigma(j) * (ring%r(j)**2 * beta * ring%w(j) * Gm / swifter_pl1P%mass * Amk)**2
         ring%Torque(j) = ring%Torque(j) + dTorque
         Torque = Torque - dTorque
      end do

!      do i = 1, seeds%N
!         ! first mode is done alone for outer part of disk (inner part of disk would be at the location of Planet
!         mode = 1
!         j0 = -1
!         do
!            y = (1._DP + 1._DP / real(mode, kind=DP)**(2._DP / 3._DP) * seeds%a(i)   !resonance locations
!            j = ringmoons_ring_bin_finder(ring, y)   !disk location of resonance
!            if (j0 == j) exit  !while the distance between Lindblad Resonances is greater than bin width
!            if (mode == 1) then 
!               if (j < ring%N - 1) then   !if within bin boundaries..
!                  Torque = (mode * seeds%Gm(i))**2 * ring%Gsigma(j) / (ring%w(j)**2 * ring%r(j)**2)
!                  Gamma_S(i) = -Torque  !Torque exerted on the satellite is equal and opposite to the torque the satellite exerts on the disk
!               else
!                  Torque = 0.0_DP
!               end if
!            end if
!            j0 = j
!            mode = mode + 1 
!         end do 
!         while 
!
!             if 0 < i < int(N):
!
!                 Torque = -1.*(mode*G*Sat_M[a])**2.*sigma[i]/(w[i]**2.*r[i]**2.) #from Planeary Rings by Esposito
!                 del_m = deltat*Torque/(R[i]*w[i])       #change in mass to LR
!                 if abs(del_m) > m[i]:    #am I trying to move more mass than is available?
!                     del_m = -1.*m[i]    #mass change is negative, so keep it that way
!                     Torque = R[i]*w[i]*del_m/deltat  #scale the torque to equal that amount moved from LR bin
!                 m[i] = m[i] + del_m
!                 sigma[i] = m[i]/deltaA[i]
!                 m[i-1] = m[i-1] - del_m
!                 sigma[i-1] = m[i-1]/deltaA[i-1]
!
!                 Gamma_S.append(-1.*Torque)  #Torque exerted on the satellite is equal and opposite to the torque the satellite exerts on the disk
!
!             y = (1. + 1./mode)**(2./3.)*Sat_r[a]   #resonance locations
!             j = int(round((y - r_I)/deltar - .5))
!         
!             if 0 < j < int(N-1):
!                 Torque = (mode*G*Sat_M[a])**2.*sigma[j]/(w[j]**2.*r[j]**2.) #from Planeary Rings by Esposito
!                 del_m = deltat*Torque/(R[j]*w[j])       #change in mass to LR
!                 if del_m > m[j-1]:    #am I trying to move more mass than is available in the bin it comes from?
!                     del_m = m[j-1]
!                     Torque = R[j]*w[j]*del_m/deltat  #scale the torque to equal that amount moved from LR bin
!
!                 m[j] = m[j] + del_m
!                 sigma[j] = m[j]/deltaA[j]
!                 m[j+1] = m[j+1] - del_m
!                 sigma[j+1] = m[j+1]/deltaA[j+1]
!
!                 Gamma_S.append(-1.*Torque)  #Torque exerted on the satellite is equal and opposite to the torque the satellite exerts on the disk
!         
!             mode = mode + 1
!
!         # ****
!         Torque = 1.*sum(Gamma_S)

!   
!    
      return
end function ringmoons_lindblad_torque
