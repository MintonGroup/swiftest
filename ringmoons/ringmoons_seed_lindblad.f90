!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_lindblad
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the lindblad torques
!                for distance from the FRL
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_seed_lindblad(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_lindblad(ring,seeds)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_lindblad
      implicit none

! Arguments
      type(ringmoons_ring), intent(inout)    :: ring
      type(ringmoons_seeds), intent(inout)      :: seeds

! Internals
      integer(I4B)                           :: i,j,j0,mode
      real(DP)                               :: y, Torque, del_M, y0
      real(DP),dimension(seeds%N)            :: Gamma_S

! Executable code

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
end subroutine ringmoons_seed_lindblad
