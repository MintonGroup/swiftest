!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_construct
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : solves
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
!  Invocation  : CALL ringmoons_seed_construct(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_construct(swifter_pl1P,ring)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_construct
      implicit none

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring), intent(inout) :: ring

! Internals
      integer(I4B)                        :: i

! Executable code

      
!#initial list of satellites
!    x = r[i_rigid]      #location at edge of disk
!    while x < r[int(N-1)]:    #run until the distance is larger than the extent of the bins for the disk
!        Sat_M.append(0.0)
!        Sat_r.append(x)         #put a placeholder for a satellite at this location
!        Sat_r_default.append(x)
!        Sat_v.append(sqrt(G*M/x))
!        Sat_e.append(eccentricity())
!        # Sat_time.append(0.0)
!        Sat_w.append(sqrt(G*M/x**3.))
!        # Sat_v.append(0.0)
!        Sat_Rad.append(r_pdisk) #inital satellite size is that of a disk particle
!        rhill.append(x*(m_pdisk/(M*3.))**(1./3.))   #hill radius for disk particle at this location
!        hill = 4.0*x*(m_pdisk/(M*3.))**(1./3.)  #4*hill radius for a disk particle at location x
!        x = x + hill    #advance location to hill (like embryo spacing)
!



      return

end subroutine ringmoons_seed_construct
