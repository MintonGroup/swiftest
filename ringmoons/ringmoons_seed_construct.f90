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
subroutine ringmoons_seed_construct(swifter_pl1P,ring,seeds)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_construct
      implicit none

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring), intent(in) :: ring
      type(ringmoons_sat), intent(inout) :: seeds

! Internals
      integer(I4B)                        :: N,iFRL,iend
      integer(I4B),parameter              :: MAXSEEDS = 100000 ! Maximum possible number of seeds
      real(DP),dimension(MAXSEEDS)        :: a 
      real(DP)                            :: rhill
      real(DP),parameter                  :: spacing_factor = 4.0_DP

! Executable code
      
!#initial list of satellites
      iFRL = ring%iFRL 
      iend = ring%N
      N = 1
      a(N) = ring%r(iFRL)
      do
         rhill = a(N) * (ring%Gm_pdisk / (3 * swifter_pl1P%mass))**(1._DP / 3._DP)
         a(N + 1) = a(i) + spacing_factor * rhill
         if (a(N + 1)  > ring%router(iend)) exit
         N = N + 1
      end do
      allocate(seeds%a(N))
      allocate(seeds%m(N))
      seeds%N = N
      seeds%a(:) = a(1:N)
      seeds%m(:) = ring%Gm_pdisk
      write(*,*) N,' seeds created'
      read(*,*)
       
      return

end subroutine ringmoons_seed_construct
