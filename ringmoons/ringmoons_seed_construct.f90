!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_construct
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Constructs the satellite seeds. Given the number of seeds, this will find the mass of each seed needed to generate
!                seeds spaced by spacing_factor times the Hill's radius between the Fluid Roche Limit (FRL) and the outer bin of the
!                ring
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
      type(ringmoons_ring), intent(inout) :: ring
      type(ringmoons_seeds), intent(inout) :: seeds

! Internals
      integer(I4B)                        :: i,j,seed_bin
      real(DP)                            :: a, Gm, Gmleft
      logical(LGT)                        :: open_space

! Executable code
  
      
      ! Make seeds small enough to fit into each bin 
      do i = ring%iFrl,ring%N
         if (ring%Gm(i) < INITIAL_MASS_FACTOR * ring%Gm_pdisk) cycle
         open_space = .true.
         do j = 1, seeds%N
            if (.not.seeds%active(j)) cycle
            if ((i >= seeds%fz_bin_inner(j)) .and. (i <= seeds%fz_bin_outer(j))) then
               open_space = .false. ! There is already a seed with a feeding zone here
               exit
            end if
         end do
         if (open_space) then
            a = ring%r(i)
            Gm = INITIAL_MASS_FACTOR * ring%Gm_pdisk !3 * swifter_pl1P%mass * ((ring%router(i) - ring%rinner(i)) / (1.25_DP * FEEDING_ZONE_FACTOR * a))**3
            call ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,Gm)
         end if
      end do             

end subroutine ringmoons_seed_construct
