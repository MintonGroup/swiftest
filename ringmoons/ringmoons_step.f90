!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_step
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step the fluid ring in time
!
!  Input
!    Arguments   swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!                ring           : Ringmoons data structure
!                dtin           : input time step
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                ring           : Updated ring
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_step(swifter_pl1P,ring,dtin,lfirst)
!
!  Notes       : Adapted from Andy Hesselbrock's RING-MOONS Python scripts
!
!**********************************************************************************************************************************
subroutine ringmoons_step(swifter_pl1P,ring,seeds,dtin,lfirst,Merror,Lerror)

! Modules
     use module_parameters
     use module_swifter
     use module_ringmoons
     use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_step
     implicit none

! Arguments
      type(swifter_pl), pointer                       :: swifter_pl1P
      real(DP), intent(in)                            :: dtin
      type(ringmoons_ring),intent(inout)              :: ring
      type(ringmoons_seeds),intent(inout)             :: seeds
      logical(LGT), intent(inout)                     :: lfirst
      real(DP),intent(out)                            :: Merror,Lerror

! Internals
      integer(I4B) :: i,loop
      real(DP) :: dtstab,dtleft,dt,seedmass
      real(DP),save :: Mtot_orig,Mtot_now,Ltot_orig,Ltot_now

! Executable code
      !if (lfirst) then
      dtleft = dtin
      !TESTING
         if (lfirst) then
            Mtot_orig = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm)
            Ltot_orig = sum(seeds%Gm(:) * sqrt(swifter_pl1P%mass * seeds%a(:))) 
            Ltot_orig = Ltot_orig + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
            Ltot_orig = Ltot_orig + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
            lfirst = .false.
         end if
         call ringmoons_viscosity(ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         Mtot_now = swifter_pl1P%mass + sum(ring%Gm) + sum(seeds%Gm)
         Ltot_now = sum(seeds%Gm(:) * sqrt(swifter_pl1P%mass * seeds%a(:))) 
         Ltot_now = Ltot_now + sum(ring%Gm(:) * ring%Iz(:) * ring%w(:))
         Ltot_now = Ltot_now + swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
         Merror =  (Mtot_now - Mtot_orig) / Mtot_orig
         Lerror = (Ltot_now - Ltot_orig) / Ltot_orig
      !^^^^^^^^  
      do loop = 1, LOOPMAX
         call ringmoons_viscosity(ring)
         dtstab = ring%stability_factor / maxval(ring%nu)
         dt = min(dtleft,dtstab)
         call ringmoons_sigma_solver(ring,dt)
         ring%Gm = ring%Gsigma * ring%deltaA
         call ringmoons_seed_grow(swifter_pl1P,ring,seeds,dt)
         call ringmoons_planet_accrete(swifter_pl1P,ring,seeds)
         dtleft = dtleft - dt
         if (dtleft <= 0.0_DP) exit
      end do 


      RETURN

END SUBROUTINE ringmoons_step
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
