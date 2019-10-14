!**********************************************************************************************************************************
!
!  Unit Name   : helio_kickvb
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Kick barycentric velocities of planets
!
!  Input
!    Arguments : npl        : number of planets
!                helio_pl1P : pointer to head of helio planet structure linked-list
!                dt         : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_kickvb(npl, helio_pl1P, dt)
!
!  Notes       : Adapted from Martin Duncan and Hal Levison's Swift routine kickvh.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_kickvb(npl, helio_pl1P, dt)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_random_access, EXCEPT_THIS_ONE => helio_kickvb
     USE module_interfaces, EXCEPT_THIS_ONE => helio_kickvb
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: npl
     REAL(DP), INTENT(IN)     :: dt
     TYPE(helio_pl), POINTER  :: helio_pl1P

! Internals
     INTEGER(I4B)              :: i
     TYPE(swifter_pl), POINTER :: swifter_plP
     TYPE(helio_pl), POINTER   :: helio_plP

! Executable code
     ! OpenMP parallelization added by D. Minton
     !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE(STATIC) &
     !$OMP SHARED(npl,dt) 
     DO i = 2, npl
          CALL get_point(i,helio_plP)
          swifter_plP => helio_plP%swifter
          swifter_plP%vb(:) = swifter_plP%vb(:) + helio_plP%ah(:) * dt
     END DO
     !$OMP END PARALLEL DO

     RETURN

END SUBROUTINE helio_kickvb
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
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
