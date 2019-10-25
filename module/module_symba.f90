!**********************************************************************************************************************************
!
!  Unit Name   : module_symba
!  Unit Type   : module
!  Project     : SWIFTEST
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Symplectic Massive Body Algorithm
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
MODULE module_symba

     USE module_parameters
     USE module_helio
     IMPLICIT NONE

     INTEGER(I4B), PARAMETER :: NENMAX = 32767
     INTEGER(I4B), PARAMETER :: NTENC = 3
     REAL(DP), PARAMETER     :: RHSCALE = 6.5_DP
     REAL(DP), PARAMETER     :: RSHELL = 0.48075_DP

     ! Added by D. Minton
     !TYPE symba_ptr_arr
     !   TYPE(symba_pl), POINTER :: thisP   ! pointer to current SyMBA planet
     !END TYPE symba_ptr_arr

     !TYPE symba_ptr_arr_tp
     !   TYPE(symba_tp), POINTER :: thisP   ! pointer to current SyMBA planet
     !END TYPE symba_ptr_arr_tp
     !^^^^^^^^^^^^^^^^^^^
     type symba_pl
          logical(LGT), dimension(:),     allocatable :: lmerged ! flag indicating whether body has merged with another this time step
          integer(I4B), dimension(:),     allocatable :: nplenc  ! number of encounters with other planets this time step
          integer(I4B), dimension(:),     allocatable :: ntpenc  ! number of encounters with test particles this time step
          integer(I4B), dimension(:),     allocatable :: levelg  ! level at which this body should be moved
          integer(I4B), dimension(:),     allocatable :: levelm  ! deepest encounter level achieved this time step
          integer(I4B), dimension(:),     allocatable :: nchild  ! number of children in merger list
          integer(I4B), dimension(:),     allocatable :: isperi  ! perihelion passage flag
          real(DP),     dimension(:),     allocatable :: peri    ! perihelion distance
          real(DP),     dimension(:),     allocatable :: atp     ! semimajor axis following perihelion passage
          type(helio_pl)                     :: helio   ! HELIO planet structure
     end type symba_pl

     type symba_tp
          integer(I4B), dimension(:),     allocatable :: nplenc  ! number of encounters with planets this time step
          integer(I4B), dimension(:),     allocatable :: levelg  ! level at which this particle should be moved
          integer(I4B), dimension(:),     allocatable :: levelm  ! deepest encounter level achieved this time step
          type(helio_tp)                     :: helio   ! HELIO test particle structure
     end type symba_tp

     type symba_plplenc
          logical(LGT), dimension(:),     allocatable :: lvdotr ! relative vdotr flag
          integer(I4B), dimension(:),     allocatable :: status ! status of the interaction
          integer(I4B), dimension(:),     allocatable :: level  ! encounter recursion level
          !TODO: Pointer or arrays?
          integer(I4B), dimension(:),     allocatable :: id1     ! external identifier first planet in encounter
          integer(I4B), dimension(:),     allocatable :: id2     ! external identifier second planet in encounter
          !type(symba_pl), POINTER :: pl1P   ! pointer to first planet in encounter
          !type(symba_pl), POINTER :: pl2P   ! pointer to second planet in encounter
     end type symba_plplenc

     type symba_pltpenc
          logical(LGT), dimension(:),     allocatable :: lvdotr ! relative vdotr flag
          integer(I4B), dimension(:),     allocatable :: status ! status of the interaction
          integer(I4B), dimension(:),     allocatable :: level  ! encounter recursion level
          !TODO: Pointer or arrays?
          integer(I4B), dimension(:),     allocatable :: idpl    ! external identifier planet in encounter
          integer(I4B), dimension(:),     allocatable :: idtp    ! external identifier test particle in encounter

          !type(symba_pl), POINTER :: plP    ! pointer to planet in encounter
          !type(symba_tp), POINTER :: tpP    ! pointer to test particle in encounter
     end type symba_pltpenc

     type symba_merger
          integer(I4B), dimension(:),     allocatable :: id     ! external identifier
          integer(I4B), dimension(:),     allocatable :: status ! status
          integer(I4B), dimension(:),     allocatable :: ncomp  ! number of component bodies in this one during this merger
          real(DP),     dimension(:,:),   allocatable :: xh     ! heliocentric position
          real(DP),     dimension(:,:),   allocatable :: vh     ! heliocentric velocity
     end type symba_merger 
END MODULE module_symba
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
