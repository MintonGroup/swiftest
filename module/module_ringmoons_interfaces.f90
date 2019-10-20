!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swifter
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in Swifter package
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
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
MODULE module_ringmoons_interfaces

      IMPLICIT NONE


      INTERFACE
         SUBROUTINE ringmoons_io_init_ring(GMP,R_Planet,ring)
            USE module_parameters
            USE module_ringmoons
            IMPLICIT NONE
            real(DP),intent(in)     :: GMP,R_Planet 
            TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         END SUBROUTINE ringmoons_io_init_ring
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_step(lfirst, t, npl, nplmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt,ring)
            USE module_parameters
            USE module_symba
            USE module_ringmoons
            IMPLICIT NONE
            LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
            INTEGER(I4B), INTENT(IN)                         :: npl, nplmax
            REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt
            REAL(DP), INTENT(INOUT)                          :: eoffset
            TYPE(symba_pl), POINTER                          :: symba_pl1P
            TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         END SUBROUTINE ringmoons_step 
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_pde_solver(dtin,ring)
         USE module_parameters
         USE module_ringmoons
         IMPLICIT NONE
         real(DP),intent(in) :: dtin
         TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         END SUBROUTINE ringmoons_pde_solver
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_allocate(ring)
         USE module_parameters
         USE module_ringmoons
         IMPLICIT NONE
         real(DP),intent(in) :: dtin
         TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         END SUBROUTINE ringmoons_allocate
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_deallocate(ring)
         USE module_parameters
         USE module_ringmoons
         IMPLICIT NONE
         real(DP),intent(in) :: dtin
         TYPE(ringmoons_ring),INTENT(INOUT) :: ring
         END SUBROUTINE ringmoons_deallocate
      END INTERFACE





END MODULE module_ringmoons_interfaces
