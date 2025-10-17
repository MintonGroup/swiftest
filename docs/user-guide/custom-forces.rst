#######################
Including custom forces
#######################

.. rubric:: by David A. Minton

Swiftest includes a mechanism that allows users to add their own custom forces to act on bodies in their simulations. This will require writing a Fortran function and compiling a local copy of the Swiftest executable. Be sure to review instructions on how to compile from the :doc:`Getting Started Guide <../getting-started-guide/index>`.

First, create a file called ``swiftest_user.f90`` using the following template:

.. code:: fortran

   submodule(swiftest) s_swiftest_user
      use swiftest
   contains
      module subroutine swiftest_user_kick_getacch_body(self, nbody_system, param, t, lbeg)
         !! author: Your Name
         !!
         !! Description of the function
         !!
         !! Adapted from David E. Kaufmann's Swifter routine whm_user_kick_getacch.f90
         !! 
         implicit none
         ! Arguments
         class(swiftest_body),         intent(inout) :: self   
            !! Swiftest particle data structure
         class(swiftest_nbody_system), intent(inout) :: nbody_system 
            !! Swiftest nbody_system_object
         class(swiftest_parameters),   intent(inout) :: param  
            !! Current run configuration parameters user parameters
         real(DP),                     intent(in)    :: t      
            !! Current time
         logical,                      intent(in)    :: lbeg   
            !! Logical flag that determines whether or not this is the beginning or end of the step

         ! Internals: Add your own variables here
         integer(I4B) :: i

         ! Your code goes here
         ! do i = 1, self%nbody
         ! ...
         ! end do

         return
      end subroutine swiftest_user_kick_getacch_body

   end submodule s_swiftest_user

By default, CMake will look for a ``swiftest_user.f90`` file here and automatically include it into the build if it exists. However, if you wish to you can place this file in a different directory and pass ``-DSWIFTEST_USER_DIR=/full/path/to/folder`` to CMake (see build instructions `here <https://swiftest.readthedocs.io/en/latest/getting-started-guide/index.html#building-the-executable-using-cmake>`__ for details of how to set CMake options).

The subroutine ``swiftest_user_kick_getacch_body`` is called during each kick step of the integrator. Do not alter the name or arguments of this subroutine.  The arguments are as follows:

- ``self``: Swiftest particle data structure. This contains all of a given type of body in the simulation. That is, this subroutine will be called on all massive bodies (``pl``) and test particles (``tp``) during the evaluation of the "kick" step. This extra force is the last set of accelerations that are computed, after accelerations due to ``pl``-body interactions and, if enabled, non-spherical central body forces and general relativity. Therefore, the user-supplied accelerations should be added to the existing acceleration values in ``self%ah``, which is an (3,n) array of accelerations. 

.. note:: Fortran stores arrays in column-major order with 1-indexing, in constrast to C and Numpy which use row-major order with 0-indexing. Therefore, the x, y, and z components of the acceleration of body ``i`` are stored in ``self%ah(1,i)``, ``self%ah(2,i)``, and ``self%ah(3,i)``, respectively. 