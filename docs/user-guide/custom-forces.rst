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

- ``self``: Swiftest particle data structure. This contains all of a given type of body in the simulation. That is, this subroutine will be called on all massive bodies (``pl``) and test particles (``tp``) during the evaluation of the "kick" step. This extra force is the last set of accelerations that are computed, after accelerations due to ``pl``-body interactions and, if enabled, non-spherical central body forces and general relativity. Therefore, the user-supplied accelerations should be added to the existing acceleration values in ``self%ah``, which is an (3,n) array of accelerations. For a complete list of all variables and methods that are a part of the ``swiftest_body`` class, see the API reference for `swiftest_body <../_static/fortran_docs/type/swiftest_body.html>`_.

- ``nbody_system``: Swiftest nbody_system object. This contains all of the bodies in the simulation (including the one represented by the ``self`` argument) as well as the central body (``cb``).  For a complete list of all variables and methods that are part of the ``swiftest_nbody_system`` class, see the API reference for `swiftest_nbody_system <../_static/fortran_docs/type/swiftest_nbody_system.html>`_.

- ``param``: Current run configuration parameters user parameters. This contains all of the simulation parameters as set by the user in the input file. For a complete list of all variables and methods athat are part of the ``swiftest_parameters`` class, see the API reference for `swiftest_parameters <../_static/fortran_docs/type/swiftest_parameters.html>`_.

- ``t``: Current time in the simulation.

- ``lbeg``: Logical flag that determines whether or not this is the beginning or end of the step. This can be used if your custom force depends on whether it is being evaluated at the beginning or end of the kick step.

.. note:: Fortran stores arrays in column-major order with 1-indexing, in constrast to C and Numpy which use row-major order with 0-indexing. Therefore, the x, y, and z components of the acceleration of body ``i`` are stored in ``self%ah(1,i)``, ``self%ah(2,i)``, and ``self%ah(3,i)``, respectively. All arrays representing cartesian vectors in Swiftest follow this convention.

Basic example
=============
As a simple example, we will create a custom force that applies a constant acceleration in the +z direction to all bodies in the simulation. This is not meant to be a physically realistic force, but it will illustrate how to implement a custom force. The code snippet below shows only what goes in the line marked **Your code goes here** in the template above.

.. code:: fortran

         ! Apply a constant acceleration in the +z direction of 1.0e-10 DU/TU^2 
         do i = 1, self%nbody
            self%ah(3,i) = self%ah(3,i) + 1.0e-10_DP  
         end do


When the code is compiled in Swiftest and extra forces are enabled, either by passing ``extra_force = True`` to the ``Simulation.run()`` method in Python or by setting ``EXTRA_FORCE = YES`` in the Fortran input file, this force will be applied to all bodies during each kick step of the integrator.

Using derived types
===================

Because Swiftest is object-oriented, the arguments ``self``, ``nbody_system``, and ``param`` are *class* variables, which means that they can be extended with new variables and methods while still being passed. Your custom code must be written with this in mind, as this function will be called for the collection of massive bodies ``pl`` and ``tp`` separately. In addition, these derived types are defined for each of the integrators, so for SyMBA simulations, the massive bodies and test particles will be of type ``symba_pl`` and ``symba_tp``, respectively. Modern Fortran includes a construct specifically designed to access the variables and methods of a derived type using the ``select type`` and ``class is`` statements. To illustrate how to use this, we will show a simple example in which the user-defined force acts only on test particles in a simulation.   

.. code:: fortran

         select type(body => self)
         class is (swiftest_tp)
            ! Apply a constant acceleration in the +z direction of 1.0e-10 DU/TU^2 to test particles only
            do i = 1, body%nbody
               body%ah(3,i) = body%ah(3,i) + 1.0e-10_DP  
            end do
         class is (swiftest_pl)
            ! Do nothing for massive bodies
         end select

We can also uses this same approach to apply a mass-dependent force to massive bodies. The variables ``mass`` (as well as the ``Gmass``, which is the product of the gravitational constant and the mass) are defined only for massive bodies, so we must use the ``select type`` construct to access them, as they are not defined in the ``swiftest_body`` class. The following code snippet shows how to apply a force proportional to the mass of each body in the simulation.

.. code:: fortran

         select type(body => self)
         class is (swiftest_pl)
            ! Apply a mass-dependent acceleration in the +z direction to massive bodies only
            do i = 1, body%nbody
               body%ah(3,i) = body%ah(3,i) + 1.0e-10_DP * body%Gmass(i)  
            end do
         class is (swiftest_tp)
            ! Do nothing for test particles
         end select

-------------------------------
Using the nbody_system argument
-------------------------------

You can also access all particles in the system using the ``nbody_system`` argument. For example, the following code snippet shows how to apply a force to all test particles that depends on the total mass of all massive bodies in the system.

.. code:: fortran

         real(DP) :: total_mass
         integer(I4B) :: j

         ! First, compute the total mass of all massive bodies in the system
         total_mass = 0.0_DP
         do j = 1, nbody_system%pl%nbody
            total_mass = total_mass + nbody_system%pl%mass(j)
         end do

         ! Now apply a force to all test particles that depends on the total mass
         select type(body => self)
         class is (swiftest_tp)
            do i = 1, body%nbody
               body%ah(3,i) = body%ah(3,i) + 1.0e-10_DP * total_mass  
            end do
         class is (swiftest_pl)
            ! Do nothing for massive bodies
         end select

-----------------------------
Selecting a specific particle
-----------------------------

Suppose you want to apply a custom force to only one body in the simulation. You could select it by index, but this is not very robust, as the index of a body may change if bodies are added or removed from the simulation. A better approach is to select a body by either its integer id or string name. Suppose we want to apply a custom force only on planet Jupiter. We generate the initial conditions of the simulation using the Python interface:

.. ipython:: python
    :okwarning:
    :suppress:

    import os
    import tempfile
    tmpdir=tempfile.TemporaryDirectory()
    os.environ['OMP_NUM_THREADS'] = '1'
    sim = swiftest.Simulation(simdir=tmpdir.name)


.. ipython:: python
    :okwarning:
    
    import swiftest
    sim = swiftest.Simulation()
    sim.add_solar_system_body(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"])
    print(f"Jupiter's id is {sim.init_cond.sel(name='Jupiter').id.values[()]}")


Therefore we could apply a force only to Jupiter using the following code snippet:


.. code:: fortran

         integer(I4B) :: jupiter_id
         integer(I4B) :: i

         ! Assume Jupiter has an id of 1 (as it is the first body added after the Sun)
         jupiter_id = 1

         select type(body => self)
         class is (swiftest_pl)
            do i = 1, body%nbody
               if (body%id(i) == jupiter_id) then
                  ! Apply a custom acceleration to Jupiter only
                  body%ah(3,i) = body%ah(3,i) + 1.0e-10_DP  
               end if
            end do
         end select

Now suppose we create the initial conditions differently, such that Jupiter is not the first body added after the Sun:

.. ipython:: python
    :okwarning:

    import swiftest
    sim = swiftest.Simulation()
    sim.add_solar_system_body(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"])
    print(f"Jupiter's id is {sim.init_cond.sel(name='Jupiter').id.values[()]}")

In this case, Jupiter will be id ``5`` instead of ``1``. This change would require recompiling the Fortran code to change the value of ``jupiter_id``. Alternatively, we could make this more robust by searching for Jupiter's id using its name, however matching strings in Fortran is somewhat cumbersome and computationally expensive, and the acceleration evaluation is called many times during a simulation, so it is critical not to introduce unnecessary overhead. A better approach is to search for Jupiter's id only once, during the first call to the subroutine, and store it in a static variable for later use. This can be done in Fortran by specifing the ``save`` attribute for a variable. This will allow us to store Jupiter's id after we find it the first time, and then use it in subsequent calls without having to search for it again. The following code snippet shows how to do this:


.. code:: fortran

         integer(I4B), save :: jupiter_id = -1
         integer(I4B) :: i
         integer(I4B) :: j

         ! If jupiter_id is -1, this means we haven't found it yet
         if (jupiter_id == -1) then
            ! Search through all massive bodies to find Jupiter's id
            do j = 1, nbody_system%pl%nbody
               if (trim(adjustl(nbody_system%pl%info(j)%name)) == "Jupiter") then ! This will match the substring, ignoring leading/trailing spaces
                  jupiter_id = nbody_system%pl%id(j)
               end if
            end do
         end if

         select type(body => self)
         class is (swiftest_pl)
            do i = 1, body%nbody
               if (body%id(i) == jupiter_id) then
                  ! Apply a custom acceleration to Jupiter only
                  body%ah(3,i) = body%ah(3,i) + 1.0e-10_DP  
               end if
            end do
         end select

This version of the code will work regardless of the order in which bodies are added to the simulation.

