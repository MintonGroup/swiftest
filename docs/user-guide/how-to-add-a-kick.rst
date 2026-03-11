##########################################################
A guide to adding a new kick to Swiftest
##########################################################

.. rubric:: by Kaustub Anand

Here, we provide a guide on how to add a kick or additional feature, including new variables and particle features, to Swiftest.
This guide is written in the context of adding a new kick or acceleration feature that requires new properties for the particles.
However, the general steps can be extended to adding other features as well. The guide is split up into two sections, the Fortran 
and Python sides of Swiftest. The former handles the calculations and data, while the latter helps with user interface.

Important nomenclature:
=======================

- Nbody system = ``nbody_system`` = nbody system object
- Central Body = ``cb`` = ``nbody_system%cb``
- Particles are of type ``body`` with subclasses:
    - Massive Particles = ``pl`` = type ``body`` of class ``swiftest_pl``
    - Test particles = ``tp`` = type ``body`` of class ``swiftest_tp``
- Input parameters = ``param`` = simulation paramters object

- *file.f90* is a file or folder.
- ``code`` is a block of code.

Note: small details such as ensuring that internal variables are defined for a given code block and correct syntax are not mentioned and left to the user.

Fortran
========
The Fortran side of Swiftest is where all the calculations for the simulation are performed. Going into *swiftest/src/* To add a kick:

- Check and add flags for the kick to the parameter object ``param`` in *base_module.f90*:

    - define and add a flag for the kick as an attribute for the ``base_parameters`` class in *base_module.f90*.
        - For example, a flag for the Yarkovsky effect would be ``lyarkovsky``

        .. code-block:: Fortran

            type, abstract :: base_parameters
                .
                .
                logical :: lyarkovsky = .false. 
                    !! Turn on Yarkovsky effect
                .
                .   

- Write the kick in a new file *swiftest_\*.f90*:

    - Create a new fortran file in *swiftest/src/swiftest* following the naming convention of *swiftest_\*.f90* where * is the name of the kick. 
        - You can also add a subroutine for a kick inside another related file.
        - For example, a Yarkovsky related kick could have a new file called *swiftest_yarkovsky.f90* or be added to a file with other radiation kicks such as *swiftest_radiation.f90*.
    - The file has to be set up as a ``submodule`` that contains various ``module subroutines`` for individual kicks. Below is an example of a yarkovsky kick for massive particles.
        
        .. code-block:: Fortran

            submodule (swiftest) s_swiftest_yarkovsky

            contains

                module subroutine swiftest_getacch_yarkovsky_pl(self, nbody_system, param)
                    !!
                    !! Calculate the Yarkovsky effect on massive bodies. 
                    !!
                    implicit none
                    ! Arguments
                    class(swiftest_pl),           intent(inout) :: self
                        !! Swiftest body object
                    class(swiftest_nbody_system), intent(inout) :: nbody_system
                        !! Swiftest nbody system object
                    class(swiftest_parameters),   intent(in)    :: param
                        !! Current run configuration parameters
                    .
                    .
                end subroutine swiftest_getacch_yarkovsky_pl
                    .
                    .
                    .
            end submodule s_swiftest_yarkovsky
    - ``subroutines`` follow the naming convention of ``swiftest_getacch_NAME_PARTICLETYPE``.
    - Write the formulas for calculating the kick in the appropriate ``subroutine`` add the kick to the ``i``-th particle's acceleration vector.
        .. code-block:: Fortran
            
            .
            .
            module subroutine swiftest_getacch_yarkovsky_pl(self, nbody_system, param)
                implicit none
                .
                .
                a_kick(:) = a_kick_mag * a_kick_direction(:)
                pl%ah(:, i) = pl%ah(:, i) + a_kick(:)
                .
                .
            end subroutine swiftest_getacch_yarkovsky_pl
            .
            .

- Add the kick and necessary particle parmeters to *swiftest_module.f90*:
    
    - Define and add any new particle parameters required by the kick to the appropriate particle type(s) (``swiftest_pl``, ``swiftest_tp``, or ``swiftest_body`` for both particles).
        - For example, we can add emissivity and albedo to massive particles (``pl``) for the Yarkovsky kicks.

        .. code-block:: Fortran

            type, abstract, extends(swiftest_body) :: swiftest_pl
                .
                .
                .
                real(DP),                dimension(:), allocatable   :: albedo
                    !! Bond albedo for radiation acceleration calculations
                real(DP),                dimension(:), allocatable   :: emissivity
                    !! Emissivity for Yarkovsky acceleration calculations
                .
                .
                .
            contains
                .
                .

    - Add the kick with a simple name as a procedure attribute for the appropriate particle type(s). 
        - For example, we can add the Yarkovsky effect to ``swiftest_pl`` only.

        .. code-block:: Fortran

            type, abstract, extends(swiftest_body) :: swiftest_pl
                .
                .
                .
            contains
                .
                .
            procedure :: accel_yarkovsky => swiftest_yarkovsky_getacch_pl 
                !! Compute the heliocentric accelerations of bodies due to the Yarkovsky effect
    
    - Add the kick submodule definition in an appropriate interface.

        .. code-block:: fortran
            
                .
                .
                .
            interface
                .
                .
                .
                module subroutine swiftest_yarkovsky_getacch_pl(self, nbody_system, param)
                    !!
                    !! Calculate the Yarkovsky effect on massive bodies. 
                    !!
                    implicit none
                    ! Arguments
                    class(swiftest_pl),           intent(inout) :: self
                        !! Swiftest body object
                    class(swiftest_nbody_system), intent(inout) :: nbody_system
                        !! Swiftest nbody system object
                    class(swiftest_parameters),   intent(in)    :: param
                        !! Current run configuration parameters
                end subroutine swiftest_yarkovsky_getacch_pl
                .
                .
                .
            end interface
                .
                .
                .

- Add the kicks to the code flow in *whm_kick.f90* and *helio_kick.f90*:

    - All the current integrators (RMVS, WHM, SyMBA, HELIO) run though *whm_kick.f90* and *helio_kick.f90*. 
    - Add the appropriate flag check and kick function per particle type in ``whm_kick_getacch_pl(...)``, ``whm_kick_getacch_tp(...)``, ``helio_kick_getacch_pl(...)``, and/or ``helio_kick_getacch_tp(...)``.
        
        .. code-block:: fortran

            module subroutine helio_kick_getacch_pl(self, nbody_system, param, t, lbeg)
                .
                .
                if (param%lyarkovsky) call pl%accel_yarkovsky(nbody_system, param)
                .
                .
                return
            end subroutine helio_kick_getacch_pl

- In the home directory of Swiftest, add the new files to *src/CMakeLists.txt* under the right section for fortran compilation:

    - A kick will typically go under ``STRICT_MATH_FILES``

    .. code-block:: Fortran

        # Add the source files
        SET(STRICT_MATH_FILES
            .
            .
            ${SRC}/swiftest/swiftest_yarkovsky.f90
            .
            .
        ) 

Now we can tackle the data and parameter I/O. If you have not added any new particle parameters/features (In our example, ``albedo`` and ``emissivity``) or new variables to the data, you can ignore this section.
Swiftest uses hdf5 and NetCDF files for data handling. We will split the I/O section into general Swiftest data checks and then NetCDF I/O checks.

- Initial I/O checks

- NetCDF I/O checks

Python
=======
