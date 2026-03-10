##########################################################
A guide to adding a new kick to Swiftest
##########################################################

.. rubric:: by Kaustub Anand

Here, we provide a guide on how to add a kick or additional feature, including new variables and particle features, to Swiftest.
This guide is written in the context of adding a new kick or acceleration feature that requires new properties for the particles.
However, the general steps can be extended to adding other features as well. The guide is split up into two sections, the Fortran 
and Python sides of Swiftest. The former handles the calculations and data, while the latter helps with user interface.

Nomenclature:
#############
- Massive Particles = ``pl`` = type ``body`` of class ``swiftest_pl``
- Test particles = ``tp`` = type ``body`` of class ``swiftest_tp``
- Central Body = ``cb`` = ``nbody_system%cb``
- Parameter input file = ``param``

- *file.f90* is a file or folder.
- ``code`` is a block of code.

Note: small details such as ensuring that internal variables are defined for a given code block and correct syntax are not mentioned and left to the user.

Fortran
========
The Fortran side of Swiftest is where all the calculations for the simulation are performed. To add a kick:

- Check and add flags to the parameter file/object ``param`` in *base_module.f90*:

    - define and add a flag for the kick as an attribute for the ``base_parameters`` class in *base_module.f90*.
        - For example, a flag for the Yarkovsky effect would be ``lyarkovsky``
        - ``logical :: lyarkovsky = .false. !! Turn on Yarkovsky effect``

- Write the kick in a new file *swiftest_\*.f90*:

    - Create a new fortran file in *swiftest/src/swiftest* following the naming convention of *swiftest_\*.f90* where * is the name of the kick. 
        - For example, a radiation related kick could be called *swiftest_radiation.f90*. 
    - make submodule .......
    - Write the formulas for calculating the kick in *swiftest_\*.f90*. 
        - ``a_kick(:) = a_kick_mag * a_kick_direction(:)``
    - Add the kick to the ``i``-th particle's acceleration vector .
        - ``pl%ah(:, i) = pl%ah(:, i) + a_kick(:)``

- Add the kick to *swiftest_module.f90*:

    - Add the kick with a simple name as a procedure attribute for the appropriate particle type(s) (``swiftest_pl``, ``swiftest_tp``, or ``swiftest_body`` for both particles). 
        - For example, we can add the Yarkovsky effect to ``swiftest_pl`` only

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

            interface
                .
                .
                .
                module subroutine swiftest_yarkovsky_getacch_pl(self, nbody_system, param)
                    !! author: Kaustub P. Anand and David A. Minton
                    !!
                    !! Calculate the Yarkovsky effect on massive bodies. 
                    !! Based on Ferich, et al, 2022 (https://iopscience.iop.org/article/10.3847/1538-4365/ac8d60) and Veras, et al, 2015 (https://academic.oup.com/mnras/article/451/3/2814/1180328)
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



Python
=======