##########################################
User guide for using the Yarkovsky effect
##########################################

.. rubric:: by Kaustub Anand

The Yarkovsky effect is a well known radiation-based force that alters the orbits of 10cm - 10km sized bodies, primarily asteroids, as they orbit around the Sun. 
We can now model the Yarkovsky effect as an added force in Swiftest simulations. This is only applicable to massive particles and our implementation is based on that of REBOUNDx, described in  
`Ferich et al. (2022) <https://doi.org/10.3847/1538-4365/ac8d60>`__. 

Set the Yarkovsky parameter flag
====================================
Let's start by setting up a simulation object and turning on the Yarkovsky effect in ``swiftest.Simulation()``. We will set the unit system, time parameters, and then add the Sun as the central body.

.. code-block:: python

    import swiftest

    sim = swiftest.Simulation(tstop = 1e6, dt = 0.01, tstep_out = 1e3, yarkovsky = True, DU = 'AU', TU = 'yr', MU = 'kg')
    sim.add_solar_system_body(name = ["Sun"])

Necessary characteristics
============================================================
Massive particles need to be added with a rotation vector and 4 additional characteristics. They are as follows with the Swiftest parameter name in brackets:

- Emissivity (``emissivity``): The emissivity of the body.
- Albedo (``albedo``): The bond albedo of the body.
- Thermal Inertia (``gamma``): The thermal inertia of the body where :math:`\Gamma = \sqrt{k \rho C}`.
    - :math:`\Gamma` is the thermal inertia
    - :math:`k` is the thermal conductivity
    - :math:`\rho` is the density
    - :math:`C` is the specific heat capacity
- Rotational Constant (``rot_k``): A rotational constant between 0 and 0.25 that depends on its rotation rate.
    - If the rotation rate approaches the critical break-up spin, :math:`rot_k -> 0`
    - If the rotation period approaches the orbital period, :math:`rot_k -> 0.25`
    - A potential equation to calculate :math:`rot_k = [1 - (\omega_{critical}/\omega)^2] * 0.25`

Except for ``gamma``, all the quantities above are unitless. The thermal inertia has SI units of :math:`Jm^{-2}K^{-1}s^{-1/2}`. 
We do not have a unit conversion for Kelvin (:math:`K`) in Swiftest. To convert the thermal inertia from SI to `sim` units we can follow
:math:`\Gamma_{sim} = \Gamma_{SI} *` ``(sim.TU2S)``:math:`^{2.5} /` ``sim.MU2KG``

Add in massive bodies 
================================================
Now that we know our required characteristics and the units are correct, we can add bodies to the simulation. For our example, we will add asteroid 11470 Davidminton and assume typical values for some it's physical properties because of lack of data. 
We will assume a rotation period of 6 hours, an albedo of 0.07, emissivity of 0.9, rot_k of 0.25, density of 1300 :math:`kg/m^3`, and a thermal inertia of 100 :math:`Jm^{-2}K^{-1}s^{-1/2}`.

.. code-block:: python

    sim.add_body(name = '11470 Davidminton', radius = 2.007e3 / sim.DU2M,       # m to AU
                            mass = 4.4022e13,                                   # kg
                            rot = np.array([0.0, 0.0, 5.2596e5]),               # deg/yr
                            a = 2.556363849,                                    # AU
                            e = 0.2087687, 
                            inc = 6.007875886,                                  # deg
                            capom = 0.0, omega = 0.0, capm = 0.0,               # deg
                            albedo = 0.07, 
                            emissivity = 0.9, 
                            gamma = 100.0 * (sim.TU2S**(5.0/2)) / sim.MU2KG,    # sim units 
                            rot_k = rot_k)

Run the simulation
==========================
Set any other relevant parameters or add in more bodies. We can now run the simulation.

.. code-block:: python

    sim.run()