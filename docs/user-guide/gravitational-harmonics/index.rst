##########################
Gravitational Harmonics
##########################

Here, we show how to use Swiftest's Gravitational Harmonics capabilities. 
Swiftest uses `SHTOOLS <https://shtools.github.io/SHTOOLS/>` to compute the gravitational harmonics coefficients for a given body and calculate it's associated acceleration kick. 
The coefficients can be computed in a number of ways: 
- Using the axes measurements of the body ( :func: `clm_from_ellipsoid <swiftest.shgrav.clm_from_ellipsoid>`)
- Using a surface relief grid ( :func: `clm_from_relief <swiftest.shgrav.clm_from_relief>`). 
  
   * This function is still in development and may not work as expected.
  
- Manually entering the coefficients when adding the central body ( :func: `add_body <swiftest.Simulation.add_body>`)

This is based on ``/spherical_harmonics_cb`` in ``swiftest/examples``.



.. .. toctree::
..    :maxdepth: 2
..    :hidden: