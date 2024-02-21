---
title: 'Swiftest: An \textit{N}-body Integrator for Gravitational Systems'
tags:
  - Python
  - Fortran
  - Astronomy
  - Dynamics
  - N-body
  - Planetary Systems
authors:
  - name: Carlisle Wishard
    affiliation: "1,4"
    orcid: 0009-0001-0733-3268
    equal-contrib: true 
  - name: Jennifer Pouplin
    affiliation: "1,3"
  - name: Jacob Elliott
    affiliation: "1"
  - name: Dana Singh
    affiliation: "5"
  - name: Kaustub Anand
    affiliation: "1, 2"
  - name: David Minton
    affiliation: "1"
    orcid: 0000-0003-1656-9704
    equal-contrib: true 
    corresponding: true 
affiliations:
  - name: Department of Earth, Atmospheric, and Planetary Sciences, Purdue University, USA
    index: 1
  - name: School of Aeronautics and Astronautics, Purdue University, USA
    index: 2
  - name: Department of Physics and Astronomy, Purdue University, USA
    index: 3
  - name: Evansville Museum of Arts, History & Science, USA
    index: 4
  - name : SAIC, Princeton, NJ, USA
    index: 5

date: 06 October 2023
bibliography: paper.bib
---

# Summary
`Swiftest` is a software package designed to model the long-term dynamics of system of bodies in orbit around a dominant central body, such a planetary system around a star, or a satellite system around a planet. The main body of the program is written in Modern Fortran, taking advantage of the object-oriented capabilities included with Fortran 2003 and the parallel capabilities included with Fortran 2008 and Fortran 2018. `Swiftest` also includes a Python package that allows the user to quickly generate input, run simulations, and process output from the simulations. `Swiftest` uses a NetCDF output file format which makes data analysis with the `Swiftest` Python package a streamlined and flexible process for the user.

# Statement of Need
Building off a strong legacy, including its predecessors `Swifter` [@kaufmann_swifter_nodate] and `Swift` [@Levison:1994], `Swiftest` takes the next step in modeling the dynamics of planetary systems by improving the performance and ease of use of software, and by introducing a new collisional fragmentation model. Currently, `Swiftest` includes the four main symplectic integrators included in its predecessors: 

WHM
  : Wisdom-Holman method.  `WHM` is a symplectic integrator that is best suited for cases in which the orbiting bodies do not have close encounters with each other. For details see @Wisdom:1991.

RMVS
  : Regularized Mixed Variable Symplectic. `RMVS` is an extension of `WHM` that handles close approaches between test particles and massive bodies. For details, see @Levison:1994.

HELIO
  : Democratic Heliocentric method.  This is a basic symplectic integrator that uses democratic heliocentric coordinates instead of the Jacobi coordinates used by `WHM`. Like `WHM` it is not suited for simulating close encounters between orbiting bodies. For details, see @Duncan:1998. 

SyMBA
  : Symplectic Massive Body Algorithm. This is an extension of `HELIO` that handles close approaches between massive bodies and any of the other objects in the simulation. It also includes semi-interacting massive bodies that can gravitationally influence fully massive bodies but not each other. This algorithm is described in the @Duncan:1998. See also @Levison:2000.

All of the integrators that are currently implemented are based symplectic integrators, which are designed to model massive bodies (e.g. planets) or test particles (e.g. asteroids) in orbit of a single dominant central body (e.g. the Sun) for very long periods of times (e.g. Gy). The core components of `Swiftest` were inherited from the `Swift`/`Swifter` codebase, but have been somewhat modified for performance. The `Swift`-family of integrators are most comparable to other symplectic integrators, which are among the most popular numerical methods used to study the long-term dynamics of planetary systems. 

The `SyMBA` integrator included in `Swiftest` is most similar to the hybrid symplectic integrator `MERCURY6` [@Chambers:1999], the `MERCURIUS` integrator of `REBOUND` [@Rein:2012;@Rein:2015], and the GPU-enabled hybrid symplectic integrators such as  `QYMSYM` [@Moore:2011] and `GENGA II` [@Grimm:2022], with some important distinctions. The hybrid symplectic integrators typically employ a symplectic method, such as the original WHM method in Jacobi coordinates or the modified method that uses the Democratic Heliocentric coordinates, only when bodies are far from each other relative to their gravitational spheres of influence (some multiple of a Hill's sphere). When bodies approach each other, the integration is smoothly switched to a non-symplectic method, such as Bulirsch-Stoer or IAS15. In contrast, `SyMBA` is a multi-step method that recursively subdivides the step size of bodies undergoing close approach with each other. 

In addition `Swiftest` contains a number of significant enhancements relative to its predecessors:

- It includes a fully symplectic implementation of the post-Newtonian correction (General Relativity) as described in @Saha:1994.
- Output data is stored in NetCDF4 (HDF5) format, making data analysis, cross-platform compatibility, and portability far more robust than the flat binary file formats used by its predecessors.
- It makes use of CPU-based parallelization via OpenMP. It also makes use of SIMD instructions, but these are available only when compiled for target host machine rather than when installed from PyPI via pip.
- Simulations can be run by means of a standalone binary driver program, similar to its predecessors, or from a compiled Python module. Both the standalone driver and Python module link to the same compiled library, so simulations run by either method are identical.
- It comes with an extensive set of Python scripts to help generate simulation initial conditions and post-process simulation results.
- It includes the `Fraggle` collisional fragmentation model that is used to generate collisional fragments on trajectories that conserve linear and angular momentum and lose the appropriate amount of collisional energy for inelastic collisions between massive bodies in `SyMBA` simulations. The collisional outcome is determined using standard methods based on the work of @Leinhardt:2012 and @Stewart:2012.

# Performance
Modeling the behavior of thousands of fully interacting bodies over long timescales is computationally expensive, with typical runs taking weeks or months to complete. The addition of collisional fragmentation can quickly generate hundreds or thousands of new bodies in a short time period, creating further computational challenges for traditional \textit{n}-body integrators. As a result, enhancing computational performance was a key aspect of the development of `Swiftest`. Here we show a comparison between the performance of `Swift`, `Swifter-OMP` (a parallel version of `Swifter`), and `Swiftest` on simulations with 1k, 2k, 8k, and 16k fully interacting bodies. The number of cores dedicated to each run is varied from 1 to 24 to test the parallel performance of each program.

\autoref{fig:performance} shows the results of this performance test. We can see that `Swiftest` outperforms `Swifter-OMP` and `Swift` in each simulation set, even when run in serial. When run in parallel, `Swiftest` shows a significant performance boost when the number of bodies is increased. The improved performance of `Swiftest` compared to `Swifter-OMP` and `Swift` is a critical step forward in \textit{n}-body modeling, providing a powerful tool for modeling the dynamical evolution of planetary systems.

![Performance testing of `Swiftest` on systems of (a) 1k, (b) 2k, (c) 8k, and (d) 16k fully interacting massive bodies. All simulations were run using the \textit{SyMBA} integrator included in `Swift`, `Swifter-OMP` (and earlier attempt to parallelize `Swifter`), and `Swiftest`. Speedup is measured relative to `Swiftest` for 1 CPU (dashed), with an ideal 1:1 speedup shown as an upper limit (dotted). The performance of `Swifter-OMP` is shown in green while the performance of `Swiftest` is shown in blue. All simulations were run on the Purdue University Rosen Center for Advanced Computing Brown Community Cluster. Brown contains 550 Dell compute nodes, with each node containing 2 12-core Intel Xeon Gold Sky Lake processors, resulting in 24 cores per node. Each node has 96 GB of memory. \label{fig:performance}](performance.png)

# Acknowledgements
`Swiftest` was developed at Purdue University and was funded under the NASA Emerging Worlds and Solar System Workings programs. Active development by the Purdue Swiftest Team is ongoing and contributions from the community are highly encouraged.

# References
