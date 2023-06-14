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
    orcid: 0009-0001-0733-3268
    equal-contrib: true 
    corresponding: true 
    affiliation: 1
  - name: David Minton
    equal-contrib: true 
    affiliation: 1
  - name: Jennifer Pouplin
    affiliation: "1"
  - name: Jake Elliott
    affiliation: "1"
  - name: Dana Singh
    affiliation: "1"
affiliations:
 - name: Department of Earth, Atmospheric, and Planetary Sciences, Purdue University, USA
   index: 1
date: 05 April 2023
bibliography: paper.bib
---

# Summary

The dynamical evolution of planetary systems is dominated by gravitational interactions between massive bodies. Determining the orbits of massive bodies over long time scales is the first step towards understanding the formation and evolution of planets, moons, asteroids, comets, and more. To model these systems, which often include hundreds or thousands of gravitationally interacting bodies, a numerical tool called an \textit{n}-body integrator is often employed. 

# Statement of Need

`Swiftest` is a software package designed to model gravitationally dominated systems. The main body of the program is written in Modern Fortran, taking advantage of the object-oriented capabilities included with Fortran 2003 and the parallel capabilities included with Fortran 2008 and Fortran 2018. `Swiftest` also includes a Python package that allows the user to quickly generate input and process output from the main integrator. `Swiftest` uses a NetCDF output file format which makes data analysis with the `Swiftest` Python package a streamlined and flexible process for the user. 

Building off a strong legacy, including its predecessors `Swifter` [@Duncan:1998] and `Swift` [@Levison:1994], `Swiftest` takes the next step in modeling gravitationally dominated systems by including collisional fragmentation. Our collisional fragmentation algorithm, `Fraggle` (based on the work of @Leinhardt:2012), is designed to resolve collisions between massive bodies and generate collisional debris. `Swiftest` fully incorporates this debris into the gravitational system, evolving these new bodies along with pre-existing bodies. This allows for a more complete model of the orbital evolution of the system and the growth of massive bodies. 

The combination of modern programming practices, flexible data processing tools, and the latest advances in the field of collisional dynamics make `Swiftest` the ideal tool for studying the formation of planetary systems, the growth of planetary moons, the evolution of asteroid families, and beyond.

# Performance

Modeling the behavior of thousands of fully interacting bodies over long timescales is computationally expensive, with typical runs taking weeks or months to complete. The addition of collisional fragmentation can quickly generate hundreds or thousands of new bodies in a short time period, creating further computational challenges for traditional \textit{n}-body integrators. As a result, enhancing computational performance was a key aspect of the development of `Swiftest`. Here we show a comparison between the performance of `Swift`, `Swifter-OMP` (a parallel version of `Swifter`), and `Swiftest` on simulations with 1k, 2k, 8k, and 16k fully interacting bodies. The number of cores dedicated to each run is varied from 1 to 24 to test the parallel performance of each program.

\autoref{fig:performance} shows the results of this performance test. We can see that `Swiftest` outperforms `Swifter-OMP` and `Swift` in each simulation set, even when run in serial. When run in parallel, `Swiftest` shows a significant performance boost when the number of bodies is increased. The improved performance of `Swiftest` compared to `Swifter-OMP` and `Swift` is a critical step forward in \textit{n}-body modeling, providing a powerful tool for modeling the dynamical evolution of planetary systems.

![Performance testing of `Swiftest` on systems of (a) 1k, (b) 2k, (c) 8k, and (d) 16k fully interacting massive bodies. All simulations were run using the \textit{SyMBA} integrator included in `Swift`, `Swifter-OMP`, and `Swiftest`. Speedup is measured relative to `Swift` (dashed), with an ideal 1:1 speedup relative to `Swiftest`  in serial shown as an upper limit (dotted). The performance of `Swifter-OMP` is shown in green while the performance of `Swiftest` is shown in blue. All simulations were run on the Purdue University Rosen Center for Advanced Computing Brown Community Cluster. Brown contains 550 Dell compute nodes, with each node containing 2 12-core Intel Xeon Gold Sky Lake processors, resulting in 24 cores per node. Each node has 96 GB of memory. \label{fig:performance}](performance.png)

# Acknowledgements

`Swiftest` was developed at Purdue University and was funded under the NASA Emerging Worlds and Solar System Workings programs. Active development by the Purdue Swiftest Team is ongoing and contributions from the community are highly encouraged.

# References
