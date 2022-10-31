---
title: 'Swiftest: An N-body Integrator for Gravitational Systems'
tags:
  - Python
  - Fortran
  - Astronomy
  - Dynamics
  - N-body
  - Planetary Systems
authors:
  - name: Carlisle Wishard
    equal-contrib: true 
    affiliation: 1
  - name: David Minton
    equal-contrib: true 
    affiliation: 1
  - name: Jennifer Pouplin
    corresponding: true 
    affiliation: 2
  - name: Jake Elliott
    corresponding: true 
    affiliation: 2
  - name: Dana Singh
    corresponding: true 
    affiliation: 2
affiliations:
 - name: Department of Earth, Atmospheric, and Planetary Sciences, Purdue University, USA
   index: 1
 - name: Independent Researcher, USA
   index: 2
date: 31 October 2022
bibliography: paper.bib
---

# Summary

The dynamical evolution of planetary systems is dominated by gravitational interactions between massive bodies. Determining the orbits of massive bodies over long time scales is the first step towards understanding the formation and evolution of planets, moons, asteroids, comets, and more. To model these systems, which often include hundreds or thousands of gravitationally interacting bodies, a numerical tool called an N-body integrator is often employed. 

# Statement of need

`Swiftest` is a software package designed to model gravitationally dominated systems. The main body of the program is written in Modern Fortran, taking advantage of the object-oriented capabilities included with Fortran 2003 and the parallel capabilities included with Fortran 2008 and Fortran 2018. `Swiftest` also includes a Python package that allows the user to quickly generate input and process output from the main integrator. `Swiftest` uses a NetCDF output file format which makes data analysis with the `Swiftest` Python package a streamlined and flexible process for the user. 

Building off a strong legacy, including its predecessors `Swifter` [@Duncan:1998] and `Swift` [@Levison:1994], `Swiftest` takes the next step in modeling gravitationally dominated systems by including collisional fragmentation. Our collisional fragmentation algorithm, `Fraggle` (based on the work of @Leinhardt:2012), is designed to resolve collisions between massive bodies and generate collisional debris. `Swiftest` fully incorporates this debris into the gravitational system, evolving these new bodies along with pre-existing bodies. This allows for a more complete model of the orbital evolution of the system and the growth of massive bodies. 

The combination of modern programming practices, flexible data processing tools, and the latest advances in the field of collisional dynamics make `Swiftest` the ideal tool for studying the formation of planetary systems, the growth of planetary moons, the evolution of asteroid families, and beyond.

# References