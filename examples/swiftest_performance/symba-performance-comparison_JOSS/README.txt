Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
This file is part of Swiftest.
Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with Swiftest. 
If not, see: https://www.gnu.org/licenses. 

README.txt

Swiftest Example : Fragmentation
Author           : David Minton
Date             : March 4, 2023

Included in the symba-performance-comparison_JOSS example directory are the following files and directories:

	- README.txt 		    : This file
	Jul27_Test/ 		    : Initial conditions directory
	  ├─ Initial_Conditions_README : Description of initial conditions files
  	  ├─ Jul27_1k_fully.in         : Initial conditions for 1k fully interacting bodies
	  ├─ Jul27_2k_fully.in 	       : Initial conditions for 2k fully interacting bodies 
	  ├─ Jul27_4k_fully.in	       : Initial conditions for 4k fully interacting bodies 
	  ├─ Jul27_8k_fully.in 	       : Initial conditions for 8k fully interacting bodies 
	  ├─ Jul27_16k_fully.in	       : Initial conditions for 16k fully interacting bodies 
	  ├─ Jul27_32k_fully.in	       : Initial conditions for 32k fully interacting bodies 
	swift/
	swifter-omp/
	swiftest/
	- parallel-performance-plots.py    : A Python Script that generates plots that compare the time test results for the different versions of SyMBA

This example is intended to be run with Swiftest, Swifter-OMP, and Swift SyMBA. 
