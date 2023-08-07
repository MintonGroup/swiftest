"""
 Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY
 without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import os

# Build the pybindings extension that allows us to run the Fortran driver as a Python module. 
root_dir = 'pybindings'
include_dirs = "/Users/daminton/git/swiftest/apple_install/usr/local/include;/Users/daminton/git/swiftest/apple_install/usr/local/include"
include_dirs = include_dirs.split()
include_dirs.append(root_dir)
link_flags = "-lswiftest  /Users/daminton/git/swiftest/apple_install/usr/local/lib/libnetcdff.a /Users/daminton/git/swiftest/apple_install/usr/local/lib/libnetcdf.a -L/Users/daminton/git/swiftest/apple_install/usr/local/lib -lhdf5_hl -lhdf5 -lm -lz -lbz2 -lxml2 -lcurl"
link_flags = link_flags.split()

pybindings_extension = [Extension('swiftest.bindings',
                         [os.path.join(root_dir,'pybindings.pyx')],
                         extra_compile_args=['-fPIC', '-O3','-fopenmp'],
                         extra_link_args=link_flags,
                         libraries=['gfortran','omp'],
                         include_dirs=include_dirs,
                         )]

setup(name='swiftest',
      version='2023.08.00',
      author='David A. Minton',
      author_email='daminton@purdue.edu',
      url='https://github.itap.purdue.edu/MintonGroup/swiftest',
      ext_modules = cythonize(pybindings_extension),
      packages=find_packages())
