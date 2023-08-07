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

# Build the pydriver extension that allows us to run the Fortran driver as a Python module. 
root_dir = 'pydriver'
include_dir = os.path.join(root_dir,'include')
lib_dir = os.path.join(root_dir,'lib')
pydriver_extension = [Extension('swiftest.pydriver',
                         [os.path.join(root_dir,'pydriver.pyx')],
                         extra_compile_args=['-fPIC', '-O3'],
                         library_dirs=[lib_dir],
                         libraries=['swiftest','netcdff','netcdf','hdf5_hl','hdf5','m','z'],
                         include_dirs=[include_dir],
                         )]

setup(name='swiftest',
      version='2023.08.00',
      author='David A. Minton',
      author_email='daminton@purdue.edu',
      url='https://github.itap.purdue.edu/MintonGroup/swiftest',
      ext_modules = cythonize(pydriver_extension),
      packages=find_packages())
