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

from skbuild import setup

with open('version.txt') as version_file:
    version = version_file.read().strip()
    
setup(name='swiftest',
      version=version,
      author='David A. Minton',
      author_email='daminton@purdue.edu',
      url='https://github.itap.purdue.edu/MintonGroup/swiftest',
      python_requires=">3.8",
      license="GPLv3",
      cmake_args=[
            '-DCMAKE_BUILD_TYPE=RELEASE',
            '-DMACHINE_CODE_VALUE="generic"',
            '-DUSE_COARRAY:BOOL=OFF',
            '-DUSE_OPENMP:BOOL=ON',
            '-DBUILD_SHARED_LIBS:BOOL=OFF'
            ],
      install_requires= [
            'numpy>=1.24.3',
            'pandas>=1.5.3',
            'scipy>=1.10.1',
            'xarray>=2022.11.0',
            'dask>=2022.1',
            'bottleneck>=1.3.5',
            'h5netcdf>=1.0.2',
            'netcdf4>=1.6.2',
            'matplotlib>=3.7.1',
            'astropy>=5.1',
            'astroquery>=0.4.6',
            'tqdm>=4.65.0',
      ],
      packages=['swiftest'])
