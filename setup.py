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

cmake_args = [
            '-DCMAKE_BUILD_TYPE=RELEASE',
            '-DUSE_COARRAY:BOOL=OFF',
            '-DUSE_OPENMP:BOOL=ON',
            ]

with open('version.txt') as version_file:
    version = version_file.read().strip()
    
setup(name='swiftest',
      version=version,
      author='David A. Minton',
      author_email='daminton@purdue.edu',
      url='https://github.itap.purdue.edu/MintonGroup/swiftest',
      python_requires=">3.8",
      license="GPLv3",
      classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha/Development',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
      ],
      keywords='astronomy astrophysics planetary nbody integrator symplectic wisdom-holman',
      cmake_args=cmake_args,
      install_requires= [
            'numpy>=1.24.3',
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
      packages=['swiftest'],
      test_suite="swiftest.tests",
      )
