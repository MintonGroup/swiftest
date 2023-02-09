#!/usr/bin/env python3

"""
 Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
 This file is part of Swiftest.
 Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with Swiftest. 
 If not, see: https://www.gnu.org/licenses. 
"""

"""
Generates and runs a set of Swiftest input files and a set of Swifter input files from initial conditions with the 
SyMBA integrator. 

Input
------
None.

Output
------
param.swiftest.in      : An ASCII file containing the parameters for the Swiftest simulation.
param.swifter.in       : An ASCII file containing the parameters for the Swifter simulation.
cb.swiftest.in         : An ASCII file containing the central body data for the Swiftest simulation.
pl.swifter.in          : An ASCII file containing the massive body data for the Swifter simulation.
pl.swiftest.in         : An ASCII file containing the massive body data for the Swiftest simulation.
tp.swifter.in          : An ASCII file containing the test particle data for the Swifter simulation.
tp.swiftest.in         : An ASCII file containing the test particle data for the Swiftest simulation.

"""

import numpy as np
import swiftest
import swiftest.io as swio
import astropy.constants as const
import sys
import xarray as xr

# Both codes use the same tp input file
tpin = "tp.in"

swifter_input  = "param.swifter.in"
swifter_pl     = "pl.swifter.in"
swifter_tp     = "tp.swifter.in"
swifter_bin    = "bin.swifter.dat"
swifter_enc    = "enc.swifter.dat"

swiftest_input = "param.swiftest.in"
swiftest_pl    = "pl.swiftest.in"
swiftest_tp    = "tp.swiftest.in"
swiftest_cb    = "cb.swiftest.in"
swiftest_bin   = "bin.swiftest.nc"
swiftest_enc   = "enc.swiftest.dat"
swiftest_dis  = "discard.swiftest.dat"

sim = swiftest.Simulation()

sim.param['T0'] = 0.0
sim.param['DT'] = 1.0 
sim.param['TSTOP'] = 365.25e1
sim.param['ISTEP_OUT']  = 10
sim.param['ISTEP_DUMP'] = 10
sim.param['CHK_QMIN_COORD'] = "HELIO"
sim.param['CHK_QMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_QMIN_RANGE'] = f"{swiftest.RSun / swiftest.AU2M} 1000.0"
sim.param['CHK_RMIN'] = swiftest.RSun / swiftest.AU2M
sim.param['CHK_RMAX'] = 1000.0
sim.param['CHK_EJECT'] = 1000.0
sim.param['OUT_STAT'] = "UNKNOWN"
sim.param['GR'] = 'NO'
sim.param['CHK_CLOSE'] = 'YES'
sim.param['RHILL_PRESENT'] = 'YES'
sim.param['GMTINY'] = 1.0e-12
sim.param['IN_FORM'] = 'XV'

sim.param['MU2KG'] = swiftest.MSun
sim.param['TU2S'] = swiftest.JD2S
sim.param['DU2M'] = swiftest.AU2M

bodyid = {
   "Sun": 0,
   "Mercury": 1,
   "Venus": 2,
   "Earth": 3,
   "Mars": 4,
   "Jupiter": 5,
   "Saturn": 6,
   "Uranus": 7,
   "Neptune": 8,
}

for name, id in bodyid.items():
   sim.add(name, idval=id)

ds = sim.ds
cb = ds.sel(id=0)
pl = ds.where(ds.id > 0, drop=True)
npl = pl.id.size

ntp = 16
dims = ['time', 'id', 'vec']
tp = []
t = np.array([0.0])
sim.param['OUT_FORM'] = "XV"
clab, plab, tlab = swio.make_swiftest_labels(sim.param)

# For each planet, we will initialize a pair of test particles. One on its way in, and one on its way out. We will also initialize two additional particles that don't encounter anything
tpidlist = np.arange(9,9+ntp)
tpnames = []
for i in tpidlist:
   tpnames.append(f"SmallBody{i+1:0.4g}")
tpxv1 = np.empty((6))
tpxv2 = np.empty((6))

p1 = []
p2 = []
p3 = []
p4 = []
p5 = []
p6 = []

for i in pl.id:
    pli = pl.sel(id=i)
    rstart = 2 * np.double(pli['radius'])  # Start the test particles at a multiple of the planet radius away
    vstart = 1.5 * np.sqrt(2 * np.double(pli['Gmass'])  / rstart)  # Start the test particle velocities at a multiple of the escape speed
    xvstart = np.array([rstart / np.sqrt(2.0), rstart / np.sqrt(2.0), 0.0, vstart, 0.0, 0.0])
    # The positions and velocities of each pair of test particles will be in reference to a planet
    plvec = np.array([np.double(pli['xhx']),
                      np.double(pli['xhy']),
                      np.double(pli['xhz']),
                      np.double(pli['vhx']),
                      np.double(pli['vhy']),
                      np.double(pli['vhz'])])
    tpxv1 = plvec + xvstart
    tpxv2 = plvec - xvstart
    p1.append(tpxv1[0])
    p1.append(tpxv2[0])
    p2.append(tpxv1[1])
    p2.append(tpxv2[1])
    p3.append(tpxv1[2])
    p3.append(tpxv2[2])
    p4.append(tpxv1[3])
    p4.append(tpxv2[3])
    p5.append(tpxv1[4])
    p5.append(tpxv2[4])
    p6.append(tpxv1[5])
    p6.append(tpxv2[5])


tvec = np.vstack([p1, p2, p3, p4, p5, p6])
tpframe = np.expand_dims(tvec.T, axis=0)
tpxr = xr.DataArray(tpframe, dims = dims, coords = {'time' : t, 'id' : tpidlist, 'vec' : tlab})
tpxr = tpxr.assign_coords(name=('id', tpnames))

tp = [tpxr]
tpda = xr.concat(tp,dim='time')
tpds = tpda.to_dataset(dim = 'vec')

sim.ds = xr.combine_by_coords([sim.ds, tpds])
swio.swiftest_xr2infile(sim.ds, sim.param)

sim.param['PL_IN'] = swiftest_pl
sim.param['TP_IN'] = swiftest_tp
sim.param['CB_IN'] = swiftest_cb
sim.param['BIN_OUT'] = swiftest_bin
sim.param['ENC_OUT'] = swiftest_enc
sim.param['DISCARD_OUT'] = swiftest_dis  
sim.save(swiftest_input)

sim.param['PL_IN'] = swifter_pl
sim.param['TP_IN'] = swifter_tp
sim.param['BIN_OUT'] = swifter_bin
sim.param['ENC_OUT'] = swifter_enc
sim.param['OUT_TYPE'] = "REAL8"
sim.save(swifter_input, codename="Swifter")
