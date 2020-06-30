import numpy as np
from matplotlib import pyplot as plt
import os
import sys
from scipy.io import FortranFile

binfilename = "bin.swiftest.dat"

with FortranFile(binfilename, 'r') as f:
    while True:  # Loop until you read the end of file
        try:
            t = f.read_reals(np.float64)  # Try first part of the header
            print(t)
        except:
            break
        if t[0] > tstop: break
        npl = f.read_ints()
        ntp = f.read_ints()
        iout_form = f.read_ints()
        name = f.read_ints()
        px = f.read_reals(np.float64)
        py = f.read_reals(np.float64)
        pz = f.read_reals(np.float64)
        vx = f.read_reals(np.float64)
        vy = f.read_reals(np.float64)
        vz = f.read_reals(np.float64)
        mass = f.read_reals(np.float64)
        radius = f.read_reals(np.float64)
        a = f.read_reals(np.float64)
        e = f.read_reals(np.float64)
        inc = f.read_reals(np.float64)
        yield t, name, mass, radius, np.c_[a, e]