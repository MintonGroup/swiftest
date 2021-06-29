import numpy as np
import swiftest.io as swio
import astropy.constants as const
import sys
import xarray as xr

# Both codes use the same tp input file
tpin = "tp.in"

swifter_input  = "param.swifter.in"
swifter_pl     = "pl.swifter.in"
swifter_bin    = "bin.swifter.dat"
swifter_enc    = "enc.swifter.dat"

swiftest_input = "param.swiftest.in"
swiftest_pl    = "pl.swiftest.in"
swiftest_cb    = "cb.swiftest.in"
swiftest_bin   = "bin.swiftest.dat"
swiftest_enc   = "enc.swiftest.dat"

#Values from JPL Horizons
AU2M = np.longdouble(const.au.value)
GMSunSI = np.longdouble(const.GM_sun.value)
Rsun = np.longdouble(const.R_sun.value)
GC = np.longdouble(const.G.value)
JD = 86400
year = np.longdouble(365.25 * JD)
c = np.longdouble(299792458.0)

MU2KG    = np.longdouble(GMSunSI / GC) #Conversion from mass unit to kg
DU2M     = np.longdouble(AU2M)         #Conversion from radius unit to centimeters
TU2S     = np.longdouble(JD)         #Conversion from time unit to seconds
GU       = np.longdouble(GC / (DU2M**3 / (MU2KG * TU2S**2)))
GMSun = np.longdouble(GMSunSI / (DU2M**3 / TU2S**2))

t_0	  = 0 # simulation start time
deltaT	= 1.00 * JD / TU2S   # simulation step size
end_sim = year / TU2S # simulation end time
t_print = deltaT  #output interval to print results

iout = int(np.ceil(t_print / deltaT))
rmin = Rsun / DU2M
rmax = 1000.0

sys.stdout = open(swiftest_input, "w")
print(f'! VERSION      Swiftest input file generated using init_cond.py')
print(f'T0             {t_0} ')
print(f'TSTOP          {end_sim}')
print(f'DT             {deltaT}')
print(f'CB_IN          {swiftest_cb}')
print(f'PL_IN          {swiftest_pl}')
print(f'TP_IN          {tpin}')
print(f'IN_TYPE        ASCII')
print(f'ISTEP_OUT      {iout:d}')
print(f'ISTEP_DUMP     {iout:d}')
print(f'BIN_OUT        {swiftest_bin}')
print(f'OUT_TYPE       REAL8')
print(f'OUT_FORM       XV')
print(f'OUT_STAT       REPLACE')
print(f'CHK_CLOSE      yes')
print(f'CHK_RMIN       {rmin}')
print(f'CHK_RMAX       {rmax}')
print(f'CHK_EJECT      {rmax}')
print(f'CHK_QMIN       {rmin}')
print(f'CHK_QMIN_COORD HELIO')
print(f'CHK_QMIN_RANGE {rmin} {rmax}')
print(f'ENC_OUT        {swiftest_enc}')
print(f'EXTRA_FORCE    no')
print(f'BIG_DISCARD    no')
print(f'ROTATION       no')
print(f'GR             no')
print(f'MU2KG          {MU2KG}')
print(f'DU2M           {DU2M}')
print(f'TU2S           {TU2S}')
sys.stdout = sys.__stdout__
param = swio.read_swiftest_param(swiftest_input)

# Dates to fetch planet ephemerides from JPL Horizons
tstart = '2021-06-15'
ds = swio.solar_system_pl(param, tstart)
cb = ds.sel(id=0)
pl = ds.where(ds.id > 0, drop=True)
npl = pl.id.size

ntp = 18
dims = ['time', 'id', 'vec']
tp = []
t = np.array([0.0])
clab, plab, tlab = swio.make_swiftest_labels(param)

# For each planet, we will initialize a pair of test particles. One on its way in, and one on its way out. We will also initialize two additional particles that don't encounter anything
tpnames = np.arange(101, 101 + ntp)
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
    rstart = 2 * np.double(pli['Radius'])  # Start the test particles at a multiple of the planet radius away
    vstart = 1.5 * np.sqrt(2 * np.double(pli['Mass'])  / rstart)  # Start the test particle velocities at a multiple of the escape speed
    xvstart = np.array([rstart / np.sqrt(2.0), rstart / np.sqrt(2.0), 0.0, vstart, 0.0, 0.0])
    # The positions and velocities of each pair of test particles will be in reference to a planet
    plvec = np.array([np.double(pli['px']),
                      np.double(pli['py']),
                      np.double(pli['pz']),
                      np.double(pli['vx']),
                      np.double(pli['vy']),
                      np.double(pli['vz'])])
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
tpxr = xr.DataArray(tpframe, dims = dims, coords = {'time' : t, 'id' : tpnames, 'vec' : tlab})

tp = [tpxr]
tpda = xr.concat(tp,dim='time')
tpds = tpda.to_dataset(dim = 'vec')

ds = xr.combine_by_coords([ds, tpds])
swio.swiftest_xr2infile(ds, param)

# Swifter PL file
plfile = open(swifter_pl, 'w')
print(npl + 1, file=plfile)
print(0,GMSun,file=plfile)
print('0.0 0.0 0.0',file=plfile)
print('0.0 0.0 0.0',file=plfile)
for i in pl.id:
    pli = pl.sel(id=i)
    print(f"{int(i)} {pli['Mass'].values[0]} {pli['Rhill'].values[0]}", file=plfile)
    print(f"{pli['Radius'].values[0]}", file=plfile)
    print(f"{pli['px'].values[0]} {pli['py'].values[0]} {pli['pz'].values[0]}", file=plfile)
    print(f"{pli['vx'].values[0]} {pli['vy'].values[0]} {pli['vz'].values[0]}", file=plfile)
plfile.close()

# Swifter parameter file
sys.stdout = open(swifter_input, "w")
print(f"! VERSION     Swifter input file generated using init_cond.py")
print(f"T0            {t_0} ")
print(f"TSTOP         {end_sim}")
print(f"DT            {deltaT}")
print(f"PL_IN         {swifter_pl}")
print(f"TP_IN         {tpin}")
print(f"IN_TYPE       ASCII")
print(f"ISTEP_OUT     {iout:d}")
print(f"ISTEP_DUMP    {iout:d}")
print(f"BIN_OUT       {swifter_bin}")
print(f"OUT_TYPE      REAL8")
print(f"OUT_FORM      XV")
print(f"OUT_STAT      UNKNOWN")
print(f"J2            {param['J2']}")
print(f"J4            {param['J4']}")
print(f"CHK_CLOSE     yes")
print(f"CHK_RMIN      {rmin}")
print(f"CHK_RMAX      {rmax}")
print(f"CHK_EJECT     {rmax}")
print(f"CHK_QMIN      {rmin}")
print(f"CHK_QMIN_COORD HELIO")
print(f"CHK_QMIN_RANGE {rmin} {rmax}")
print(f"ENC_OUT        {swifter_enc}")
print(f"EXTRA_FORCE    no")
print(f"BIG_DISCARD    no")
print(f"RHILL_PRESENT  yes")
sys.stdout = sys.__stdout__

