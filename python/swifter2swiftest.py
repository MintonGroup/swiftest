import numpy as np
import sys
import swiftestio as swio

"""
  Converts initial conditions files from Swifter to Swiftest
"""
G = 6.6743E-11 # Universal gravitational constant in SI units

def read_param(param):
    # Read and parse old Swifter parameter file
    swifterfile =  param['swifterfile']
    swiftestfile = param['swiftestfile']

    # Read param.in file
    print('Reading Swifter file ' + param['swifterfile'])
    fold = open(swifterfile, 'r')
    swifterlines = fold.readlines()
    fold.close()
    for line in swifterlines:
        fields = line.split()
        if len(fields) > 0:
            if ('PL_IN' == fields[0].upper()): param['PL_OLD'] = fields[1]
            if ('TP_IN' == fields[0].upper()): param['TP_OLD'] = fields[1]
            if ('CHK_CLOSE' == fields[0].upper()): param['CHK_CLOSE'] = fields[1].upper()
            if ('RHILL_PRESENT'== fields[0].upper()): param['RHILL_PRESENT'] = fields[1].upper()
            if ('J2' == fields[0].upper()): param['J2'] = float(fields[1])
            if ('J4' == fields[0].upper()): param['J4'] = float(fields[1])

    print('Reading Swiftest file ' + param['swiftestfile'])
    fnew = open(swiftestfile, 'r')
    swiftestlines = fnew.readlines()
    for line in swiftestlines:
        fields = line.split()
        if len(fields) > 0:
            if ('MU2KG' == fields[0].upper()): param['MU2KG'] = float(fields[1])
            if ('DU2M' == fields[0].upper()): param['DU2M'] = float(fields[1])
            if ('TU2S' == fields[0].upper()): param['TU2S'] = float(fields[1])
            if ('CB_IN' == fields[0].upper()): param['CB_NEW'] = fields[1]
            if ('PL_IN' == fields[0].upper()): param['PL_NEW'] = fields[1]
            if ('TP_IN' == fields[0].upper()): param['TP_NEW'] = fields[1]
            if ('CHK_CLOSE' == fields[0].upper()): param['CHK_CLOSE_NEW'] = fields[1].upper()

if __name__ == '__main__':

    param = {'swifterfile' : "param.in",
              'swiftestfile': "param.in"}

    read_param(param)
    MU2KG = param['MU2KG']
    DU2M  = param['DU2M']
    TU2S  = param['TU2S']

    for key, val in param.items():
        print(key, val)

    GU = G / (DU2M ** 3 / (MU2KG * TU2S ** 2))

    cbrad = 6.95700e8 / DU2M

    plnew = open(param['PL_NEW'], 'w')
    print(f'Writing out new PL file: {param["PL_NEW"]}')

    with open(param['PL_OLD'], 'r') as plold:
        line = plold.readline()
        line = line.split("!")[0]  # Ignore comments
        i_list = [i for i in line.split(" ") if i.strip()]
        npl = int(i_list[0])
        print(npl)
        print(npl - 1, file=plnew)
        line = plold.readline()
        #line = line[0].split("!")[0] # Ignore comments
        i_list = [i for i in line.split(" ") if i.strip()]
        GMsun = float(i_list[1]) # Store central body GM for later
        line = plold.readline()  # Ignore the two zero vector lines
        line = plold.readline()
        for n in range(1, npl):  # Loop over planets
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            name = int(i_list[0])
            GMpl = float(i_list[1])
            #print(f'{name} {GMpl} -> {GMpl / GU}')
            print(name, GMpl,file=plnew)
            if param['CHK_CLOSE'] == 'YES':
                line = plold.readline()
                i_list = [i for i in line.split(" ") if i.strip()]
                plrad = float(i_list[0])
                print(plrad)
                print(plrad,file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            xh = float(i_list[0])
            yh = float(i_list[1])
            zh = float(i_list[2])
            print(xh, yh, zh)
            print(xh, yh, zh, file=plnew)
            line = plold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            vx = float(i_list[0])
            vy = float(i_list[1])
            vz = float(i_list[2])
            print(vx, vy, vz)
            print(vx, vy, vz, file=plnew)

    plold.close()
    plnew.close()
    tpnew = open(param['TP_NEW'], 'w')

    print(f'Writing out new TP file: {param["TP_NEW"]}')
    with open(param['TP_OLD'], 'r') as tpold:
        line = tpold.readline()
        line = line.split("!")[0]  # Ignore comments
        i_list = [i for i in line.split(" ") if i.strip()]
        ntp = int(i_list[0])
        print(ntp)
        print(ntp, file=tpnew)
        for n in range(0, ntp):  # Loop over test particles
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            name = int(i_list[0])
            print(name)
            print(name, file=tpnew)
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            xh = float(i_list[0])
            yh = float(i_list[1])
            zh = float(i_list[2])
            print(xh, yh, zh)
            print(xh, yh, zh, file=tpnew)
            line = tpold.readline()
            i_list = [i for i in line.split(" ") if i.strip()]
            vx = float(i_list[0])
            vy = float(i_list[1])
            vz = float(i_list[2])
            print(vx, vy, vz)
            print(vx, vy, vz, file=tpnew)

    tpold.close()
    tpnew.close()

    print(f'Writing out new CB file: {param["CB_NEW"]}')
    # Write out new central body file
    cbnew = open(param['CB_NEW'], 'w')

    print(f'GMsun = {GMsun} -> {GMsun / GU}')
    print(GMsun, file=cbnew)
    print(cbrad)
    print(cbrad, file=cbnew)
    print(param['J2'])
    print(param['J2'], file=cbnew)
    print(param['J4'])
    print(param['J4'], file=cbnew)

    cbnew.close()

