import numpy as np
import sys

"""
  Converts initial conditions files from Swifter to Swiftest
"""
G = 6.6743E-11 # Universal gravitational constant in SI units

def read_config(config):
    # Read and parse old Swifter parameter file
    swifterfile =  config['swifterfile']
    swiftestfile = config['swiftestfile']

    # Read config.in file
    print('Reading Swifter file ' + config['swifterfile'])
    fold = open(swifterfile, 'r')
    swifterlines = fold.readlines()
    fold.close()
    for line in swifterlines:
        fields = line.split()
        if len(fields) > 0:
            if ('PL_IN' == fields[0].upper()): config['PL_OLD'] = fields[1]
            if ('TP_IN' == fields[0].upper()): config['TP_OLD'] = fields[1]
            if ('CHK_CLOSE' == fields[0].upper()): config['CHK_CLOSE'] = fields[1].upper()
            if ('RHILL_PRESENT'== fields[0].upper()): config['RHILL_PRESENT'] = fields[1].upper()
            if ('J2' == fields[0].upper()): config['J2'] = float(fields[1])
            if ('J4' == fields[0].upper()): config['J4'] = float(fields[1])

    print('Reading Swiftest file ' + config['swiftestfile'])
    fnew = open(swiftestfile, 'r')
    swiftestlines = fnew.readlines()
    for line in swiftestlines:
        fields = line.split()
        if len(fields) > 0:
            if ('MU2KG' == fields[0].upper()): config['MU2KG'] = float(fields[1])
            if ('DU2M' == fields[0].upper()): config['DU2M'] = float(fields[1])
            if ('TU2S' == fields[0].upper()): config['TU2S'] = float(fields[1])
            if ('CB_IN' == fields[0].upper()): config['CB_NEW'] = fields[1]
            if ('PL_IN' == fields[0].upper()): config['PL_NEW'] = fields[1]
            if ('TP_IN' == fields[0].upper()): config['TP_NEW'] = fields[1]
            if ('CHK_CLOSE' == fields[0].upper()): config['CHK_CLOSE_NEW'] = fields[1].upper()

if __name__ == '__main__':

    config = {'swifterfile' : "param.in",
              'swiftestfile': "config.in"}

    read_config(config)
    MU2KG = config['MU2KG']
    DU2M  = config['DU2M']
    TU2S  = config['TU2S']

    for key, val in config.items():
        print(key, val)

    GU = G / (DU2M ** 3 / (MU2KG * TU2S ** 2))

    cbrad = 6.95700e8 / DU2M

    plnew = open(config['PL_NEW'], 'w')
    print(f'Writing out new PL file: {config["PL_NEW"]}')

    with open(config['PL_OLD'], 'r') as plold:
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
            print(f'{name} {GMpl} -> {GMpl / GU}')
            print(name, GMpl / GU,file=plnew)
            if config['CHK_CLOSE'] == 'YES':
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
    tpnew = open(config['TP_NEW'], 'w')

    print(f'Writing out new TP file: {config["TP_NEW"]}')
    with open(config['TP_OLD'], 'r') as tpold:
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

    print(f'Writing out new CB file: {config["CB_NEW"]}')
    # Write out new central body file
    cbnew = open(config['CB_NEW'], 'w')

    print(f'GMsun = {GMsun} -> {GMsun / GU}')
    print(GMsun / GU, file=cbnew)
    print(cbrad)
    print(cbrad, file=cbnew)
    print(config['J2'])
    print(config['J2'], file=cbnew)
    print(config['J4'])
    print(config['J4'], file=cbnew)

    cbnew.close()

