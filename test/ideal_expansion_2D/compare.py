import h5py
import numpy as np

refname = 'sol.0100.ref'
filename = 'sol-000100.cgns'

vars = ['SpeciesDensity1','EnergyInternal','MomentumX','MomentumY']

dataref = [];
datafile = []

for v in vars:
    grp_key = "/Base/Zone 1/FS/" + v + "/ data"
    with h5py.File(refname, 'r') as f:
        dataref = np.array(f[grp_key])

    with h5py.File(filename, 'r') as f:
        datafile = np.array(f[grp_key])


    diff = np.absolute(dataref-datafile)
    reldiff = diff/dataref
    maxdiff = np.amax(diff)
    maxreldiff = np.amax(reldiff)
    print(v)
    print("  Max Difference:          {0:.8e}".format(maxdiff))
    print("  Max Relative Difference: {0:.2f}%".format(maxreldiff))
