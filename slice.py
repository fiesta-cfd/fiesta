import numpy as np
import h5py

s = 2

def foo(name,node):
    if isinstance(node,h5py.Dataset):
         with h5py.File('slice-2000.h5','a') as f2:
            dset = f[node.name]
            tmpd = dset[::s,::s]
            dset2 = f2.create_dataset(node.name,tmpd.shape,dtype='f4',data=dset[::s,::s])
            print(node.name,dset.shape,dset2.name,dset2.shape,dset2.id)

with h5py.File('sol-2000.h5','r') as f:
    x = f["Grid/Dimension0"].shape[0]
    y = f["Grid/Dimension0"].shape[1]
    print(x,y)
    
    f.visititems(foo)

with h5py.File('test.h5','w') as f3:
    f3.create_group('mygroup')
    f3.create_group('mygroup2')
