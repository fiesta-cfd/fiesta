import numpy as np
import h5py
import xml.etree.ElementTree as ET

a = [10,20]
b = [200,100]

s = [2,2]

fname = 'sol-2000'

### xdmf
tree = ET.parse(fname+".xmf")
root = tree.getroot()

print("### Grid ###")
topo = root.find("./Domain/Grid/Topology")
dstr = topo.get('NumberOfElements')
dims = [int(s) for s in dstr.split() if s.isdigit()]
ndim = len(dims)

print(dims)

for data in root.findall("./Domain/Grid/Geometry/DataItem"):
    print(data.get('Dimensions'))

print("### Vectors ###")
for data in root.findall(".//*[@AttributeType='Vector']/DataItem"):
    print(data.get('Dimensions'))

print("### Scalars ###")
for data in root.findall(".//*[@AttributeType='Scalar']/DataItem"):
    print(data.get('Dimensions'))

topo.set('NumberOfElements','200 300')
with open('slice-'+fname+'.xmf','wb') as f:
    f.write('<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'.encode('utf-8'))
    tree.write(f,'utf-8')

#### HDF5
def foo(name,node):
    if isinstance(node,h5py.Dataset):
         with h5py.File('slice-'+fname+'.h5','a') as f2:
            dset = f[node.name]
            tmpd = dset[0:dims[0]:s[0],0:dims[1]:s[1]]
            dset2 = f2.create_dataset(node.name,tmpd.shape,dtype='f4',data=dset[0:dims[0]:s[0],0:dims[1]:s[1]])
            print(node.name,dset.shape,dset2.name,dset2.shape,dset2.id)

with h5py.File(fname+'.h5','r') as f:
    x = f["Grid/Dimension0"].shape[0]
    y = f["Grid/Dimension0"].shape[1]
    print(x,y)
    
    f.visititems(foo)
