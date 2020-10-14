import numpy as np
import h5py
import xml.etree.ElementTree as ET
import glob
import os

a = [0,0]
n = [500,150]

s = [2,2]

fs = glob.glob('sol*.xmf')

a.reverse()
n.reverse()
s.reverse()
for f in fs:
    fname = os.path.splitext(f)[0]
    #fname = 'sol-2000'
    print(fname)

    ### read xdmf
    tree = ET.parse(fname+".xmf")
    root = tree.getroot()

    ### get dimensions
    topo = root.find("./Domain/Grid/Topology")
    dstr = topo.get('NumberOfElements')
    dims = [int(s) for s in dstr.split() if s.isdigit()]

    ### compute slice parameters
    ndims = len(dims)
    cn = [ai-bi for ai,bi in zip(n,s)]
    gb = [ai+ni for ai,ni in zip(a,n)]
    cb = [ai+ni for ai,ni in zip(a,cn)]

    gdims = [int((ni)/si) for ni,si in zip(n,s)]
    cdims = [int((ni)/si) for ni,si in zip(cn,s)]

    ### create dimension strings
    cdimsStr = [str(c) for c in cdims]
    gdimsStr = [str(c) for c in gdims]
    grdDims = " ".join(gdimsStr)
    scaDims = " ".join(cdimsStr)
    vecDims = " ".join(cdimsStr)+" "+str(ndims)

    ### Grid ###"
    topo.set('NumberOfElements',grdDims)
    for data in root.findall("./Domain/Grid/Geometry/DataItem"):
        data.set('Dimensions',grdDims)
        data.text = data.text.replace('sol-','slice-sol-')

    ### Vectors ###"
    for data in root.findall(".//*[@AttributeType='Vector']/DataItem"):
        data.set('Dimensions',vecDims)
    for data in root.findall(".//*[@AttributeType='Vector']/DataItem/DataItem"):
        data.set('Dimensions',scaDims)
        data.text = data.text.replace('sol-','slice-sol-')

    ### Scalars ###"
    for data in root.findall(".//*[@AttributeType='Scalar']/DataItem"):
        data.set('Dimensions',scaDims)
        data.text = data.text.replace('sol-','slice-sol-')

    with open('slice-'+fname+'.xmf','wb') as f:
        f.write('<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'.encode('utf-8'))
        tree.write(f,'utf-8')

    #### HDF5
    def foo(name,node):
        if isinstance(node,h5py.Dataset):
             with h5py.File('slice-'+fname+'.h5','a') as f2:
                dset = f[node.name]
                if "Dimension" in node.name:
                    tmpd = dset[a[0]:gb[0]:s[0],a[1]:gb[1]:s[1]]
                    dset2 = f2.create_dataset(node.name,tmpd.shape,dtype='f4',data=dset[a[0]:gb[0]:s[0],a[1]:gb[1]:s[1]])
                else:
                    tmpd = dset[a[0]:cb[0]:s[0],a[1]:cb[1]:s[1]]
                    dset2 = f2.create_dataset(node.name,tmpd.shape,dtype='f4',data=dset[a[0]:cb[0]:s[0],a[1]:cb[1]:s[1]])

    with h5py.File(fname+'.h5','r') as f:
        f.visititems(foo)
