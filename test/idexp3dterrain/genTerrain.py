#! python3

import numpy as np
import h5py
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from random import random
from random import seed

np.set_printoptions(linewidth=1000)
levels = 7
#amps = [0,0,0,0,0,0,0,0,0,0,0]
#amps = [1.5,0.5,0,0,0,0,0,0,0,0,0]
amp = 1.5
sd = 468

sz = 2**levels+1
z = np.zeros((sz,sz))
seed(sd)

z[0,0] = 1.5
z[0,-1] = 2.5
z[-1,0] = 1.5
z[-1,-1] = 5.0

print('Generating Terrain...')
for s in range(levels):
    amp = amp/((s+1)**(1.0/2.0))
    #amp = amps[s]
    skip = 2**(levels-s)
    for j in range(2**s):
        for i in range(2**s):
            skh = int(skip/2)
            li = i*skip
            lj = j*skip
            ci = int(li+skip/2)
            cj = int(lj+skip/2)
            c = z[li,lj]
            c += z[li+skip,lj]
            c += z[li,lj+skip]
            c += z[li+skip,lj+skip]
            z[ci,cj] = (c/4.0) + 2*amp*random()-amp

            di = ci-skh
            dj = cj
            c  = z[di+skh,dj]
            c += z[di,dj-skh]
            c += z[di,dj+skh]
            if di-skh >= 0:
                c += z[di-skh,dj]
                z[di,dj] = (c/4.0) + 2*amp*random()-amp
            else:
                z[di,dj] = (c/3.0) + 2*amp*random()-amp

            di = ci+skh
            dj = cj
            c = z[di-skh,dj]
            c += z[di,dj-skh]
            c += z[di,dj+skh]
            if di+skh < sz:
                c += z[di+skh,dj]
                z[di,dj] = (c/4.0) + 2*amp*random()-amp
            else:
                z[di,dj] = (c/3.0) + 2*amp*random()-amp

            di = ci
            dj = cj-skh
            c  = z[di+skh,dj]
            c += z[di,dj+skh]
            c += z[di-skh,dj]
            if dj-skh >= 0:
                c += z[di,dj-skh]
                z[di,dj] = (c/4.0) + 2*amp*random()-amp
            else:
                z[di,dj] = (c/3.0) + 2*amp*random()-amp

            di = ci
            dj = cj+skh
            c = z[di-skh,dj]
            c += z[di,dj-skh]
            c += z[di+skh,dj]
            if dj+skh < sz:
                c += z[di,dj+skh]
                z[di,dj] = (c/4.0) + 2*amp*random()-amp
            else:
                z[di,dj] = (c/3.0) + 2*amp*random()-amp

print("Writing data...")
hf = h5py.File('terrain.h5','w')
hf.create_dataset('Height',data=z)
hf.close()

print("Plotting...")
x=np.linspace(0,10,sz)
y=np.linspace(0,10,sz)
x,y=np.meshgrid(x,y)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(x,y,z,cmap=cm.terrain,edgecolor='none',rstride=2,cstride=2)
ax.set_zlim(0,10)
plt.show()
