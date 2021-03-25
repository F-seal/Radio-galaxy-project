import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

data = np.load('./npy_data/other/Total_flux.npy')

x = np.radians(data[:,0]-180)
y = np.radians(data[:,1])
z = data[:,2]


cm = plt.cm.get_cmap('RdYlBu')
fig,ax = plt.subplots(1,1,figsize=(20,20))
ax = plt.subplot(111,projection='aitoff')
ax = plt.scatter(x,y,s=0.001,c=z)
fig.colorbar(ax,fraction=0.02,pad=0.05)
plt.grid(True)
plt.savefig('./result/flux.png')