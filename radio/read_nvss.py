import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd

N_dat = 1773484

f = open('./data/nvss.dat')
data = np.full([N_dat,3],np.nan)

for i in range(N_dat):
    temp = f.readline()
    data[i,0] = float(temp[40:42])* 15. + float(temp[43:45])*15./60. + float(temp[46:51])*15./3600
    sign = (temp[52]=='-')
    data[i,1] = float(temp[53:55]) + float(temp[56:58])/60. +float(temp[59:63])/3600.
    if sign:
        data[i,1] = -1* data[i,1]
    data[i,2] = float(temp[75:83])

np.save('./np_data/catalog/nvss.npy',data)

df = pd.DataFrame(data=data,columns=['RA','DEC','Flux'])
df.to_csv('./data/nvss.scv')
