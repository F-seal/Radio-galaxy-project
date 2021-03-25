from RadioData import *

l = np.arange(5,80,5)
cl_v = np.empty(15)
e_v = np.empty(15) 
cl_n = np.empty(15)
e_n = np.empty(15)

for i in range(15):
    flux_cut = 5 * i + 5
    cl_v[i],e_v[i] = RadioData(flux_cut=flux_cut,catalog_name='vlass').get_dipole()
    cl_n[i],e_n[i] = RadioData(flux_cut=flux_cut,catalog_name='nvss').get_dipole()

# plt.errorbar(l,cl_v,yerr=e_v,label='VLASS')
plt.errorbar(l,cl_n,yerr=e_n,label='NVSS')
plt.title(r'${C_l}_1$ at different flux cut')
plt.ylabel(r'$C_l$')
plt.xlabel('Flux cut (mJy)')
plt.legend()
plt.show()
plt.close()