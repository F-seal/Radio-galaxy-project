import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.coordinates import SkyCoord
import os
import scipy

class RadioData():
    def __init__(self,NSIDE=64,flux_cut=1,
                catalog_path='./np_data/catalog/',
                mask_path='./mask/',catalog_name='vlass'):
    # Read catalog and mask, calculate useful constants, and
        self.mask_filename = mask_path + 'mask_'+str(NSIDE) + '.fits'
        catalog_filename = catalog_path + catalog_name + '.npy'
        self.NSIDE = NSIDE
        self.NPIX = hp.nside2npix(NSIDE)
        self.name = catalog_name.upper()
        self.flux_cut = flux_cut
        self.pixel_area = 4*np.pi/self.NPIX
        self.mask = hp.read_map(self.mask_filename,dtype=int)
        self.b_mask = np.array(self.mask,dtype=bool)
        self.total_area = np.sum(self.mask) * self.pixel_area
        self.fsky = self.total_area /(4*np.pi)
        self.map_o,self.flux = RadioData.get_map(self.NSIDE,catalog_filename,flux_cut)
        self.total_number = np.sum(self.map_o*self.mask)
        self.shot_noise = self.total_area/self.total_number
        mean = self.total_number/np.sum(self.mask)
        self.map = (self.map_o-mean)/mean
        self.result_dir = './'+self.name+'/'+str(self.flux_cut) +'/'
        if not os.path.exists(self.result_dir):
            os.makedirs(self.result_dir)
        self.map_filename = self.result_dir + 'map_' + str(self.NSIDE) + '.fits'
        hp.write_map(self.map_filename,self.map,overwrite=True)

    @classmethod
    def get_map(cls,NSIDE,filename,LOW_FLUX):
    # Get number density map with certain flux cut from a catalog 
        Flux = np.load(filename)
        RA_0,DEC_0,flux_0 = Flux[:,0],Flux[:,1],Flux[:,2]
        flux_flag =  (flux_0 > LOW_FLUX)
        RA,DEC,flux = RA_0[flux_flag],DEC_0[flux_flag],flux_0[flux_flag]
        indices = hp.ang2pix(NSIDE,RA,DEC,lonlat=True)
        idx, counts = np.unique(indices, return_counts=True)
        NPIX = hp.nside2npix(NSIDE)
        map = np.full(NPIX,0)
        map[idx] = counts
        return map,flux

    def get_cl(self,method='spice',lmax=-1):
    # Calculate cl from number density map
        out_name = self.result_dir + 'cl_'+ str(self.NSIDE) + '.dat'
        if method == 'spice':
            exe = 'spice -mapfile '+ self.map_filename +' -maskfile '+ self.mask_filename +' -clfile '\
                + out_name +' -nlmax '+ str(lmax)
            print(exe)
            os.system(exe)
        if method == 'Peebles':
            if lmax == -1:
                lmax = 3*self.NSIDE+1
            cl_box = np.zeros(lmax)
            for l in range(lmax):
                a = np.zeros(2*l+1,dtype=complex)
                j = np.zeros(2*l+1,dtype=float)
                cl = 0
                for m in range(-1*l,l+1): 
                    for i in range(self.NPIX):
                        if self.mask[i] == 1:
                            theta,phi = hp.pix2ang(self.NSIDE,i)
                            a[m] +=  scipy.special.sph_harm(m, l, theta, phi) * self.map[i] * self.pixel_area
                            j[m] +=  abs(scipy.special.sph_harm(m, l, theta, phi))**2 * self.pixel_area
                    cl += abs(a[m])**2 / j[m]
                cl_box[l-1] = cl
            np.savetxt(out_name,cl_box)

    def plot_cl(self,show=False,lmax=50):
    # Plot cl as function of l
        if not os.path.exists(self.result_dir + 'cl_'+ str(self.NSIDE) + '.dat'):
            self.get_cl(lmax=lmax)
        f = open(self.result_dir + 'cl_'+ str(self.NSIDE) + '.dat')
        f.readline()
        l = np.arange(lmax)
        cl = np.zeros(lmax)
        for i in range(lmax):
            cl[i] = float(f.readline().split('  ')[-1][:-1]) - self.shot_noise
        plt.scatter(l,cl*1000)
        plt.title(r'$C_l$ of VLASS with flux cut '+str(self.flux_cut)+'mJy')
        plt.xlabel(r'$l$')
        plt.ylabel(r'$C_l*1000$')

        if show:
            plt.show()
        plt.savefig(self.result_dir+'cl.png')
        plt.close()

    def get_dipole(self):
    # calculate cl_1 and the error of cl_1
        if not os.path.exists(self.result_dir + 'cl_'+ str(self.NSIDE) + '.dat'):
            self.get_cl(lmax=1)
        f = open(self.result_dir + 'cl_'+ str(self.NSIDE) + '.dat')
        f.readline()
        f.readline()
        cl = float(f.readline().split('  ')[-1][:-1]) - self.shot_noise
        error = np.sqrt(2/(3*self.fsky))*(cl+self.shot_noise)
        return cl,error

    def get_distribution(self,show=False):
    # get histogram of the distribution of sourcecs
        map = self.map_o[self.b_mask]
        plt.hist(map,bins=20)
        plt.xlim(0,np.max(map))
        plt.xlabel('Number of sources per pixel')
        plt.ylabel('Number of pixels')
        if show:
            plt.show()
        plt.savefig(self.result_dir+'distribution_'+str(self.NSIDE)+'.png')
        plt.close()

    def print_info(self):
    # Print infomation
        print('*'*50)
        print('Catalog name:',self.name)
        print('Flux cut at '+str(self.flux_cut)+'mJy')
        print('Nisde:',self.NSIDE)
        print('Sky coverage:',self.fsky)
        print('Number of sources:',self.total_number)
        print('Shot noise:',self.shot_noise)

    def save_map(self,contrast=True,density=True,mask=False,show=False):
    # Save the maps
        map = self.map
        map_o =np.array(self.map_o,dtype=np.float64)
        if mask:
            map_o[~self.b_mask] = hp.UNSEEN 
            map[~self.b_mask] = hp.UNSEEN
            map = hp.ma(map)
            map_o = hp.ma(map_o)
        if contrast:
            hp.mollview(map,title='Number density contrast map at flux cut '+str(self.flux_cut)+'mJy' )
            if show:
                plt.show()
            if mask:
                plt.savefig(self.result_dir+'contrast_map_masked_'+str(self.NSIDE)+'.png')
            else:
                plt.savefig(self.result_dir+'contrast_map_'+str(self.NSIDE)+'.png')
            plt.close()
        if density:
            hp.mollview(map_o,title='Number density map at flux cut '+str(self.flux_cut)+'mJy')
            if show:
                plt.show()
            if mask:
                plt.savefig(self.result_dir+'number_density_map_masked_'+str(self.NSIDE)+'.png')
            else:
                plt.savefig(self.result_dir+'number_density_map_'+str(self.NSIDE)+'.png')
            plt.close()



if __name__ == '__main__':
    rd = RadioData(flux_cut=50,catalog_name='vlass')
    # rd.get_cl()
    rd.print_info()
    # rd.get_distribution(show=True)
    # rd.save_map(show=True,mask=True)
    # rd.plot_cl(show=True)
    # print(rd.get_dipole())