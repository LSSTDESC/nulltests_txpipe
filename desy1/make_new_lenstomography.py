# Make new lens tomography file so that the desy1 example runs
# Mainly changes to column names so that its compatible with TXPipe as of 09/14/23

import h5py
import numpy as np 
import pylab as mplot
import astropy.io.fits as pf

RM_des = pf.open('/global/cfs/cdirs/lsst/groups/WL/users/yomori/repo/nulltests_txpipe/desy1$ l /global/cfs/cdirs/lsst/groups/WL/projects/txpipe-sys-tests/des-y1/DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits')
print(RM_des[1].header.keys)

with h5py.File('/global/cfs/cdirs/lsst/groups/WL/projects/txpipe-sys-tests/des-y1/lens_tomography_catalog_desy1_RM_091423.h5', 'w') as f:

    f.create_group("provenance")
    f.create_group("tomography")

    z = RM_des[1].data['ZREDMAGIC']
    tomo_array = np.zeros(len(z)) - 1
    tomo_array[(z>=0.15)*(z<0.3)] = 0
    tomo_array[(z>=0.3)*(z<0.45)] = 1
    tomo_array[(z>=0.45)*(z<0.6)] = 2
    tomo_array[(z>=0.6)*(z<0.75)] = 3
    tomo_array[(z>=0.75)*(z<0.9)] = 4
    
    f['tomography/bin'] = tomo_array.copy()
    f['tomography/counts'] = np.array([len(tomo_array[tomo_array==0]),len(tomo_array[tomo_array==1]),len(tomo_array[tomo_array==2]),len(tomo_array[tomo_array==3]),len(tomo_array[tomo_array==4])])
    f['tomography/counts_2d'] = [len(tomo_array[tomo_array>=0])]
    f['tomography/lens_weight'] = RM_des[1].data['weight']   
    
    f['tomography'].attrs['zbin_edges'] = np.array([0.15, 0.3 , 0.45, 0.6 , 0.75, 0.9 ])
    f['tomography'].attrs['nbin'] = 5