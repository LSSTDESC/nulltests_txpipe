'''
This code takes the public star catalogs and converts
them into h5 format that can be ingested in TXPipe.
'''
import os,sys
import h5py 
import numpy as np
from astropy.io import fits

def get_T_e1e2(d,colprefix,idx):
    ixx = d[1].data[colprefix+'_ixx'][idx]
    iyy = d[1].data[colprefix+'_iyy'][idx]
    ixy = d[1].data[colprefix+'_ixy'][idx]
    T   = ixx + iyy
    e1  = 0.5*(ixx-iyy)/T
    e2  = 0.5*2*ixy/T
    # The factor of 0.5 comes from eq 3 of 
    # https://arxiv.org/abs/1705.06745 which is 
    # a conversion factor from e1,e2 -> g1,g2
    # see also eq.4.11 of https://arxiv.org/abs/astro-ph/9912508
    return T,e1,e2

fieldi    = int(sys.argv[1])


if fieldi==0:
    fieldlist = ['GAMA09H','GAMA15H','HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']

    ra  = np.array([])
    dec = np.array([])
    T_meas   = np.array([])
    T_model  = np.array([])
    e1_meas  = np.array([])
    e2_meas  = np.array([])
    e1_model = np.array([])
    e2_model = np.array([])
    idxr     = np.array([])
    idxu     = np.array([])

    for field in fieldlist:
        d = fits.open('/global/cfs/cdirs/lsst/groups/WL/users/yomori/scratch/HSC/hsc_stars/%s_stars.fits'%field)

        # If 'icalib_psf_used'==True, the star was used to model the psf
        # If 'icalib_psf_used'==False the star is used as a test ("reserved star")
        # according to HSC paper the 20% of the stars are reserved
        idx_u = np.where(d[1].data['icalib_psf_used']==True)[0]
        idx_r = np.where(d[1].data['icalib_psf_used']==False)[0]

        print("Processing: %s"%field)
        print("Total number of stars: %d"%(len(idx_u)+len(idx_r) ) )
        print("Number of stars used : %d"%(len(idx_u)))
        print("Number of stars reserved : %d"%(len(idx_r)))

        print("Computing T,e1,e2 from ixx,iyy,ixy")
        rT_model, re1_model, re2_model = get_T_e1e2(d,'ishape_sdss_psf',idx_r)
        rT_meas, re1_meas, re2_meas    = get_T_e1e2(d,'ishape_sdss',idx_r);
        #re1_meas = re1_meas-np.mean(re1_meas)
        #re2_meas = re2_meas-np.mean(re2_meas)

        uT_model, ue1_model, ue2_model = get_T_e1e2(d,'ishape_sdss_psf',idx_u)
        uT_meas, ue1_meas, ue2_meas    = get_T_e1e2(d,'ishape_sdss',idx_u);
        #ue1_meas = ue1_meas-np.mean(ue1_meas)
        #ue2_meas = ue2_meas-np.mean(ue2_meas)

        ra  = np.concatenate([ra , d[1].data['ira'][idx_r] , d[1].data['ira'][idx_u] ])        
        dec = np.concatenate([dec, d[1].data['idec'][idx_r], d[1].data['idec'][idx_u] ])
        T_meas   = np.concatenate([T_meas,rT_meas,uT_meas])
        T_model  = np.concatenate([T_model,rT_model,uT_model])
        e1_meas  = np.concatenate([e1_meas,re1_meas,ue1_meas])
        e2_meas  = np.concatenate([e2_meas,re2_meas,ue2_meas])
        e1_model = np.concatenate([e1_model,re1_model,ue1_model])
        e2_model = np.concatenate([e2_model,re2_model,ue2_model])
        idxr     = np.concatenate([idxr,np.ones(len(idx_r)),np.zeros(len(idx_u))])
        idxu     = np.concatenate([idxu,np.ones(len(idx_r)),np.zeros(len(idx_u))])

    f=h5py.File('./star_catalog_hscy1_allfields.h5', 'w')
    f.create_group("stars")
    f['stars/ra']          = ra
    f['stars/dec']         = dec
    f['stars/measured_T']  = T_meas
    f['stars/measured_e1'] = e1_meas
    f['stars/measured_e2'] = e2_meas
    f['stars/model_T']     = T_model
    f['stars/model_e1']    = e1_model
    f['stars/model_e2']    = e2_model
    f['stars/calib_psf_reserved'] = idxr
    f['stars/calib_psf_used']     = idxu
    f.close()



 
else:

    fieldlist = [None,'GAMA09H','GAMA15H','HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']
    field     = fieldlist[fieldi] 

    d = fits.open('/global/cfs/cdirs/lsst/groups/WL/users/yomori/scratch/HSC/hsc_stars/%s_stars.fits'%field)

    # If 'icalib_psf_used'==True, the star was used to model the psf 
    # If 'icalib_psf_used'==False the star is used as a test ("reserved star")
    # according to HSC paper the 20% of the stars are reserved
    idx_u = np.where(d[1].data['icalib_psf_used']==True)[0]
    idx_r = np.where(d[1].data['icalib_psf_used']==False)[0]

    print("Processing: %s"%field)
    print("Total number of stars: %d"%(len(idx_u)+len(idx_r) ) )
    print("Number of stars used : %d"%(len(idx_u)))
    print("Number of stars reserved : %d"%(len(idx_r)))

    print("Computing T,e1,e2 from ixx,iyy,ixy")
    ixx_meas  = d[1].data['ishape_sdss_ixx'][idx_r]
    iyy_meas  = d[1].data['ishape_sdss_iyy'][idx_r] 
    ixy_meas  = d[1].data['ishape_sdss_ixy'][idx_r]
    rT_meas    = ixx_meas + iyy_meas
    re1_meas   = (ixx_meas-iyy_meas)/rT_meas
    re2_meas   = 2*ixy_meas/rT_meas

    ixx_model = d[1].data['ishape_sdss_psf_ixx'][idx_r]
    iyy_model = d[1].data['ishape_sdss_psf_iyy'][idx_r] 
    ixy_model = d[1].data['ishape_sdss_psf_ixy'][idx_r]
    rT_model   = ixx_model + iyy_model
    re1_model  = (ixx_model-iyy_model)/rT_model
    re2_model  = 2*ixy_model/rT_model

    print("Computing T,e1,e2 from ixx,iyy,ixy")
    ixx_meas  = d[1].data['ishape_sdss_ixx'][idx_u]
    iyy_meas  = d[1].data['ishape_sdss_iyy'][idx_u]
    ixy_meas  = d[1].data['ishape_sdss_ixy'][idx_u]
    uT_meas    = ixx_meas + iyy_meas
    ue1_meas   = (ixx_meas-iyy_meas)/uT_meas
    ue2_meas   = 2*ixy_meas/uT_meas

    ixx_model = d[1].data['ishape_sdss_psf_ixx'][idx_u]
    iyy_model = d[1].data['ishape_sdss_psf_iyy'][idx_u]
    ixy_model = d[1].data['ishape_sdss_psf_ixy'][idx_u]
    uT_model   = ixx_model + iyy_model
    ue1_model  = (ixx_model-iyy_model)/uT_model
    ue2_model  = 2*ixy_model/uT_model


    f=h5py.File('./star_catalog_hscy1_%s.h5'%field, 'w')
    f.create_group("stars")
    f['stars/ra']          = np.concatenate([d[1].data['ira'][idx_r],d[1].data['ira'][idx_u]])
    f['stars/dec']         = np.concatenate([d[1].data['idec'][idx_r],d[1].data['idec'][idx_u]])
    f['stars/measured_T']  = np.concatenate([rT_meas,uT_meas])
    f['stars/measured_e1'] = np.concatenate([-re1_meas,-ue1_meas])
    f['stars/measured_e2'] = np.concatenate([re2_meas,ue2_meas])
    f['stars/model_T']     = np.concatenate([rT_model,uT_model])
    f['stars/model_e1']    = np.concatenate([-re1_model,-ue1_model])
    f['stars/model_e2']    = np.concatenate([re2_model,ue2_model]) 
    f['stars/calib_psf_reserved'] = np.concatenate([np.ones(len(idx_r)),np.zeros(len(idx_u))])
    f['stars/calib_psf_used']     = np.concatenate([np.zeros(len(idx_u)),np.ones(len(idx_r))])
    f.close()

'''
# COLUMNS IN THE DES STAR CATALOG
'calib_psf_reserved',
'calib_psf_used',
'ccd',
'ccd_x',
'ccd_y',
'dec',
'extendedness',
'fov_x',
'fov_y',
'mag',
'measured_T',
'measured_e1',
'measured_e2',
'model_T',
'model_e1',
'model_e2',
'ra'
''';


'''
# COLUMNS IN THE HSC STAR CATALOG https://arxiv.org/abs/1705.06745
dtype=(numpy.record,
('object_id', '>i8'),
('parent_id', '>i8'),
('ira', '>f8'),                             # right ascension (J2000.0) measured in i−band
('idec', '>f8'),                            # declination ascension (J2000.0) measured in i−band
('imag_psf', '>f4'),                        # flux measured by a fit to the PSF model (mag
('imag_psf_err', '>f4'),                    # uncertainty for flux.psf
('iflux_psf', '>f8'),                       # flux measured by a fit to the PSF model (erg s^{-1} cm^{-2} Hz^{-1})
('iflux_psf_err', '>f8'),                   # uncertainty for flux.psf
('iflux_psf_flags', 'i1'),                  # set if the flux.psf measurement failed
('ishape_sdss_ixx', '>f4'),                 # Adaptive moments in arcsec 
('ishape_sdss_iyy', '>f4'),                 # Adaptive moments in arcsec   
('ishape_sdss_ixy', '>f4'),                 # Adaptive moments in arcsec 
('ishape_sdss_ixx_var', '>f4'),
('ishape_sdss_iyy_var', '>f4'),
('ishape_sdss_ixy_var', '>f4'),
('ishape_sdss_psf_ixx', '>f4'),             # Adaptive moments of PSF evaluated at object position in arcsec
('ishape_sdss_psf_iyy', '>f4'),             # Adaptive moments of PSF evaluated at object position in arcsec
('ishape_sdss_psf_ixy', '>f4'),             # Adaptive moments of PSF evaluated at object position in arcsec
('tract', '>i4'),                           # tract ID
('icalib_psf_used', 'i1'),                  # Propagated from visits (True if used, False if not)
('merge_peak_g', 'i1'),                     # peak detected in grizy−band
('merge_peak_r', 'i1'),                     # peak detected in grizy−band
('merge_peak_i', 'i1'),                     # peak detected in grizy−band
('merge_peak_z', 'i1'),                     # peak detected in grizy−band
('merge_peak_y', 'i1'),                     # peak detected in grizy−band
('icountinputs', '>i2'),                    # number of grizy−band visits contributing at center
('ideblend_has_stray_flux', 'i1'),          # Avoid spurious detections and those contaminated by blends
('iflags_pixel_bright_object_center', 'i1'),# source center is close to/source footprint includes BRIGHT_OBJECT pixels
('iflags_pixel_bright_object_any', 'i1'),   # source center is close to/source footprint includes BRIGHT_OBJECT pixels
('iblendedness_abs_flux', '>f4'),           # measure of how flux is affected by neighbors
('iflags_negative', 'i1'),                  # set if source was detected as significantly negative 
('ideblend_too_many_peaks', 'i1'),          # Source had too many peaks; only the brightest were included
('ideblend_parent_too_big', 'i1'),          # Parent footprint covered too many pixels
('icentroid_naive_flags', 'i1'),            # set if the centroid.naive measurement did not fully succeed
('iflags_pixel_interpolated_any', 'i1'),    # A pixel flagged as interpolated is close to object center
('iflags_pixel_saturated_any', 'i1'),       # A pixel flagged as saturated is close to object center
('iflags_pixel_cr_any', 'i1'),              # A pixel flagged as a cosmic ray hit is close to object center
('iflags_pixel_suspect_any', 'i1')]))       # source's footprint includes suspect pixels

T_psf = I_xx + I_yy (i.e. psf size))
''';
