'''
This code takes the public star catalogs and converts
them into h5 format that can be ingested in TXPipe.
'''
import h5py 
import numpy as np
from astropy.io import fits

fieldi    = int(sys.argv[1])

fieldlist = ['GAMA09H','GAMA15H',' HECTOMAP', 'VVDS', 'WIDE12H', 'XMM']
field     = fieldlist[fieldi] 

d = fits.open('/global/cfs/cdirs/lsst/groups/WL/users/yomori/scratch/HSC/hsc_stars/%s_stars.fits'%field)

# If 'icalib_psf_used'==True, the star was used to model the psf 
# If 'icalib_psf_used'==False the star is used as a test ("reserved star")
# according to HSC paper the 20% of the stars are reserved
idx_v = np.where(d[1].data['icalib_psf_used']==True)[0]
idx_r = np.where(d[1].data['icalib_psf_used']==False)[0]

print("Total number of stars: %d"%(len(idx_v)+len(idx_t) ) )
print("Number of stars used : %d"%(len(idx_v)))
print("Number of stars reserved : %d"%(len(idx_r)))

print("Computing T,e1,e2 from ixx,iyy,ixy")
ixx_meas  = d[1].data['ishape_sdss_ixx'][idx_r]
iyy_meas  = d[1].data['ishape_sdss_iyy'][idx_r] 
ixy_meas  = d[1].data['ishape_sdss_ixy'][idx_r]
T_meas    = ixx_meas + iyy_meas
e1_meas   = (ixx_meas-iyy_meas)/T_meas
e2_meas   = 2*ixy_meas/T_meas

ixx_model = d[1].data['ishape_sdss_psf_ixx'][idx_r]
iyy_model = d[1].data['ishape_sdss_psf_iyy'][idx_r] 
ixy_model = d[1].data['ishape_sdss_psf_ixy'][idx_r]
T_model   = ixx_model + iyy_model
e1_model  = (ixx_model-iyy_model)/T_model
e2_model  = 2*ixy_model/T_model

f=h5py.File('./star_catalog_hscy1_%s.h5'%field, 'w')
f.create_group("stars")
f['stars/ra']          = d[1].data['ira'][idx_r]
f['stars/dec']         = d[1].data['idec'][idx_r]
f['stars/measured_T']  = d[1].data['T_meas'][idx_r]
f['stars/measured_e1'] = d[1].data['e1_meas'][idx_r]
f['stars/measured_e2'] = d[1].data['e2_meas'][idx_r]
f['stars/model_T']     = d[1].data['T_model'][idx_r]
f['stars/model_e1']    = d[1].data['e1_model'][idx_r]
f['stars/model_e2']    = d[1].data['e2_model'][idx_r]
f['stars/calib_psf_reserved'] = np.ones(idx_r)#d[1].data['ira'][idx_r]
f,close()



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