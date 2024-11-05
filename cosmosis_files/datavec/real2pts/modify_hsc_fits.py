from astropy.io import fits
import numpy as np
import shutil

a=fits.open('hsc_y3_real_withnz.fits')

#xip binning
xipangmin = np.tile((10**(np.log10(3)+0.128*np.arange(14)))[:-1],10)
xipangmax = np.tile((10**(np.log10(3)+0.128*np.arange(14)))[1:] ,10)

#Fudging of xim binning but this still doesnt work unless you set the tolerance off 
ximangmin = (10**(np.log10(13.09548)+0.128*np.arange(11)))[:-1]*0.9999996156845705
ximangmax = (10**(np.log10(13.09548)+0.128*np.arange(11)))[1:]*0.9999996156845705
ximangmin = np.tile(ximangmin,10)
ximangmax = np.tile(ximangmax,10)

# Step 1: Open the FITS file and access the table
filename = "hsc_y3_real_withnz_mod.fits" 
shutil.copy('hsc_y3_real_withnz.fits', filename)

with fits.open(filename, mode="update") as hdul:
    table_hdu = hdul[2]  # Assuming the table is in the first extension
    original_table = table_hdu.data
    original_columns = table_hdu.columns

    hdul[1].header['COVDATA'] = 'T'
    
    # Step 2: Define new columns
    # For example, adding a new column 'new_col1' and 'new_col2' with dummy data
    new_col1 = fits.Column(name='ANGLEMIN', format='D', array=xipangmin)
    new_col2 = fits.Column(name='ANGLEMAX', format='D', array=xipangmax)

    # Step 3: Append the new columns to the existing columns
    new_columns = original_columns + fits.ColDefs([new_col1, new_col2])
    new_table_hdu = fits.BinTableHDU.from_columns(new_columns)
    new_table_hdu.name = "xip"
    # Step 4: Update the HDU list and write back to the file
    hdul[2] = new_table_hdu
    hdul[2].header['QUANT1']='G+R'
    hdul[2].header['QUANT2']='G+R'
    hdul[2].header['KERNEL_1']='nz_source'
    hdul[2].header['KERNEL_2']='nz_source'
    hdul[2].header['WINDOWS']='SAMPLE'
    hdul[2].header['N_ZBIN_1']=4
    hdul[2].header['N_ZBIN_2']=4
    hdul[2].header['N_ANG']   = 13
    hdul[2].header['EXTNAME']   = 'xip'
    hdul[2].header['2PTDATA'] = 'T'
    hdul.flush()  # Save changes to the original FITS file

with fits.open(filename, mode="update") as hdul:
    table_hdu = hdul[3]  # Assuming the table is in the first extension
    original_table = table_hdu.data
    original_columns = table_hdu.columns

    # Step 2: Define new columns
    # For example, adding a new column 'new_col1' and 'new_col2' with dummy data
    new_col1 = fits.Column(name='ANGLEMIN', format='D', array=ximangmin)
    new_col2 = fits.Column(name='ANGLEMAX', format='D', array=ximangmax)

    # Step 3: Append the new columns to the existing columns
    new_columns = original_columns + fits.ColDefs([new_col1, new_col2])
    new_table_hdu = fits.BinTableHDU.from_columns(new_columns)
    new_table_hdu.name = "xim"
    # Step 4: Update the HDU list and write back to the file
    hdul[3] = new_table_hdu
    hdul[3].header['QUANT1']='G+R'
    hdul[3].header['QUANT2']='G+R'
    hdul[3].header['KERNEL_1']='nz_source'
    hdul[3].header['KERNEL_2']='nz_source'
    hdul[3].header['WINDOWS']='SAMPLE'
    hdul[3].header['N_ZBIN_1']=4
    hdul[3].header['N_ZBIN_2']=4
    hdul[3].header['N_ANG']   = 10
    hdul[3].header['EXTNAME']   = 'xim'
    hdul[3].header['2PTDATA'] = 'T'
    
    hdul.flush()  # Save changes to the original FITS file

